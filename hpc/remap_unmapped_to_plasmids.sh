#!/bin/bash
#SBATCH --job-name=kpsc_plasmid_remap
#SBATCH --output=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Phase B: extract unmapped reads from each existing BAM and remap to
# HS11286_extended.fasta to obtain read counts at blaCTX-M-15 and blaNDM-1.
#
# The existing BAMs were mapped to HS11286 (chromosome + native plasmids).
# MK552109.1 (blaCTX-M-15) and MZ606384.2 (blaNDM-1) were not in that
# reference, so reads from those plasmids are in the unmapped fraction.
#
# Output per sample: <accession>.plasmid_counts.tsv (one line)
#   <acc>  <ctxm_count>  <ndm_count>
#
# Usage:
#   cd ~/CNVRock
#   N=$(wc -l < assets/kpsc_bam_accessions.txt)
#   sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh

set -euo pipefail

module load samtools/1.20
module load bwa/0.718

REPO_DIR="/home/lshlt19/CNVRock"
BAM_DIR="$REPO_DIR/data/raw/bam"
REF="$REPO_DIR/assets/HS11286_extended.fasta"
ACCS="$REPO_DIR/assets/kpsc_bam_accessions.txt"
OUT_DIR="$REPO_DIR/data/inputs/plasmid_remap_counts"

# Gene loci (±500 bp padding, matching existing blaKPC-2 bin convention)
# blaCTX-M-15: MK552109.1:119392-120264 → padded 118892-120764
# blaNDM-1:    MZ606384.2:90937-91746   → padded 90437-92246
CTXM_REGION="MK552109.1:118892-120764"
NDM_REGION="MZ606384.2:90437-92246"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCS")
if [[ -z "$ACC" ]]; then
    echo "No accession for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

BAM="$BAM_DIR/${ACC}.bam"
OUT="$OUT_DIR/${ACC}.plasmid_counts.tsv"

if [[ -f "$OUT" ]]; then
    echo "$ACC already done — skipping."
    exit 0
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found: $BAM" >&2
    exit 1
fi

TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

echo "Task $SLURM_ARRAY_TASK_ID: $ACC"

# Step 1: Extract unmapped reads as FASTQ
samtools view -b -f 4 "$BAM" | \
    samtools fastq - > "$TMP/unmapped.fq"

N_READS=$(( $(wc -l < "$TMP/unmapped.fq") / 4 ))
echo "  Unmapped reads: $N_READS"

if [[ "$N_READS" -eq 0 ]]; then
    echo -e "${ACC}\t0\t0" > "$OUT"
    echo "  No unmapped reads — writing zeros."
    exit 0
fi

# Step 2: Remap to extended reference (new plasmid contigs will attract their reads)
bwa mem -t "$SLURM_CPUS_PER_TASK" "$REF" "$TMP/unmapped.fq" 2>/dev/null | \
    samtools view -b -F 4 -q 10 | \
    samtools sort -o "$TMP/plasmid.bam"
samtools index "$TMP/plasmid.bam"

# Step 3: Count reads at each gene locus
ctxm_count=$(samtools view -c -F 4 "$TMP/plasmid.bam" "$CTXM_REGION" 2>/dev/null || echo 0)
ndm_count=$(samtools  view -c -F 4 "$TMP/plasmid.bam" "$NDM_REGION"  2>/dev/null || echo 0)

echo -e "${ACC}\t${ctxm_count}\t${ndm_count}" > "$OUT"
echo "  CTX-M=$ctxm_count  NDM=$ndm_count"
