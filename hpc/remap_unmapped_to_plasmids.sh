#!/bin/bash
#SBATCH --job-name=kpsc_plasmid_remap
#SBATCH --output=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Extract unmapped reads from each existing BAM and remap to
# HS11286_plasmids_only.fasta (plasmid contigs only, no chromosome) to obtain
# read counts at plasmid resistance genes.
#
# Using plasmid-only reference avoids MAPQ=0 from chromosomal paralogues
# (e.g. blaSHV competing with blaCTX-M, aac6-Ib competing with aac6-Ib-cr).
#
# Which genes are counted is driven by assets/plasmid_refs/plasmid_gene_coords.tsv.
# blaKPC-2 is excluded — it is already in the original BAMs (NC_016846.1).
#
# Output per sample: <accession>.plasmid_counts.tsv
#   Header: sample_id <tab> gene1 <tab> gene2 ...
#   Data:   <acc>     <tab> count <tab> count ...
#
# Usage:
#   cd ~/CNVRock
#   # Delete old counts if gene list has changed:
#   rm -f data/inputs/plasmid_remap_counts/*.plasmid_counts.tsv
#   N=$(wc -l < assets/kpsc_bam_accessions.txt)
#   sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh

set -euo pipefail

module load samtools/1.20
module load bwa/0.718

REPO_DIR="/home/lshlt19/CNVRock"
BAM_DIR="$REPO_DIR/data/raw/bam"
# Plasmid-only reference excludes the HS11286 chromosome and native plasmid (NC_016845.1,
# NC_016846.1) so that reads cross-mapping to chromosomal paralogues (e.g. blaSHV for
# CTX-M, aac6-Ib for aac6-Ib-cr) are not assigned MAPQ=0 and then filtered out.
REF="$REPO_DIR/assets/HS11286_plasmids_only.fasta"
COORDS="$REPO_DIR/assets/plasmid_refs/plasmid_gene_coords.tsv"
ACCS="$REPO_DIR/assets/kpsc_bam_accessions.txt"
OUT_DIR="$REPO_DIR/data/inputs/plasmid_remap_counts"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

# ── Read gene loci from coords TSV (skip header and blaKPC-2) ────────────────
GENE_NAMES=()
GENE_REGIONS=()
while IFS=$'\t' read -r gene contig start end rest; do
    [[ "$gene" == "gene" ]] && continue        # header row
    [[ "$gene" == "blaKPC-2" ]] && continue   # already in original BAMs
    GENE_NAMES+=("$gene")
    GENE_REGIONS+=("${contig}:${start}-${end}")
done < "$COORDS"

if [[ "${#GENE_NAMES[@]}" -eq 0 ]]; then
    echo "ERROR: no non-KPC genes found in $COORDS" >&2
    exit 1
fi

echo "Genes to count: ${GENE_NAMES[*]}"

# ── Get sample accession for this array task ──────────────────────────────────
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

# ── Extract unmapped reads as FASTQ ──────────────────────────────────────────
samtools view -b -f 4 "$BAM" | \
    samtools fastq - > "$TMP/unmapped.fq"

N_READS=$(( $(wc -l < "$TMP/unmapped.fq") / 4 ))
echo "  Unmapped reads: $N_READS"

# ── Remap to extended reference ───────────────────────────────────────────────
if [[ "$N_READS" -gt 0 ]]; then
    bwa mem -t "$SLURM_CPUS_PER_TASK" "$REF" "$TMP/unmapped.fq" 2>/dev/null | \
        samtools view -b -F 4 -q 10 | \
        samtools sort -o "$TMP/plasmid.bam"
    samtools index "$TMP/plasmid.bam"
    HAVE_BAM=1
else
    HAVE_BAM=0
fi

# ── Count reads at each gene locus ────────────────────────────────────────────
HEADER="sample_id"
VALUES="$ACC"
for i in "${!GENE_NAMES[@]}"; do
    gene="${GENE_NAMES[$i]}"
    region="${GENE_REGIONS[$i]}"
    HEADER="${HEADER}\t${gene}"
    if [[ "$HAVE_BAM" -eq 1 ]]; then
        cnt=$(samtools view -c -F 4 "$TMP/plasmid.bam" "$region" 2>/dev/null || echo 0)
    else
        cnt=0
    fi
    VALUES="${VALUES}\t${cnt}"
    echo "  ${gene}: ${cnt}"
done

# Write header + data (header helps merge script identify columns)
echo -e "$HEADER" > "$OUT"
echo -e "$VALUES" >> "$OUT"
echo "Done: $ACC"
