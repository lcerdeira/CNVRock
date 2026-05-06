#!/bin/bash
#SBATCH --job-name=kpsc_plasmid_remap
#SBATCH --output=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/plasmid_remap_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Extract unmapped reads from each existing BAM and remap to individual
# per-gene plasmid references to obtain read counts at plasmid resistance genes.
#
# Per-gene references avoid MAPQ=0 from inter-plasmid cross-mapping:
# using a combined reference causes closely-related genes (e.g. blaCTX-M-15
# vs blaCTX-M-14, or any plasmid sharing repeat regions) to compete, assigning
# MAPQ=0 to reads that would otherwise map uniquely to a single gene's plasmid.
#
# Each gene is mapped independently to its own representative plasmid FASTA
# (assets/plasmid_refs/<gene>.fasta, pre-indexed with bwa index).
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
#   # Build per-gene BWA indices once (if not already done):
#   for f in assets/plasmid_refs/*.fasta; do bwa index "$f"; done
#   # Delete old counts if gene list has changed:
#   rm -f data/inputs/plasmid_remap_counts/*.plasmid_counts.tsv
#   N=$(wc -l < assets/kpsc_bam_accessions.txt)
#   sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh

set -euo pipefail

module load samtools/1.20
module load bwa/0.718

REPO_DIR="/home/lshlt19/CNVRock"
BAM_DIR="$REPO_DIR/data/raw/bam"
PLASMID_DIR="$REPO_DIR/assets/plasmid_refs"
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

# ── Validate that per-gene plasmid FASTAs and BWA indices exist ───────────────
for gene in "${GENE_NAMES[@]}"; do
    fa="$PLASMID_DIR/${gene}.fasta"
    if [[ ! -f "${fa}.bwt" ]]; then
        echo "ERROR: BWA index missing for $fa — run: bwa index $fa" >&2
        exit 1
    fi
done

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

# ── Extract unmapped reads as FASTQ (once, reused for all genes) ──────────────
samtools view -b -f 4 "$BAM" | \
    samtools fastq - > "$TMP/unmapped.fq"

N_READS=$(( $(wc -l < "$TMP/unmapped.fq") / 4 ))
echo "  Unmapped reads: $N_READS"

# ── Map to each gene's plasmid individually and count ────────────────────────
HEADER="sample_id"
VALUES="$ACC"
for i in "${!GENE_NAMES[@]}"; do
    gene="${GENE_NAMES[$i]}"
    region="${GENE_REGIONS[$i]}"
    HEADER="${HEADER}\t${gene}"

    if [[ "$N_READS" -gt 0 ]]; then
        gene_fa="$PLASMID_DIR/${gene}.fasta"
        bwa mem -t "$SLURM_CPUS_PER_TASK" "$gene_fa" "$TMP/unmapped.fq" 2>/dev/null | \
            samtools view -b -F 4 -q 10 | \
            samtools sort -o "$TMP/${gene}.bam"
        samtools index "$TMP/${gene}.bam"
        cnt=$(samtools view -c -F 4 "$TMP/${gene}.bam" "$region" 2>/dev/null || echo 0)
        rm -f "$TMP/${gene}.bam" "$TMP/${gene}.bam.bai"
    else
        cnt=0
    fi

    VALUES="${VALUES}\t${cnt}"
    echo "  ${gene}: ${cnt}"
done

# ── Write output ──────────────────────────────────────────────────────────────
echo -e "$HEADER" > "$OUT"
echo -e "$VALUES" >> "$OUT"
echo "Done: $ACC"
