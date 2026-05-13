#!/bin/bash
#SBATCH --job-name=kpsc_exp_rerun
#SBATCH --output=/home/lshlt19/CNVRock/logs/exp_rerun_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/exp_rerun_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

# Rerun GATK CollectReadCounts on BAMs that are missing their .bai index.
#
# Used during the expansion-cohort GATK pass when 182/8,159 BAMs failed with:
#   A USER ERROR has occurred: Traversal by intervals was requested but some
#   input files are not indexed.
#
# This script:
#   1. Loads samtools (which the main collect_expansion_readcounts.sh did not)
#   2. Runs `samtools index` if the .bai is missing
#   3. Runs GATK CollectReadCounts with the same flags as the main script
#
# Build the rerun manifest first (list of samples still missing count files):
#   > assets/kpsc_rerun_paths_to_bams.tsv
#   for bam in data/raw/bam/*.bam; do
#     sample=$(basename "$bam" .bam)
#     out="data/raw/readcounts_expansion/${sample}.counts.tsv"
#     [[ -f "$out" ]] && continue
#     echo -e "${sample}\t/home/lshlt19/CNVRock/${bam}" \
#       >> assets/kpsc_rerun_paths_to_bams.tsv
#   done
#
# Submit:
#   N=$(wc -l < assets/kpsc_rerun_paths_to_bams.tsv)
#   sbatch --array=1-${N}%50 hpc/rerun_unindexed_bams.sh

set -euo pipefail
module load gatk/4.6.0.0 java/20.0.1 samtools/1.20

MANIFEST=/home/lshlt19/CNVRock/assets/kpsc_rerun_paths_to_bams.tsv
INTERVALS=/home/lshlt19/CNVRock/assets/kpsc-whole-chrom.interval_list
REFERENCE=/home/lshlt19/CNVRock/assets/HS11286.fasta
OUT_DIR=/home/lshlt19/CNVRock/data/raw/readcounts_expansion

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[[ -z "$LINE" ]] && exit 0
SAMPLE_ID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f2)
OUT="$OUT_DIR/${SAMPLE_ID}.counts.tsv"
[[ -f "$OUT" ]] && echo "$SAMPLE_ID already done" && exit 0
[[ ! -f "$BAM" ]] && echo "ERROR: BAM missing: $BAM" >&2 && exit 1

# Index BAM if missing
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Indexing $BAM ..."
    samtools index "$BAM"
fi

SCRATCH=/tmp/${SAMPLE_ID}_rerun_$$
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

gatk CollectReadCounts \
    --java-options            "-Xmx5g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality 40 \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --disable-sequence-dictionary-validation \
    --tmp-dir                 "$SCRATCH" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"
echo "Done: $SAMPLE_ID"
