#!/bin/bash
#SBATCH --job-name=kpsc_kpc2_rc
#SBATCH --output=/home/lshlt19/CNVRock/logs/kpc2_rc_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/kpc2_rc_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

# Collect blaKPC-2 read counts from original BAMs using GATK CollectReadCounts.
# Uses original HS11286.fasta (7 contigs matching BAM headers) and the
# kpc2-only interval list so GATK doesn't fail on contig mismatches.
#
# Usage:
#   N=$(wc -l < assets/kpsc-paths-to-bams.tsv)
#   sbatch --array=1-${N}%100 hpc/collect_kpc2_readcounts.sh

set -euo pipefail

module load gatk/4.6.0.0 java/20.0.1

REPO_DIR="/home/lshlt19/CNVRock"
MANIFEST="$REPO_DIR/assets/kpsc-paths-to-bams.tsv"
INTERVALS="$REPO_DIR/assets/kpsc-kpc2-only.interval_list"
REFERENCE="$REPO_DIR/assets/HS11286.fasta"
OUT_DIR="$REPO_DIR/models/experiments/21/plasmid_readcounts"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
if [[ -z "$LINE" ]]; then
    echo "No entry for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

SAMPLE_ID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f2)
OUT="$OUT_DIR/${SAMPLE_ID}.plasmid_counts.tsv"

if [[ -f "$OUT" ]]; then
    echo "$SAMPLE_ID already done — skipping."
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ID"

gatk CollectReadCounts \
    --java-options            "-Xmx3g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality 40 \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"
echo "Done: $SAMPLE_ID → $OUT"
