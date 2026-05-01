#!/bin/bash
#SBATCH --job-name=kpsc_plasmid_rc
#SBATCH --output=/home/lshlt19/CNVRock/logs/plasmid_rc_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/plasmid_rc_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --array=1-1000%100   # overridden at submission time

# Collect read counts over plasmid gene intervals from existing BAMs.
#
# Phase A: blaKPC-2 works on BAMs already mapped to HS11286.fasta (which includes
#          NC_016846.1). No re-mapping needed.
# Phase B: blaCTX-M-15 and blaNDM-1 will only have reads if BAMs were mapped to
#          HS11286_extended.fasta (new BAMs going forward).
#
# Usage:
#   cd ~/CNVRock
#   N=$(wc -l < assets/kpsc-paths-to-bams.tsv)
#   sbatch --array=1-${N}%100 hpc/collect_plasmid_readcounts.sh

set -euo pipefail

module load gatk/4.6.0.0 java/20.0.1

REPO_DIR="/home/lshlt19/CNVRock"
MANIFEST="$REPO_DIR/assets/kpsc-paths-to-bams.tsv"
INTERVALS="$REPO_DIR/assets/kpsc-plasmid-genes.interval_list"

# Use extended reference if available, fall back to original for blaKPC-2-only run
EXTENDED="$REPO_DIR/assets/HS11286_extended.fasta"
ORIGINAL="$REPO_DIR/assets/HS11286.fasta"
if [[ -f "$EXTENDED" ]]; then
    REFERENCE="$EXTENDED"
else
    echo "WARNING: HS11286_extended.fasta not found; using original (blaKPC-2 only)."
    REFERENCE="$ORIGINAL"
fi

OUT_DIR="$REPO_DIR/models/experiments/21/plasmid_readcounts"
mkdir -p "$OUT_DIR"

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
