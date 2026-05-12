#!/bin/bash
#SBATCH --job-name=kpsc_exp_rc
#SBATCH --output=/home/lshlt19/CNVRock/logs/exp_rc_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/exp_rc_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --array=1-5000%100   # overridden at submission time

# Extract chromosomal 1 kb bin read counts for Phase 2 expansion BAMs.
#
# Reads manifest:  assets/kpsc_expansion_paths_to_bams.tsv
#   (pre-built: sample_id <TAB> /path/to/bam)
# Outputs counts to: data/raw/readcounts_expansion/{sample_id}.counts.tsv
#
# Usage (run after download_expansion_sra.sh completes):
#   cd ~/CNVRock
#   N=$(wc -l < assets/kpsc_expansion_paths_to_bams.tsv)
#   sbatch --array=1-${N}%100 hpc/collect_expansion_readcounts.sh
#
# Or chain automatically:
#   DL_JOB=$(sbatch --parsable --array=1-${N}%50 hpc/download_expansion_sra.sh)
#   sbatch --dependency=afterok:${DL_JOB} --array=1-${N}%100 hpc/collect_expansion_readcounts.sh

set -euo pipefail

module load gatk/4.6.0.0 java/20.0.1

REPO_DIR="/home/lshlt19/CNVRock"
MANIFEST="$REPO_DIR/assets/kpsc_expansion_paths_to_bams.tsv"
# Use whole-chrom intervals (NC_016845.1 only) — compatible with both the
# basic HS11286.fasta BAMs and the HS11286_extended.fasta BAMs.
INTERVALS="$REPO_DIR/assets/kpsc-whole-chrom.interval_list"
REFERENCE="$REPO_DIR/assets/HS11286_extended.fasta"
OUT_DIR="$REPO_DIR/data/raw/readcounts_expansion"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
LINE=$(sed -n "${LINE_NUM}p" "$MANIFEST")
if [[ -z "$LINE" ]]; then
    echo "No entry for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

SAMPLE_ID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f2)
OUT="$OUT_DIR/${SAMPLE_ID}.counts.tsv"

if [[ -f "$OUT" ]]; then
    echo "$SAMPLE_ID already done — skipping."
    exit 0
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found for $SAMPLE_ID: $BAM" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ID"

SCRATCH_TMP="/tmp/${SAMPLE_ID}_gatk_$$"
mkdir -p "$SCRATCH_TMP"
trap "rm -rf $SCRATCH_TMP" EXIT

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
    --tmp-dir                 "$SCRATCH_TMP" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"
echo "Done: $SAMPLE_ID → $OUT"
