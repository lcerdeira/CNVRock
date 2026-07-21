#!/bin/bash
#SBATCH --job-name=gatk_mq20
#SBATCH --output=/home/lshlt19/CNVRock/logs/mq20_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/mq20_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal

# GATK CollectReadCounts at MQ=20, reading from existing BAMs.
#
# Used by the hybrid pipeline:
#   - chromosome store    ← MQ=20 (this script's output, cleaner signal)
#   - plasmid-family store ← MQ=0  (aspera_subset_pipeline.sh output, recovers
#                                    multi-mapping allele variants)
#
# Input list: BAMS_LIST env var (default = list every {ACC}.bam in data/raw/bam_subset/)
# Output:     data/raw/readcounts_subset_mq20/{ACC}.counts.tsv
#
# Idempotent: skips samples whose mq20 count file already exists.

set -euo pipefail

REPO_DIR=/home/lshlt19/CNVRock
export PATH="$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"
module load gatk/4.6.0.0 java/20.0.1

# Env-var overridable (KpSC defaults preserved when unset), so other organisms
# can point the recount at their own reference / BAMs without editing this file.
REFERENCE="${REFERENCE:-$REPO_DIR/assets/HS11286_extended.fasta}"
INTERVALS="${INTERVALS:-$REPO_DIR/assets/HS11286_extended_1kb.interval_list}"
BAM_DIR="${BAM_DIR:-$REPO_DIR/data/raw/bam_subset}"
COUNTS_DIR="${COUNTS_DIR:-$REPO_DIR/data/raw/readcounts_subset_mq20}"
BAMS_LIST="${BAMS_LIST:-$REPO_DIR/assets/_bams_list_for_mq20.txt}"

mkdir -p "$COUNTS_DIR" "$REPO_DIR/logs"

LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
BAM=$(sed -n "${LINE_NUM}p" "$BAMS_LIST")
[[ -z "$BAM" ]] && { echo "No BAM for task $SLURM_ARRAY_TASK_ID"; exit 0; }

ACC=$(basename "$BAM" .bam)
OUT="$COUNTS_DIR/${ACC}.counts.tsv"

if [[ -f "$OUT" ]]; then
    echo "$ACC already done"; exit 0
fi

if [[ ! -s "$BAM" || ! -s "${BAM}.bai" ]]; then
    echo "ERROR: BAM/index missing for $ACC" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC on $(hostname)"
SCRATCH=/tmp/${ACC}_mq20_$$
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

gatk CollectReadCounts \
    --java-options            "-Xmx8g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality 20 \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --tmp-dir                 "$SCRATCH" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"
echo "Done: $ACC → $OUT"
