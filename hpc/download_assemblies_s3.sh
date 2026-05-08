#!/bin/bash
#SBATCH --job-name=kpsc_asm_s3
#SBATCH --output=/home/lshlt19/CNVRock/logs/asm_s3_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/asm_s3_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --array=1-5000%100   # overridden at submission time

# Download AllTheBacteria KpSC assemblies (FASTA) from S3 for Kleborate GT.
#
# Input:  assets/kpsc_expansion_assembly_urls.tsv
#           columns: sample_accession, run_accession, aws_url
# Output: data/assemblies/{sample_accession}.fa.gz
#
# Usage:
#   cd ~/CNVRock
#   N=$(( $(wc -l < assets/kpsc_expansion_assembly_urls.tsv) - 1 ))
#   sbatch --array=1-${N}%100 hpc/download_assemblies_s3.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
URL_TSV="$REPO_DIR/assets/kpsc_expansion_assembly_urls.tsv"
ASM_DIR="$REPO_DIR/data/assemblies"

mkdir -p "$ASM_DIR" "$REPO_DIR/logs"

# Skip header row (line 1 is header, task IDs start at 1 → add 1)
LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
ROW=$(sed -n "${LINE}p" "$URL_TSV")

if [[ -z "$ROW" ]]; then
    echo "No row for task $SLURM_ARRAY_TASK_ID (line $LINE) — skipping."
    exit 0
fi

SAMPLE_ACC=$(echo "$ROW" | awk -F'\t' '{print $1}')
AWS_URL=$(echo "$ROW"    | awk -F'\t' '{print $3}')

OUT_FILE="$ASM_DIR/${SAMPLE_ACC}.fa.gz"

if [[ -f "$OUT_FILE" ]]; then
    echo "$SAMPLE_ACC already downloaded — skipping."
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ACC → $AWS_URL"

# Download with retry (3 attempts, 30 s timeout per attempt)
ATTEMPT=0
while [[ $ATTEMPT -lt 3 ]]; do
    if wget -q --timeout=120 --tries=1 "$AWS_URL" -O "$OUT_FILE.tmp"; then
        mv "$OUT_FILE.tmp" "$OUT_FILE"
        echo "Done: $SAMPLE_ACC"
        exit 0
    fi
    ATTEMPT=$(( ATTEMPT + 1 ))
    echo "Attempt $ATTEMPT failed for $SAMPLE_ACC; retrying in 10 s …"
    sleep 10
done

rm -f "$OUT_FILE.tmp"
echo "ERROR: all attempts failed for $SAMPLE_ACC ($AWS_URL)" >&2
exit 1
