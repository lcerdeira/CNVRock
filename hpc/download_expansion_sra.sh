#!/bin/bash
#SBATCH --job-name=kpsc_exp_dl
#SBATCH --output=/home/lshlt19/CNVRock/logs/exp_dl_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/exp_dl_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --array=1-5000%50   # overridden at submission time

# Download + align + GATK read counts for Phase 2 KpSC expansion samples.
#
# Pipeline per sample (no intermediate files kept):
#   1. Download R1 [+ R2] FASTQ.gz via HTTPS → local /tmp scratch
#   2. BWA mem align → sorted BAM in /tmp scratch
#   3. GATK CollectReadCounts → counts TSV in data/raw/readcounts_expansion/
#   4. Delete FASTQs and BAM (only count file is kept)
#
# This avoids storing ~34 TB of FASTQs and ~44 TB of BAMs on the HPC.
# Total permanent storage: ~88 GB (count files only).
#
# Prerequisites:
#   python3 data/setup/build_ena_url_manifest.py  (builds assets/ena_url_manifest.tsv)
#
# Usage (wave-based):
#   nohup bash hpc/submit_bam_wave.sh > ~/bam_wave.log 2>&1 &

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
export PATH="/home/lshlt19/miniconda3/envs/cnvrock/bin:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718
module load samtools/1.20

MANIFEST="$REPO_DIR/assets/ena_url_manifest.tsv"
REFERENCE="$REPO_DIR/assets/HS11286_extended.fasta"
INTERVALS="$REPO_DIR/assets/HS11286_extended_intervals.interval_list"
COUNTS_DIR="$REPO_DIR/data/raw/readcounts_expansion"

mkdir -p "$COUNTS_DIR" "$REPO_DIR/logs"

# BATCH_OFFSET allows wave-based submission
LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID + 1 ))  # +1 to skip header

ROW=$(sed -n "${LINE_NUM}p" "$MANIFEST")
if [[ -z "$ROW" ]]; then
    echo "No row for task $SLURM_ARRAY_TASK_ID (line $LINE_NUM) — skipping."
    exit 0
fi

ACC=$(echo    "$ROW" | cut -f1)
LAYOUT=$(echo "$ROW" | cut -f2)
R1_URL=$(echo "$ROW" | cut -f3)
R2_URL=$(echo "$ROW" | cut -f4)

COUNTS_OUT="$COUNTS_DIR/${ACC}.counts.tsv"
if [[ -f "$COUNTS_OUT" ]]; then
    echo "$ACC already done — skipping."
    exit 0
fi

if [[ -z "$R1_URL" ]]; then
    echo "ERROR: no URL for $ACC in manifest — skipping." >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC  (layout=$LAYOUT)"

# ── Use /tmp scratch for all intermediate files ───────────────────────────────
SCRATCH="/tmp/${ACC}_$$"
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

R1="$SCRATCH/${ACC}_1.fastq.gz"
R2="$SCRATCH/${ACC}_2.fastq.gz"
BAM="$SCRATCH/${ACC}.bam"

# ── Download FASTQs via HTTPS ─────────────────────────────────────────────────
download_file() {
    local URL="$1" DST="$2"
    local ATTEMPT=0
    while [[ $ATTEMPT -lt 3 ]]; do
        if wget -q --timeout=300 --tries=1 "$URL" -O "$DST"; then
            [[ -s "$DST" ]] && return 0
        fi
        ATTEMPT=$(( ATTEMPT + 1 ))
        echo "  wget attempt $ATTEMPT failed for $URL; retrying in 120 s …"
        rm -f "$DST"
        sleep 120
    done
    echo "ERROR: failed to download $URL" >&2
    return 1
}

download_file "$R1_URL" "$R1"
if [[ "$LAYOUT" == "PAIRED" && -n "$R2_URL" ]]; then
    download_file "$R2_URL" "$R2"
fi

# ── Align to extended reference ───────────────────────────────────────────────
if [[ "$LAYOUT" == "PAIRED" && -f "$R2" && -s "$R2" ]]; then
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" "$R2" \
        | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
else
    echo "  Aligning single-end reads for $ACC"
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" \
        | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
fi

samtools index "$BAM"
rm -f "$R1" "$R2"   # free FASTQ space immediately

# ── GATK CollectReadCounts ────────────────────────────────────────────────────
gatk CollectReadCounts \
    -I "$BAM" \
    -L "$INTERVALS" \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O "${COUNTS_OUT}.tmp" \
    --tmp-dir "$SCRATCH"

mv "${COUNTS_OUT}.tmp" "$COUNTS_OUT"
# BAM and index removed by trap on EXIT

echo "Done: $ACC → $COUNTS_OUT"
