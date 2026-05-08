#!/bin/bash
#SBATCH --job-name=kpsc_exp_dl
#SBATCH --output=/home/lshlt19/CNVRock/logs/exp_dl_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/exp_dl_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --array=1-5000%50   # overridden at submission time

# Download + align Phase 2 KpSC expansion samples (88,135 new samples).
#
# Differences from download_sra.sh:
#   - reads from kpsc_expansion_sra_accessions.txt
#   - ENA FTP URLs computed inline (no pre-built map needed)
#   - BAMs written to data/raw/bam/  (same pool as existing cohort)
#   - manifest written to assets/kpsc_expansion_paths_to_bams.tsv.partial
#
# Usage:
#   cd ~/CNVRock
#   N=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)
#   sbatch --array=1-${N}%50 hpc/download_expansion_sra.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
export PATH="/home/lshlt19/miniconda3/envs/cnvrock/bin:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718
module load samtools/1.20

ACCESSIONS="$REPO_DIR/assets/kpsc_expansion_sra_accessions.txt"
REFERENCE="$REPO_DIR/assets/HS11286_extended.fasta"
FASTQ_TMP="$REPO_DIR/data/raw/fastq"
BAM_DIR="$REPO_DIR/data/raw/bam"

mkdir -p "$FASTQ_TMP" "$BAM_DIR" "$REPO_DIR/logs"

# BATCH_OFFSET lets us split >5000-task jobs across multiple submissions:
#   sbatch --export=BATCH_OFFSET=5000 --array=1-5000%50 download_expansion_sra.sh
LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
ACC=$(sed -n "${LINE_NUM}p" "$ACCESSIONS")
if [[ -z "$ACC" ]]; then
    echo "No accession for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

BAM_OUT="$BAM_DIR/${ACC}.bam"
if [[ -f "$BAM_OUT" && -f "$BAM_OUT.bai" ]]; then
    echo "$ACC already done — skipping."
    printf '%s\t%s\n' "$ACC" "$BAM_OUT" \
        >> "$REPO_DIR/assets/kpsc_expansion_paths_to_bams.tsv.partial"
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC"

# ── Compute ENA FTP URL inline ─────────────────────────────────────────────
# Pattern: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{acc[:6]}/{sub}/{acc}/{acc}_{1,2}.fastq.gz
# sub = '' for ≤9-char accessions; '00{last_digit}' for longer ones.
PREFIX="${ACC:0:6}"
if [[ ${#ACC} -le 9 ]]; then
    FTP_PATH="ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${ACC}"
else
    SUB="00${ACC: -1}"
    FTP_PATH="ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${SUB}/${ACC}"
fi
R1_URL="ftp://${FTP_PATH}/${ACC}_1.fastq.gz"
R2_URL="ftp://${FTP_PATH}/${ACC}_2.fastq.gz"

R1="$FASTQ_TMP/${ACC}_1.fastq.gz"
R2="$FASTQ_TMP/${ACC}_2.fastq.gz"

# Download with retry
for URL in "$R1_URL $R1" "$R2_URL $R2"; do
    SRC=$(echo "$URL" | cut -d' ' -f1)
    DST=$(echo "$URL" | cut -d' ' -f2)
    ATTEMPT=0
    while [[ $ATTEMPT -lt 3 ]]; do
        if wget -q --timeout=300 --tries=1 "$SRC" -O "$DST"; then break; fi
        ATTEMPT=$(( ATTEMPT + 1 ))
        echo "  wget attempt $ATTEMPT failed for $SRC; retrying in 30 s …"
        sleep 30
    done
    if [[ ! -f "$DST" || ! -s "$DST" ]]; then
        echo "ERROR: failed to download $SRC" >&2
        exit 1
    fi
done

# ── Align to extended reference ────────────────────────────────────────────
bwa mem -t "$SLURM_CPUS_PER_TASK" \
    -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
    "$REFERENCE" "$R1" "$R2" \
    | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM_OUT"

samtools index "$BAM_OUT"

rm -f "$R1" "$R2"

printf '%s\t%s\n' "$ACC" "$BAM_OUT" \
    >> "$REPO_DIR/assets/kpsc_expansion_paths_to_bams.tsv.partial"

echo "Done: $ACC → $BAM_OUT"
