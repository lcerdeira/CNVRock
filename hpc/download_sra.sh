#!/bin/bash
#SBATCH --job-name=kpsc_sra_download
#SBATCH --output=/home/lshlt19/CNVRock/logs/sra_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/sra_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --array=1-5000%50   # overridden at submission time

# Usage:
#   cd ~/CNVRock
#   N=$(wc -l < assets/kpsc_sra_accessions.txt)
#   sbatch --array=1-${N}%50 hpc/download_sra.sh
#
# Requires:
#   assets/kpsc_sra_accessions.txt  — one SRA run accession per line
#   assets/kpsc_run_ftp_map.tsv     — run_accession, ftp_r1, ftp_r2
#   assets/HS11286.fasta + bwa index files

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718
module load samtools/1.20

ACCESSIONS="$REPO_DIR/assets/kpsc_sra_accessions.txt"
FTP_MAP="$REPO_DIR/assets/kpsc_run_ftp_map.tsv"
REFERENCE="$REPO_DIR/assets/HS11286.fasta"
FASTQ_TMP="$REPO_DIR/data/raw/fastq"
BAM_DIR="$REPO_DIR/data/raw/bam"

mkdir -p "$FASTQ_TMP" "$BAM_DIR"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS")
if [[ -z "$ACC" ]]; then
    echo "No accession for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

BAM_OUT="$BAM_DIR/${ACC}.bam"
if [[ -f "$BAM_OUT" && -f "$BAM_OUT.bai" ]]; then
    echo "$ACC already done — skipping."
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC"

R1_URL=$(awk -v acc="$ACC" 'NR>1 && $1==acc {print $2}' "$FTP_MAP")
R2_URL=$(awk -v acc="$ACC" 'NR>1 && $1==acc {print $3}' "$FTP_MAP")

if [[ -z "$R1_URL" || -z "$R2_URL" ]]; then
    echo "ERROR: no FTP URLs found for $ACC in $FTP_MAP"
    exit 1
fi

R1="$FASTQ_TMP/${ACC}_1.fastq.gz"
R2="$FASTQ_TMP/${ACC}_2.fastq.gz"

wget -q "ftp://$R1_URL" -O "$R1"
wget -q "ftp://$R2_URL" -O "$R2"

bwa mem -t "$SLURM_CPUS_PER_TASK" \
    -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
    "$REFERENCE" "$R1" "$R2" \
    | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM_OUT"

samtools index "$BAM_OUT"

rm -f "$R1" "$R2"

printf '%s\t%s\n' "$ACC" "$BAM_OUT" \
    >> "$REPO_DIR/assets/kpsc-paths-to-bams.tsv.partial"

echo "Done: $ACC → $BAM_OUT"
