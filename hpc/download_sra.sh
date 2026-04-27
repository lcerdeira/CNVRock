#!/bin/bash
#SBATCH --job-name=kpsc_sra_download
#SBATCH --output=logs/sra_download_%A_%a.out
#SBATCH --error=logs/sra_download_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
# TODO: replace 'cpu' with your LSHTM partition name (check with: sinfo)
#SBATCH --partition=cpu
#SBATCH --array=1-5000%50   # adjust upper bound to number of accessions; 50 jobs at a time

# Usage:
#   mkdir -p logs fastq bam
#   sbatch hpc/download_sra.sh
#
# Requires:
#   assets/kpsc_sra_accessions.txt    — one SRA accession per line
#   assets/HS11286.fasta (+ .bwt etc) — bwa index (run bwa index HS11286.fasta first)
#   modules: sra-tools, bwa, samtools (adjust module load lines below)

set -euo pipefail

# TODO: adjust module names to match LSHTM HPC module system
# Run 'module avail' or 'module spider sra' to find the right names
module load sra-tools 2>/dev/null || true
module load bwa       2>/dev/null || true
module load samtools  2>/dev/null || true

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ACCESSIONS="$REPO_DIR/assets/kpsc_sra_accessions.txt"
REFERENCE="$REPO_DIR/assets/HS11286.fasta"
FASTQ_DIR="$REPO_DIR/data/raw/fastq"
BAM_DIR="$REPO_DIR/data/raw/bam"

mkdir -p "$FASTQ_DIR" "$BAM_DIR"

# Pick this job's accession from the list (1-indexed)
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS")
if [[ -z "$ACC" ]]; then
    echo "No accession for task $SLURM_ARRAY_TASK_ID — skipping."
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: processing $ACC"

BAM_OUT="$BAM_DIR/${ACC}.bam"

# Skip if already done
if [[ -f "$BAM_OUT" && -f "$BAM_OUT.bai" ]]; then
    echo "$ACC already mapped — skipping."
    exit 0
fi

# Download raw reads
fasterq-dump "$ACC" \
    --outdir "$FASTQ_DIR" \
    --split-files \
    --threads "$SLURM_CPUS_PER_TASK" \
    --temp "$FASTQ_DIR/tmp_${ACC}"

R1="$FASTQ_DIR/${ACC}_1.fastq"
R2="$FASTQ_DIR/${ACC}_2.fastq"

# Map to HS11286 and produce sorted, indexed BAM
bwa mem -t "$SLURM_CPUS_PER_TASK" -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
    "$REFERENCE" "$R1" "$R2" \
    | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM_OUT"

samtools index "$BAM_OUT"

# Remove raw fastqs to save space (reads are in the BAM)
rm -f "$R1" "$R2"

echo "$ACC done → $BAM_OUT"

# Append to BAM manifest (append-safe with newline)
printf '%s\t%s\n' "$ACC" "$BAM_OUT" >> "$REPO_DIR/assets/kpsc-paths-to-bams.tsv.partial"
