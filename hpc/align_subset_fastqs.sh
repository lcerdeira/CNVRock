#!/bin/bash
#SBATCH --job-name=kpsc_subset_align
#SBATCH --output=/home/lshlt19/CNVRock/logs/align_subset_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/align_subset_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Align pre-downloaded FASTQs (uploaded from local Aspera dump) to the
# HS11286_extended reference and run GATK CollectReadCounts at 1 kb resolution
# across the chromosome + plasmids combined.
#
# Why a separate script: download_expansion_sra.sh fetches from EBI which is
# rate-limited on the HPC. For the 5,000-sample stratified subset we download
# locally via Aspera and rsync to ~/CNVRock/data/raw/fastq_subset/ instead.
#
# Inputs:
#   data/raw/fastq_subset/{ACC}_1.fastq.gz
#   data/raw/fastq_subset/{ACC}_2.fastq.gz  (paired only)
#   assets/kpsc_expansion_subset_5k.tsv     (cols: accession, layout, ...)
#   assets/HS11286_extended.fasta           (must have .fai, .dict, .bwt etc)
#   assets/HS11286_extended_1kb.interval_list   (1 kb bins across all contigs)
#
# Outputs:
#   data/raw/readcounts_subset/{ACC}.counts.tsv   (~7,365 rows + header)
#   FASTQs are deleted after successful count to save space (override KEEP_FASTQ=1).
#
# Submit:
#   N=$(( $(wc -l < assets/kpsc_expansion_subset_5k.tsv) - 1 ))
#   sbatch --array=1-${N}%50 hpc/align_subset_fastqs.sh

set -euo pipefail

REPO_DIR=/home/lshlt19/CNVRock
export PATH="/home/lshlt19/miniconda3/envs/cnvrock/bin:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718 samtools/1.20 gatk/4.6.0.0 java/20.0.1

MANIFEST="$REPO_DIR/assets/kpsc_expansion_subset_5k.tsv"
FASTQ_DIR="$REPO_DIR/data/raw/fastq_subset"
REFERENCE="$REPO_DIR/assets/HS11286_extended.fasta"
INTERVALS="$REPO_DIR/assets/HS11286_extended_1kb.interval_list"
OUT_DIR="$REPO_DIR/data/raw/readcounts_subset"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID + 1 ))  # +1 to skip header
ROW=$(sed -n "${LINE_NUM}p" "$MANIFEST")
[[ -z "$ROW" ]] && { echo "No row for task $SLURM_ARRAY_TASK_ID — skipping."; exit 0; }

ACC=$(echo    "$ROW" | cut -f1)
LAYOUT=$(echo "$ROW" | cut -f2)
OUT="$OUT_DIR/${ACC}.counts.tsv"

[[ -f "$OUT" ]] && { echo "$ACC already done — skipping."; exit 0; }

R1="$FASTQ_DIR/${ACC}_1.fastq.gz"
R2="$FASTQ_DIR/${ACC}_2.fastq.gz"
if [[ ! -s "$R1" ]]; then
    echo "ERROR: R1 missing for $ACC: $R1" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC  (layout=$LAYOUT)  on $(hostname)"

SCRATCH="/tmp/${ACC}_align_$$"
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

BAM="$SCRATCH/${ACC}.bam"

# ── BWA align ─────────────────────────────────────────────────────────────────
if [[ "$LAYOUT" == "PAIRED" && -s "$R2" ]]; then
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" "$R2" \
      | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
else
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" \
      | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
fi
samtools index "$BAM"

# ── GATK CollectReadCounts ────────────────────────────────────────────────────
gatk CollectReadCounts \
    --java-options            "-Xmx12g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality 40 \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --tmp-dir                 "$SCRATCH" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"

# ── Free FASTQ space ──────────────────────────────────────────────────────────
if [[ "${KEEP_FASTQ:-0}" != "1" ]]; then
    rm -f "$R1" "$R2"
fi

echo "Done: $ACC → $OUT"
