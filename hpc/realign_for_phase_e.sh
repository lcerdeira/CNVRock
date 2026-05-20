#!/bin/bash
#SBATCH --job-name=phaseE_realign
#SBATCH --output=/home/lshlt19/CNVRock/logs/phaseE_align_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/phaseE_align_%A_%a.err
#SBATCH --time=00:45:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Phase E: re-align each existing FASTQ (or BAM-extract) to the updated
# HS11286_extended.fasta (now carrying 8 new plasmid contigs) and emit
# fresh MQ=20 chromosome counts + MQ=0 plasmid counts including the new
# contigs.
#
# Strategy: extract reads from the persistent BAM via samtools fastq,
# realign with BWA-MEM to the new extended reference, sort/index, then
# run two GATK CollectReadCounts passes on the new interval list.
#
# Idempotent: skips samples whose Phase E count files already exist.
set -euo pipefail

REPO=/home/lshlt19/CNVRock
cd "$REPO"
export PATH="$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"
module load bwa samtools gatk/4.6.0.0 java/20.0.1 2>/dev/null || true

REF="$REPO/assets/HS11286_extended.fasta"
INT="$REPO/assets/HS11286_extended_1kb.interval_list"
BAM_DIR="$REPO/data/raw/bam_subset"
OUT_MQ20="$REPO/data/raw/readcounts_phase_e_mq20"
OUT_MQ0="$REPO/data/raw/readcounts_phase_e_mq0"
mkdir -p "$OUT_MQ20" "$OUT_MQ0"

BAM_FILES=( $(ls "$BAM_DIR"/*.bam) )
BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
[[ -z "$BAM" ]] && { echo "no BAM for task $SLURM_ARRAY_TASK_ID"; exit 0; }
ACC=$(basename "$BAM" .bam)

OUT20="$OUT_MQ20/${ACC}.counts.tsv"
OUT00="$OUT_MQ0/${ACC}.counts.tsv"
[[ -s "$OUT20" && -s "$OUT00" ]] && { echo "$ACC done"; exit 0; }

SCRATCH=/tmp/${ACC}_phaseE_$$
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

echo "Realigning $ACC…"
samtools fastq -@ 2 -1 "$SCRATCH/R1.fq.gz" -2 "$SCRATCH/R2.fq.gz" \
               -s /dev/null -0 /dev/null "$BAM"

bwa mem -t 4 "$REF" "$SCRATCH/R1.fq.gz" "$SCRATCH/R2.fq.gz" \
  | samtools sort -@ 2 -o "$SCRATCH/aln.bam" -
samtools index "$SCRATCH/aln.bam"

for MQ in 20 0; do
    OUT=$([[ $MQ == 20 ]] && echo "$OUT20" || echo "$OUT00")
    [[ -s "$OUT" ]] && continue
    gatk CollectReadCounts \
        --java-options "-Xmx10g" \
        --reference "$REF" --intervals "$INT" \
        --input "$SCRATCH/aln.bam" --format TSV \
        --read-filter MappingQualityReadFilter --minimum-mapping-quality $MQ \
        --interval-merging-rule OVERLAPPING_ONLY \
        --tmp-dir "$SCRATCH" \
        --output "${OUT}.tmp"
    mv "${OUT}.tmp" "$OUT"
done
echo "Done $ACC"
