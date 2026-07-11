#!/bin/bash
#SBATCH --job-name=lr_ont_depth
#SBATCH --output=/home/lshlt19/CNVRock/logs/lr_ont_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/lr_ont_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#
# Long-read (ONT) DEPTH validation of CNVRock chromosomal blaSHV CRR.
#
# Unlike assembly-based validation (which collapses tandem arrays), long-read
# DEPTH at the blaSHV locus captures tandem amplification the same way short
# reads do — every copy maps to the same reference position — giving an
# INDEPENDENT copy-ratio to correlate against CNVRock's short-read CRR.
#
# Per array task: download the ONT run's fastq, align to HS11286 with
# minimap2 (map-ont), compute mean depth at blaSHV and mean chromosomal
# depth, and write CRR = blaSHV_depth / chrom_depth. Cleans up scratch.
#
# Launch:
#   N=$(($(wc -l < data/results/longread_validation_ont_manifest.tsv)-1))
#   sbatch --array=1-${N}%20 hpc/longread_ont_depth.sh

set -uo pipefail
REPO=/home/lshlt19/CNVRock
BIN=/home/lshlt19/miniforge3/envs/nexus_env/bin        # minimap2 + samtools
export PATH="$BIN:$PATH"

MANIFEST=$REPO/data/results/longread_validation_ont_manifest.tsv
REF=$REPO/assets/HS11286_extended.fasta
OUTDIR=$REPO/data/results/longread_depth
mkdir -p "$OUTDIR" "$REPO/logs"
SCRATCH=/home/lshlt19/scratch/lr_ont/${SLURM_ARRAY_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# blaSHV locus on the HS11286 main chromosome
CHR=NC_016845.1
SHV_START=2549403
SHV_END=2550263

# Manifest columns: sample_accession, run_accession, fastq_ftp, base_count, read_count
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
SAMPLE=$(echo "$LINE" | cut -f1)
RUN=$(echo "$LINE"    | cut -f2)
FTP=$(echo "$LINE"    | cut -f3)
OUT=$OUTDIR/${RUN}.tsv
[ -s "$OUT" ] && { echo "already done: $RUN"; exit 0; }
[ -z "$FTP" ] && { echo "no fastq_ftp for $RUN"; exit 0; }

# ENA fastq_ftp may list multiple files (semicolon); take the first (ONT single-end)
URL="https://$(echo "$FTP" | cut -d';' -f1)"
FQ=$SCRATCH/${RUN}.fastq.gz
echo "[$(date)] downloading $RUN from $URL"
wget -q -O "$FQ" "$URL" || { echo "download failed: $RUN"; exit 1; }

echo "[$(date)] aligning $RUN"
minimap2 -t 4 -ax map-ont "$REF" "$FQ" 2>/dev/null \
  | samtools sort -@ 4 -o "$SCRATCH/$RUN.bam" - || { echo "align failed: $RUN"; exit 1; }
samtools index "$SCRATCH/$RUN.bam"

# mean depth at blaSHV (all positions, incl. zero)
SHV_DEPTH=$(samtools depth -a -r "${CHR}:${SHV_START}-${SHV_END}" "$SCRATCH/$RUN.bam" \
  | awk '{s+=$3; n++} END{ if(n>0) printf "%.4f", s/n; else print 0 }')
# mean chromosomal depth (robust normaliser) via samtools coverage (meandepth col 7)
CHR_MEAN=$(samtools coverage -r "$CHR" "$SCRATCH/$RUN.bam" | awk 'NR==2{print $7}')

CRR=$(awk -v s="$SHV_DEPTH" -v m="$CHR_MEAN" 'BEGIN{ if(m>0) printf "%.4f", s/m; else print "NA" }')
echo -e "${SAMPLE}\t${RUN}\t${SHV_DEPTH}\t${CHR_MEAN}\t${CRR}" > "$OUT"
echo "[$(date)] done $RUN: SHV_depth=$SHV_DEPTH chrom_mean=$CHR_MEAN CRR=$CRR"
