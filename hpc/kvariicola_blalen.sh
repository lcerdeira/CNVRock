#!/bin/bash
#SBATCH --job-name=kvar_blalen
#SBATCH --output=/home/lshlt19/CNVRock/logs/kvar_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/kvar_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#
# Multi-reference chromosomal calling demonstration (limitation 5a).
#
# blaSHV is called only for K. pneumoniae/quasipneumoniae because K. variicola
# carries the LEN-family homolog (blaLEN) at the syntenic locus, whose reads
# cross-map onto the HS11286 blaSHV coordinates and produce spurious calls. A
# species-appropriate reference resolves this. Per K. variicola isolate we
# recover the reads from the retained HS11286 BAM, RE-ALIGN them to a
# K. variicola reference (NC_011283.1) with minimap2 (short-read preset), and
# compute the blaLEN copy-ratio from depth.
#
#   sbatch --array=1-$(wc -l < assets/kvariicola_10k_ids.txt)%25 hpc/kvariicola_blalen.sh

set -uo pipefail
REPO=/home/lshlt19/CNVRock
export PATH="/home/lshlt19/miniforge3/envs/lralign/bin:$PATH"   # minimap2 + samtools>=1.17
REF=$REPO/assets/kvariicola_ref/kvar.fna
IDS=$REPO/assets/kvariicola_10k_ids.txt
OUTDIR=$REPO/data/results/kvariicola_multiref
mkdir -p "$OUTDIR" "$REPO/logs"
SCRATCH=/home/lshlt19/scratch/kvar/${SLURM_ARRAY_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

CHR=NC_011283.1
# two blaLEN-family loci (blaSHV CDS maps to both; mapq 0 = 2 copies)
LEN1=2835133-2835994
LEN2=2852061-2852922

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS")
BAM=$REPO/data/raw/bam_subset/${ID}.bam
OUT=$OUTDIR/${ID}.tsv
[ -s "$OUT" ] && { echo "done $ID"; exit 0; }
[ -f "$BAM" ] || { echo "no BAM for $ID"; exit 0; }

# recover reads from the HS11286 BAM and re-align to the K. variicola reference
samtools fastq -@ 4 "$BAM" 2>/dev/null > "$SCRATCH/reads.fq" || { echo "fastq fail $ID"; exit 1; }
minimap2 -t 4 -ax sr "$REF" "$SCRATCH/reads.fq" 2>/dev/null \
  | samtools sort -@ 4 -T "$SCRATCH/st" -o "$SCRATCH/$ID.bam" - || { echo "align fail $ID"; exit 1; }
samtools index "$SCRATCH/$ID.bam"

d1=$(samtools depth -a -r "${CHR}:${LEN1}" "$SCRATCH/$ID.bam" | awk '{s+=$3;n++}END{if(n)printf "%.3f",s/n;else print 0}')
d2=$(samtools depth -a -r "${CHR}:${LEN2}" "$SCRATCH/$ID.bam" | awk '{s+=$3;n++}END{if(n)printf "%.3f",s/n;else print 0}')
CHRM=$(samtools coverage -r "$CHR" "$SCRATCH/$ID.bam" | awk 'NR==2{print $7}')
CRR=$(awk -v a="$d1" -v b="$d2" -v m="$CHRM" 'BEGIN{g=(a+b)/2; if(m>0)printf "%.4f",g/m; else print "NA"}')
echo -e "${ID}\t${d1}\t${d2}\t${CHRM}\t${CRR}" > "$OUT"
echo "done $ID: blaLEN_CRR=$CRR (loci depths $d1,$d2 | chrom $CHRM)"
