#!/bin/bash
# Submit the MQ>=20 chromosomal recount for the S. aureus cohort.
#
# The MQ=0 pass (aspera_subset_pipeline.sh) gives the plasmid store, where
# multi-mapping reads must be kept. The chromosomal genes (mecA, norA, pbp4,
# gdpP) want the cleaner MQ>=20 signal, recounted from the BAMs that pass
# already produced — no re-download, no re-align.
#
# Run ONLY after the alignment array (job 6401672) has finished, so every
# BAM exists. recount_mq20_from_bams.sh is idempotent and env-var driven, so
# nothing in it needs editing — this wrapper just points it at the S. aureus
# reference, intervals, BAMs and a fresh output dir, and builds the BAM list.
#
# Usage (HPC login node):
#   bash hpc/submit_ecoli_mq20.sh

set -euo pipefail
REPO_DIR="$HOME/CNVRock"
BAM_DIR="$REPO_DIR/data/raw/ecoli_bam"
LIST="$REPO_DIR/assets/_bams_list_ecoli_mq20.txt"

ls "$BAM_DIR"/*.bam > "$LIST"
N=$(wc -l < "$LIST")
echo "E. coli BAMs to recount at MQ>=20: $N"

sbatch --array=1-"${N}"%50 \
  --export=ALL,\
REFERENCE="$REPO_DIR/assets/ecoli_ref/EC958_extended.fasta",\
INTERVALS="$REPO_DIR/assets/ecoli_ref/EC958_extended_1kb.interval_list",\
BAM_DIR="$BAM_DIR",\
COUNTS_DIR="$REPO_DIR/data/raw/ecoli_readcounts_mq20",\
BAMS_LIST="$LIST",\
MIN_MQ=20 \
  hpc/recount_mq20_from_bams.sh
