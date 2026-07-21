#!/bin/bash
#SBATCH --job-name=saureus_stores
#SBATCH --output=/home/lshlt19/CNVRock/logs/saureus_stores_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/saureus_stores_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal

# Build the two NPY stores for S. aureus (exp 47), once both count passes exist:
#   chromosome store  <- MQ>=20 counts  (saureus_readcounts_mq20)
#   plasmid store     <- MQ=0   counts  (saureus_readcounts_mq0)
#
# Both builders are fully CLI-parameterised (no hardcoded paths), so this just
# points them at the S. aureus reference contigs, panel and manifest.
#
# Run AFTER submit_saureus_mq20.sh has finished.
#   sbatch hpc/build_saureus_stores.sh

set -euo pipefail
REPO_DIR="$HOME/CNVRock"
export PATH="$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"
cd "$REPO_DIR"

REF_DIR="assets/saureus_ref"
MANIFEST="$REF_DIR/saureus_subset.tsv"

# S. aureus contigs: chromosome is NC_007793.1; plasmids are the other four.
CHROM_CONTIG="NC_007793.1"

echo "== chromosome store (MQ>=20) =="
python3 data/setup/readcounts_to_npy_kpsc.py \
    --counts-dir  data/raw/saureus_readcounts_mq20 \
    --manifest    "$MANIFEST" \
    --out-dir     data/inputs/saureus-USA300-mq20-1000bp-npy \
    --keep-contigs "$CHROM_CONTIG" \
    --workers 16

echo "== plasmid store (MQ=0, family-aggregated) =="
python3 data/setup/plasmid_genes_to_npy_kpsc.py \
    --counts-dir  data/raw/saureus_readcounts_mq0 \
    --manifest    "$MANIFEST" \
    --gene-coords "$REF_DIR/gene_coords_plasmid.tsv" \
    --families    "$REF_DIR/gene_families.tsv" \
    --out-dir     data/inputs/saureus-USA300-mq0-plasmid-1000bp-npy \
    --workers 16

echo "== done =="
ls -la data/inputs/saureus-USA300-*/
