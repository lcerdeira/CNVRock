#!/bin/bash
#SBATCH --job-name=vae_ablation
#SBATCH --output=/home/lshlt19/CNVRock/logs/vae_ablation_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/vae_ablation_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# VAE-baseline ablation: simple-median vs housekeeping vs VAE expected-depth.
set -eo pipefail
REPO=/home/lshlt19/CNVRock
cd "$REPO"
"$HOME/miniconda3/envs/cnvrock/bin/python" analysis/vae_ablation.py \
    --store-dir   "$REPO/data/inputs/KpSC-expansion-10k-mq20-1000bp-npy" \
    --results-dir "$REPO/data/results/33_kpsc_expansion_10k" \
    --meta        "$REPO/assets/kpsc_expansion_metadata_runlevel.tsv" \
    --cabbage     "$REPO/assets/cabbage_kpsc_phenotypes.tsv" \
    --kleborate   "$REPO/assets/kpsc_expansion_kleborate_gt_runlevel.tsv" \
    --out-dir     "$REPO/data/results/vae_ablation"
