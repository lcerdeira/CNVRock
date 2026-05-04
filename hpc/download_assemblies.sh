#!/bin/bash
#SBATCH --job-name=kpsc_assemblies
#SBATCH --output=/home/lshlt19/CNVRock/logs/assemblies_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/assemblies_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Download AllTheBacteria KpSC assemblies (FASTA) for Kleborate ground truth.
#
# Downloads ATB v0.2 metadata (~1 GB), maps run accessions to assembly FTPs,
# then downloads and decompresses all assemblies in parallel.
#
# Usage:
#   sbatch hpc/download_assemblies.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cnvrock

mkdir -p "$REPO_DIR/data/assemblies" "$REPO_DIR/logs"

echo "Downloading assemblies for KpSC samples …"
python3 "$REPO_DIR/data/setup/download_assemblies.py" \
    --accessions   "$REPO_DIR/assets/kpsc_sra_accessions.txt" \
    --out-dir      "$REPO_DIR/data/assemblies/" \
    --metadata-cache "$REPO_DIR/assets/atb_metadata_v02.tsv.gz" \
    --workers      "$SLURM_CPUS_PER_TASK"

echo ""
echo "Assembly download complete."
echo "Next step: sbatch hpc/get_kleborate_gt.sh"
