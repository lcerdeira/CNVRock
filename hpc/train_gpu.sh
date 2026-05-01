#!/bin/bash
#SBATCH --job-name=cnvrock_train
#SBATCH --output=/home/lshlt19/CNVRock/logs/train_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/train_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# Usage:
#   cd ~/CNVRock
#   sbatch hpc/train_gpu.sh models/experiments/21/config.yaml

set -euo pipefail

CONFIG="${1:-models/experiments/21/config.yaml}"
REPO_DIR="/home/lshlt19/CNVRock"
CONDA_DIR="/home/lshlt19/miniconda3"
ENV_BIN="$CONDA_DIR/envs/cnvrock/bin"

export PATH="$ENV_BIN:$CONDA_DIR/bin:$PATH"
module load cuda/12.4.0

echo "=== CNVRock Training ==="
echo "Config:  $CONFIG"
echo "GPU:     $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Python:  $(python --version)"
echo "Torch:   $(python -c 'import torch; print(torch.__version__, "| CUDA:", torch.cuda.is_available())')"
echo "Started: $(date)"
echo ""

cd "$REPO_DIR/models"
python train.py "$CONFIG"

echo ""
echo "Finished: $(date)"
