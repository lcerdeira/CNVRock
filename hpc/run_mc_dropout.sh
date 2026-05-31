#!/bin/bash
#SBATCH --job-name=mc_dropout_exp33
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=slurm-mc-dropout-%j.out

set -euo pipefail

source /home/lshlt19/miniconda3/etc/profile.d/conda.sh
conda activate cnvrock

cd /home/lshlt19/CNVRock

echo "[$(date)] MC Dropout job starting on $(hostname)"
echo "  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'none')"
echo "  CUDA: $(python3 -c 'import torch; print(torch.cuda.is_available())')"

python3 analysis/mc_dropout_uncertainty.py \
    --exp       33 \
    --n-mc      50 \
    --n-samples 1000 \
    --hpc-root  /home/lshlt19/CNVRock

echo "[$(date)] Done."
