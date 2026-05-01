#!/bin/bash
#SBATCH --job-name=kpsc_amrfinder_gt
#SBATCH --output=/home/lshlt19/CNVRock/logs/amrfinder_gt_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/amrfinder_gt_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

set -euo pipefail

ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

cd /home/lshlt19/CNVRock

python3 data/setup/get_amrfinder_gt.py
