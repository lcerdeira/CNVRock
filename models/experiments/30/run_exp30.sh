#!/bin/bash
# Exp 30: hold-out validation — generate split then train on 80% subset.
#
# Run from repo root on HPC login node:
#   bash models/experiments/30/run_exp30.sh
#
# The script:
#   1. Generates the 80/20 stratified hold-out split (idempotent — skips if exists)
#   2. Submits the GPU training job for exp 30

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
CONFIG="$REPO_DIR/models/experiments/30/config.yaml"
HOLDOUT_FILE="$REPO_DIR/assets/holdout_sample_ids.txt"

echo "=== Exp 30: KpSC hold-out validation ==="

# Step 1: generate split (only if not already done)
if [[ -f "$HOLDOUT_FILE" ]]; then
    N=$(wc -l < "$HOLDOUT_FILE")
    echo "Hold-out split already exists ($N samples) — skipping generation."
else
    echo "Generating stratified 80/20 hold-out split..."
    cd "$REPO_DIR"
    python data/setup/create_holdout_split.py \
        --gt   assets/kpsc_amrfinder_ground_truth.tsv \
        --store data/inputs/KpSC-HS11286-1000bp-core-npy \
        --seed 42 --holdout-frac 0.20
fi

# Step 2: submit GPU training job
echo ""
echo "Submitting training job..."
JOB_ID=$(sbatch --parsable hpc/train_gpu.sh "$CONFIG")
echo "Submitted job $JOB_ID — monitor with: tail -f logs/train_${JOB_ID}.out"
echo ""
echo "Results will appear in: data/results/30_kpsc_holdout_v1/evaluation.txt"
