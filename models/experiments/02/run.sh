#!/bin/bash
# Re-run HMM + CNV calling + evaluation on the experiment 01 checkpoint.
# No retraining — only HMM stickiness and CNV thresholds changed.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/01_pf9_conv_vae/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/02/config.yaml "$CHECKPOINT"
