#!/bin/bash
# Re-run CNV calling + evaluation on the experiment 01 checkpoint + exp 04 HMM segments.
# No retraining, no re-segmentation — only CNV caller version and thresholds changed.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/01_pf9_conv_vae/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/05/config.yaml "$CHECKPOINT"
