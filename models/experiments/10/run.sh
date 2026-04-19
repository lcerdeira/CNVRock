#!/bin/bash
# Re-run HMM, CNV calling, and evaluation on the exp 09 checkpoint.
# No retraining — only the HMM self_transition parameter changed (0.80 → 0.75).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/09_vae_dropout/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/10/config.yaml "$CHECKPOINT"
