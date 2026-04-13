#!/bin/bash
# Re-run CNV calling + evaluation only. No VAE retraining, no HMM re-segmentation.
# cnv_min_confidence is a post-segmentation gate — only the CNV caller and eval
# need to re-run. Reuses exp 07 segments (same HMM config) and exp 06 checkpoint.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/06_weighted_sampling/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/08/config.yaml "$CHECKPOINT"
