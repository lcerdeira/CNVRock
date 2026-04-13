#!/bin/bash
# Re-run HMM segmentation + CNV calling + evaluation on the experiment 06 checkpoint.
# No VAE retraining. HMM self_transition changed (0.75 → 0.80), so segmentation
# must be re-run; CNV calling and evaluation follow automatically.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/06_weighted_sampling/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/07/config.yaml "$CHECKPOINT"
