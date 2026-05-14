#!/bin/bash
# HMM self_transition reduction (0.80 → 0.65) using the exp 12 checkpoint.
# No retraining required — wrap_up.py runs inference, segmentation, CNV calling,
# and evaluation on the existing checkpoint.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up \
    experiments/13/config.yaml \
    "$REPO_ROOT/data/results/12_higher_dropout/checkpoint.pth"
