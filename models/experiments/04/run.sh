#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

CHECKPOINT="$REPO_ROOT/data/results/01_pf9_conv_vae/checkpoint.pth"

"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/04/config.yaml "$CHECKPOINT"
