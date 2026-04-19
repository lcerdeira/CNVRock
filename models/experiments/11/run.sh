#!/bin/bash
# Full VAE retraining with 03_conv_vae (encoder dropout p=0.10).
# HMM segmentation, CNV calling, and evaluation follow automatically via train.py.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/11/config.yaml
