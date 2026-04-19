#!/bin/bash
# Full VAE retraining with 04_conv_vae (encoder dropout p=0.30).
# Architecture change requires a new checkpoint — cannot reuse exp 09.
# train.py handles training, inference, HMM segmentation, CNV calling, and evaluation.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/12/config.yaml
