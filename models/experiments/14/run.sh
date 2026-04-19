#!/bin/bash
# Full VAE retraining with 05_conv_vae (encoder with residual skip connections, p=0.30).
# Architecture change requires a new checkpoint — cannot reuse exp 13.
# train.py handles training, inference, HMM segmentation, CNV calling, and evaluation.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/14/config.yaml
