#!/bin/bash
# Full VAE retrain + HMM segmentation + CNV calling + evaluation.
# max_beta raised to 4.0 — retraining required.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/06/config.yaml
