#!/bin/bash
# Fine-tune exp 14 checkpoint with sin loss regularisation.
# Runs the full pipeline: fine-tune → inference → HMM → CNV calling → evaluation.
# Note: when reusing only CNV/eval steps (no retraining), use symlinks (ln -sf)
# instead of cp to reference parent experiment outputs — avoids duplicating large files.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/18/config.yaml
