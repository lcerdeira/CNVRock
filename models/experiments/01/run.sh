#!/bin/bash
# Run this experiment locally on a Mac (uses MPS if available, otherwise CPU).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/01/config.yaml
