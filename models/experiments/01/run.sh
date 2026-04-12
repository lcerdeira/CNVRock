#!/bin/bash
# Run this experiment locally on a Mac (uses MPS if available, otherwise CPU).
set -euo pipefail

cd "$(dirname "$0")"

.venv/bin/python ../../train.py config.yaml
