#!/bin/bash
# Run this experiment locally on a Mac (uses MPS if available, otherwise CPU).
set -euo pipefail

cd "$(dirname "$0")"

python ../../train.py config.yaml
