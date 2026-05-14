#!/bin/bash
# Exp 21: First KpSC experiment — full training from scratch.
# Prerequisites: complete all data preparation steps in README.md before running.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

echo "=== Experiment 21: KpSC chromosomal CNV ==="
echo "Repo root: $REPO_ROOT"

# Guard: verify that the KpSC NPY store exists before spending GPU time
STORE="$REPO_ROOT/data/inputs/KpSC-HS11286-1000bp-core-npy"
if [[ ! -f "$STORE/counts.npy" ]]; then
    echo "ERROR: $STORE/counts.npy not found."
    echo "Complete data preparation steps in models/experiments/21/README.md first."
    exit 1
fi

# Note: gene coordinate validation is performed by _check_coordinates() inside
# 06_gene_cnv_caller.py when CNV calling runs. A clear RuntimeError is raised
# if placeholder coordinates were not updated.

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/21/config.yaml
