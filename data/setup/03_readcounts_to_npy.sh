#!/usr/bin/env bash
set -euo pipefail

bsub \
    -J "rc2npy" \
    -o "log.o" \
    -e "log.e" \
    -q basement \
    -G team342 \
    -R "select[mem>16000] rusage[mem=16000]" \
    -M 16000 \
    -n 4 \
    "PYTHONUNBUFFERED=1 .venv/bin/python -u readcounts_to_npy.py"