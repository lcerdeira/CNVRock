#!/usr/bin/env bash
set -euo pipefail

bsub \
    -J "rc2npy" \
    -o "log.o" \
    -e "log.e" \
    -q basement \
    -G team342 \
    -R "select[mem>48000] rusage[mem=48000]" \
    -M 48000 \
    -n 8 \
    "PYTHONUNBUFFERED=1 python3 -u readcounts_to_npy.py"