#!/bin/bash
set -euo pipefail

JOBID=$1

source ../.env

bsub \
    -J monitor \
    -o monitor-logs.o \
    -e monitor-logs.e \
    -n 1 \
    -M 1000 \
    -R "select[mem>=1000] rusage[mem=1000] span[hosts=1]" \
    -q basement "python -u monitor.py $JOBID"
