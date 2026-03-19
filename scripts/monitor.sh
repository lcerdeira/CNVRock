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
    -G team342 \
    -R "select[mem>=1000] rusage[mem=1000] span[hosts=1]" \
    -q basement "python3 -u ../scripts/monitor.py $JOBID"