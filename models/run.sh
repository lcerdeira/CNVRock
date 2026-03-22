#!/bin/bash
set -euo pipefail

JOBID=$(bsub \
    -J gw-cnv \
    -o log.o \
    -e log.e \
    -n 1 \
    -G team342 \
    -M 4000 \
    -R 'select[mem>=4000] rusage[mem=4000] span[hosts=1]' \
    -gpu "num=1:mode=exclusive_process" \
    -q gpu-basement \
    "python main.py" 2>&1 | grep -oP '(?<=Job <)\d+(?=>)')

echo "Submitted job $JOBID"
bash ../scripts/monitor.sh $JOBID