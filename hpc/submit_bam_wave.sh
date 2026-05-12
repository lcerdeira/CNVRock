#!/bin/bash
# submit_bam_wave.sh  — Wave-based BAM download submission for expansion cohort.
#
# Prerequisites:
#   1. Build URL manifest: python3 data/setup/build_ena_url_manifest.py
#   2. This script is started on the HPC login node (NOT a compute node).
#
# Usage:
#   cd ~/CNVRock
#   nohup bash hpc/submit_bam_wave.sh > ~/bam_wave.log 2>&1 &
#
# Keeps ≤5000 tasks in the queue at a time (well under MaxJobCount=10000).
# Skips accessions already done (BAM + .bai exist).

set -euo pipefail
cd ~/CNVRock

SCRIPT="hpc/download_expansion_sra.sh"
MANIFEST="assets/ena_url_manifest.tsv"

CHUNK=100       # tasks per SLURM array job
WAVE=50         # chunks per wave  (50×100 = 5000 tasks)
SLEEP_BETWEEN=2 # seconds between sbatch calls
CONCUR=10       # max concurrent tasks per array (limits EBI connections)

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: manifest not found at $MANIFEST"
    echo "Run first: python3 data/setup/build_ena_url_manifest.py"
    exit 1
fi

# Number of data rows (exclude header)
N=$(( $(wc -l < "$MANIFEST") - 1 ))
echo "Manifest rows: $N"
echo "Started: $(date)"

wait_for_jobs() {
    local wave_jobs=("$@")
    echo "  Waiting for ${#wave_jobs[@]} jobs…"
    while true; do
        local still_running=0
        for jid in "${wave_jobs[@]}"; do
            if squeue -j "$jid" -h 2>/dev/null | grep -q .; then
                still_running=1
                break
            fi
        done
        [[ $still_running -eq 0 ]] && break
        sleep 120
    done
}

OFFSET=0
WAVE_NUM=0
while [[ $OFFSET -lt $N ]]; do
    WAVE_NUM=$(( WAVE_NUM + 1 ))
    echo "--- Wave $WAVE_NUM (offset=$OFFSET): $(date) ---"
    WAVE_JOBS=()
    for (( chunk=0; chunk<WAVE; chunk++ )); do
        START=$(( OFFSET + chunk * CHUNK ))
        [[ $START -ge $N ]] && break
        REMAINING=$(( N - START ))
        SIZE=$(( REMAINING < CHUNK ? REMAINING : CHUNK ))
        JID=$(BATCH_OFFSET=$START sbatch --parsable \
            --array=1-${SIZE}%${CONCUR} \
            "$SCRIPT")
        echo "  dl: job $JID  offset=$START  size=$SIZE"
        WAVE_JOBS+=("$JID")
        sleep $SLEEP_BETWEEN
    done
    OFFSET=$(( OFFSET + WAVE * CHUNK ))
    wait_for_jobs "${WAVE_JOBS[@]}"
    echo "  Wave $WAVE_NUM done: $(date)"

    # Progress report
    N_BAMS=$(ls data/raw/bam/*.bam 2>/dev/null | wc -l || echo 0)
    echo "  BAMs so far: $N_BAMS"
done

echo ""
echo "All BAM download waves complete: $(date)"
N_BAMS=$(ls data/raw/bam/*.bam 2>/dev/null | wc -l || echo 0)
echo "Total BAMs: $N_BAMS / $N"
