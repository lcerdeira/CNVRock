#!/bin/bash
# Submit Phase 2 expansion jobs to SLURM.
#
# Three parallel streams:
#   A) download_expansion_sra.sh    — download FASTQs + align (88k BAMs)
#   B) download_assemblies_s3.sh    — download assemblies from S3 (Kleborate)
#   C) collect_expansion_readcounts.sh — GATK readcounts (depends on A)
#
# MaxArraySize=5000, so A and C are split into chunks of 5000 with BATCH_OFFSET.
# B reads from assembly_urls.tsv (1-indexed with header skip), also chunked.
#
# Usage:
#   cd ~/CNVRock
#   bash hpc/submit_phase2_expansion.sh

set -uo pipefail  # no -e: let sbatch_retry handle failures without aborting the script
cd "$(dirname "$0")/.."

# Retry sbatch up to 5 times with exponential back-off (SLURM transient errors)
sbatch_retry() {
    local attempt=0
    while [[ $attempt -lt 5 ]]; do
        out=$(sbatch "$@" 2>&1)
        if echo "$out" | grep -qE '^[0-9]+$'; then
            echo "$out"
            return 0
        fi
        attempt=$(( attempt + 1 ))
        echo "  sbatch attempt $attempt failed: $out — retrying in $(( 60 * attempt )) s …" >&2
        sleep $(( 60 * attempt ))
    done
    echo "ERROR: sbatch failed after 5 attempts" >&2
    return 1
}

REPO_DIR="$(pwd)"
N_SAMPLES=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)
N_ASM=$(( $(wc -l < assets/kpsc_expansion_assembly_urls.tsv) - 1 ))  # skip header
CHUNK=4999  # MaxArraySize=5000 means max task index is 4999

# Optional: resume from a specific offset (e.g. if first N chunks already submitted)
#   bash hpc/submit_phase2_expansion.sh 4999   ← skip first chunk (offset 0)
START_OFFSET=${1:-0}

echo "Phase 2 expansion submission"
echo "  Samples (SRA):      $N_SAMPLES"
echo "  Samples (assembly): $N_ASM"
echo "  Chunk size:         $CHUNK"
[[ $START_OFFSET -gt 0 ]] && echo "  Resuming from offset: $START_OFFSET"
echo ""

# ── Stream A: download FASTQs + align ────────────────────────────────────
echo "=== Stream A: BAM download + align ==="
DL_JOBS=()
OFFSET=$START_OFFSET
while [[ $OFFSET -lt $N_SAMPLES ]]; do
    REMAINING=$(( N_SAMPLES - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    JOB_ID=$(sbatch_retry --parsable \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%50 \
        hpc/download_expansion_sra.sh) || { echo "  SKIP chunk offset=$OFFSET (sbatch failed)"; OFFSET=$(( OFFSET + CHUNK )); continue; }
    DL_JOBS+=("$JOB_ID")
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE"
    OFFSET=$(( OFFSET + CHUNK ))
    sleep 60  # give SLURM time to register the array before next submission
done

# ── Stream B: assembly download (independent) ────────────────────────────
echo ""
echo "=== Stream B: assembly download ==="
ASM_JOBS=()
OFFSET=$START_OFFSET
while [[ $OFFSET -lt $N_ASM ]]; do
    REMAINING=$(( N_ASM - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    # Assembly URL file has a header; row i (1-indexed, no header) = line i+1 in file.
    # The download_assemblies_s3.sh already adds +1 to skip header, so pass BATCH_OFFSET as-is.
    JOB_ID=$(sbatch_retry --parsable \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%100 \
        hpc/download_assemblies_s3.sh) || { echo "  SKIP chunk offset=$OFFSET (sbatch failed)"; OFFSET=$(( OFFSET + CHUNK )); continue; }
    ASM_JOBS+=("$JOB_ID")
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE"
    OFFSET=$(( OFFSET + CHUNK ))
    sleep 60
done

# ── Stream C: read counts (depends on ALL of Stream A) ───────────────────
echo ""
echo "=== Stream C: read count extraction (after Stream A) ==="
# Build colon-separated dependency string for all DL jobs
DEP=$(IFS=:; echo "${DL_JOBS[*]}")
OFFSET=$START_OFFSET
while [[ $OFFSET -lt $N_SAMPLES ]]; do
    REMAINING=$(( N_SAMPLES - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    JOB_ID=$(sbatch_retry --parsable \
        --dependency=afterok:${DEP} \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%100 \
        hpc/collect_expansion_readcounts.sh) || { echo "  SKIP chunk offset=$OFFSET (sbatch failed)"; OFFSET=$(( OFFSET + CHUNK )); continue; }
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE  (after ${DEP})"
    OFFSET=$(( OFFSET + CHUNK ))
    sleep 60
    sleep 2
done

echo ""
echo "All jobs submitted. Check with: squeue -u lshlt19"
