#!/bin/bash
# Submit Phase 2 expansion jobs to SLURM.
#
# Three parallel streams:
#   A) download_expansion_sra.sh    — download FASTQs + align (88k BAMs)
#   B) download_assemblies_s3.sh    — download assemblies from S3 (Kleborate)
#   C) collect_expansion_readcounts.sh — GATK readcounts (depends on all of A)
#
# Chunks within Stream A are chained (each depends on the previous), and
# Stream B chunks are chained separately.  This keeps at most 2 large arrays
# in the queue at once, staying well under MaxJobCount=10000.
# A and B run concurrently because they are independent.
# C starts only after the last chunk of A completes.
#
# Usage:
#   cd ~/CNVRock
#   bash hpc/submit_phase2_expansion.sh [A_PREV_JOB_ID] [B_PREV_JOB_ID] [A_START_OFFSET] [B_START_OFFSET]
#
# Resume example (skip already-running chunks 1 of A and B):
#   bash hpc/submit_phase2_expansion.sh 5372994 ""   # A chunk1 running; submit A chunk2+ and all B
#
# Resume with non-standard offsets (e.g. two old 4999-task jobs already cover 9998 samples):
#   bash hpc/submit_phase2_expansion.sh 5376786 "" 9998  # A: continue from offset 9998, chain after job 5376786

set -uo pipefail
cd "$(dirname "$0")/.."

# Retry sbatch up to 8 times with growing back-off (SLURM transient errors)
sbatch_retry() {
    local attempt=0
    while [[ $attempt -lt 8 ]]; do
        out=$(sbatch "$@" 2>&1)
        if echo "$out" | grep -qE '^[0-9]+$'; then
            echo "$out"
            return 0
        fi
        attempt=$(( attempt + 1 ))
        local delay=$(( 30 * attempt ))
        echo "  sbatch attempt $attempt failed — retrying in ${delay}s …" >&2
        sleep $delay
    done
    echo "ERROR: sbatch failed after 8 attempts" >&2
    return 1
}

REPO_DIR="$(pwd)"
N_SAMPLES=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)
N_ASM=$(( $(wc -l < assets/kpsc_expansion_assembly_urls.tsv) - 1 ))
CHUNK=500   # 500-task arrays submit reliably; 4999-task arrays cause SLURM controller overload

# Optional: pass job IDs of already-running first chunks to chain from them
A_PREV_JOB="${1:-}"      # e.g. 5376786 (last A job already submitted)
B_PREV_JOB="${2:-}"      # e.g. 5380000 (last B job already submitted)
A_START_OFFSET="${3:-}"  # explicit offset to resume A from (default: CHUNK)
B_START_OFFSET="${4:-}"  # explicit offset to resume B from (default: CHUNK)

echo "Phase 2 expansion submission (chained chunks)"
echo "  Samples (SRA):       $N_SAMPLES"
echo "  Samples (assembly):  $N_ASM"
echo "  Chunk size:          $CHUNK"
echo "  Chunks (A):          $(( (N_SAMPLES + CHUNK - 1) / CHUNK ))"
echo "  Chains: A→A→…→A runs concurrently with B→B→…→B; C starts after last A"
[[ -n "$A_PREV_JOB" ]] && echo "  A: continuing after job $A_PREV_JOB (offset=${A_START_OFFSET:-$CHUNK})"
[[ -n "$B_PREV_JOB" ]] && echo "  B: continuing after job $B_PREV_JOB (offset=${B_START_OFFSET:-$CHUNK})"
echo ""

# ── Stream A: download FASTQs + align (chained) ──────────────────────────
echo "=== Stream A: BAM download + align ==="
DL_JOBS=()
OFFSET=0
[[ -n "$A_PREV_JOB" ]] && { DL_JOBS+=("$A_PREV_JOB"); OFFSET=${A_START_OFFSET:-$CHUNK}; }

while [[ $OFFSET -lt $N_SAMPLES ]]; do
    REMAINING=$(( N_SAMPLES - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))

    # Chain: this chunk depends on the previous A chunk
    if [[ ${#DL_JOBS[@]} -gt 0 ]]; then
        PREV="${DL_JOBS[-1]}"
        DEP_FLAG="--dependency=afterok:${PREV}"
    else
        DEP_FLAG=""
    fi

    JOB_ID=$(sbatch_retry --parsable \
        $DEP_FLAG \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%${SIZE} \
        hpc/download_expansion_sra.sh) || {
            echo "  FAILED chunk A offset=$OFFSET after 8 retries — stopping A chain" >&2
            break
        }
    DL_JOBS+=("$JOB_ID")
    echo "  Submitted A: job $JOB_ID  offset=$OFFSET  tasks=$SIZE${DEP_FLAG:+  (after $PREV)}"
    OFFSET=$(( OFFSET + CHUNK ))
    sleep 5  # brief pause; no need to wait long since chunks are chained
done

# ── Stream B: assembly download (chained, independent of A) ──────────────
echo ""
echo "=== Stream B: assembly download ==="
ASM_JOBS=()
OFFSET=0
[[ -n "$B_PREV_JOB" ]] && { ASM_JOBS+=("$B_PREV_JOB"); OFFSET=${B_START_OFFSET:-$CHUNK}; }

while [[ $OFFSET -lt $N_ASM ]]; do
    REMAINING=$(( N_ASM - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))

    if [[ ${#ASM_JOBS[@]} -gt 0 ]]; then
        PREV="${ASM_JOBS[-1]}"
        DEP_FLAG="--dependency=afterok:${PREV}"
    else
        DEP_FLAG=""
    fi

    JOB_ID=$(sbatch_retry --parsable \
        $DEP_FLAG \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%${SIZE} \
        hpc/download_assemblies_s3.sh) || {
            echo "  FAILED chunk B offset=$OFFSET after 8 retries — stopping B chain" >&2
            break
        }
    ASM_JOBS+=("$JOB_ID")
    echo "  Submitted B: job $JOB_ID  offset=$OFFSET  tasks=$SIZE${DEP_FLAG:+  (after $PREV)}"
    OFFSET=$(( OFFSET + CHUNK ))
    sleep 5
done

# ── Stream C: read counts (depends on LAST chunk of A) ───────────────────
echo ""
echo "=== Stream C: read count extraction ==="
if [[ ${#DL_JOBS[@]} -eq 0 ]]; then
    echo "  No A jobs submitted — skipping C" >&2
else
    # C chunks are chained off the last A job, then chained among themselves
    LAST_A="${DL_JOBS[-1]}"
    C_PREV=""
    OFFSET=0
    while [[ $OFFSET -lt $N_SAMPLES ]]; do
        REMAINING=$(( N_SAMPLES - OFFSET ))
        SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))

        if [[ -n "$C_PREV" ]]; then
            DEP_FLAG="--dependency=afterok:${C_PREV}"
        else
            DEP_FLAG="--dependency=afterok:${LAST_A}"
        fi

        JOB_ID=$(sbatch_retry --parsable \
            $DEP_FLAG \
            --export=BATCH_OFFSET=${OFFSET} \
            --array=1-${SIZE}%100 \
            hpc/collect_expansion_readcounts.sh) || {
                echo "  FAILED chunk C offset=$OFFSET — stopping C chain" >&2
                break
            }
        C_PREV="$JOB_ID"
        echo "  Submitted C: job $JOB_ID  offset=$OFFSET  tasks=$SIZE  (after ${DEP_FLAG#*:})"
        OFFSET=$(( OFFSET + CHUNK ))
        sleep 5
    done
fi

echo ""
echo "=== Summary ==="
echo "Stream A jobs (${#DL_JOBS[@]}): ${DL_JOBS[*]:-none}"
echo "Stream B jobs (${#ASM_JOBS[@]}): ${ASM_JOBS[*]:-none}"
echo "Stream C: chained after last A"
echo ""
echo "Check with:  squeue -u lshlt19"
