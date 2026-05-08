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

set -euo pipefail
cd "$(dirname "$0")/.."

REPO_DIR="$(pwd)"
N_SAMPLES=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)
N_ASM=$(( $(wc -l < assets/kpsc_expansion_assembly_urls.tsv) - 1 ))  # skip header
CHUNK=5000

echo "Phase 2 expansion submission"
echo "  Samples (SRA):     $N_SAMPLES"
echo "  Samples (assembly): $N_ASM"
echo "  Chunk size:         $CHUNK"
echo ""

# ── Stream A: download FASTQs + align ────────────────────────────────────
echo "=== Stream A: BAM download + align ==="
DL_JOBS=()
OFFSET=0
while [[ $OFFSET -lt $N_SAMPLES ]]; do
    REMAINING=$(( N_SAMPLES - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    JOB_ID=$(sbatch --parsable \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%50 \
        hpc/download_expansion_sra.sh)
    DL_JOBS+=("$JOB_ID")
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE"
    OFFSET=$(( OFFSET + CHUNK ))
done

# ── Stream B: assembly download (independent) ────────────────────────────
echo ""
echo "=== Stream B: assembly download ==="
ASM_JOBS=()
OFFSET=0
while [[ $OFFSET -lt $N_ASM ]]; do
    REMAINING=$(( N_ASM - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    # Assembly URL file has a header; row i (1-indexed, no header) = line i+1 in file.
    # The download_assemblies_s3.sh already adds +1 to skip header, so pass BATCH_OFFSET as-is.
    JOB_ID=$(sbatch --parsable \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%100 \
        hpc/download_assemblies_s3.sh)
    ASM_JOBS+=("$JOB_ID")
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE"
    OFFSET=$(( OFFSET + CHUNK ))
done

# ── Stream C: read counts (depends on ALL of Stream A) ───────────────────
echo ""
echo "=== Stream C: read count extraction (after Stream A) ==="
# Build colon-separated dependency string for all DL jobs
DEP=$(IFS=:; echo "${DL_JOBS[*]}")
OFFSET=0
while [[ $OFFSET -lt $N_SAMPLES ]]; do
    REMAINING=$(( N_SAMPLES - OFFSET ))
    SIZE=$(( REMAINING > CHUNK ? CHUNK : REMAINING ))
    JOB_ID=$(sbatch --parsable \
        --dependency=afterok:${DEP} \
        --export=BATCH_OFFSET=${OFFSET} \
        --array=1-${SIZE}%100 \
        hpc/collect_expansion_readcounts.sh)
    echo "  Submitted: job $JOB_ID  offset=$OFFSET  tasks=$SIZE  (after ${DEP})"
    OFFSET=$(( OFFSET + CHUNK ))
done

echo ""
echo "All jobs submitted. Check with: squeue -u lshlt19"
