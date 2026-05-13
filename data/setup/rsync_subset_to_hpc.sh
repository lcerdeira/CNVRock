#!/bin/bash
# Upload downloaded FASTQ files from the local Aspera dump to the LSHTM HPC.
#
# Designed to run from your Mac after data/setup/aspera_download_subset.sh
# (or in parallel — rsync is restartable and will pick up new files).
#
# Routes via `pryor` (jump host) because loginhpc is not directly reachable.
# Your ~/.ssh/config should already have a `loginhpc` Host entry using
# ProxyJump=pryor.
#
# Usage:
#   bash data/setup/rsync_subset_to_hpc.sh \
#        /Volumes/Orange/cnvrock-dataset/fastq \
#        /home/lshlt19/CNVRock/data/raw/fastq_subset
#
# Optional env vars:
#   BW_KBPS   rsync bandwidth limit in KB/s (default unlimited)
#   DELETE    set to 1 to delete local files after successful upload (default 0)

set -euo pipefail

SRC="${1:?usage: $0 <local_src_dir> <remote_dest_dir>}"
DEST="${2:?usage: $0 <local_src_dir> <remote_dest_dir>}"

BW_KBPS="${BW_KBPS:-0}"
DELETE="${DELETE:-0}"

if [[ ! -d "$SRC" ]]; then
    echo "ERROR: source dir not found: $SRC" >&2
    exit 1
fi

# Ensure remote directory exists
ssh -o ConnectTimeout=10 loginhpc "mkdir -p '$DEST'"

# Build rsync flags
RSYNC_FLAGS=(-avz --partial --inplace --info=progress2 --human-readable
             --include='*.fastq.gz' --include='*/' --exclude='*')
if (( BW_KBPS > 0 )); then
    RSYNC_FLAGS+=(--bwlimit="$BW_KBPS")
fi
if (( DELETE == 1 )); then
    RSYNC_FLAGS+=(--remove-source-files)
fi

# Run via SSH using ProxyJump-aware config
rsync -e "ssh -o ConnectTimeout=10" "${RSYNC_FLAGS[@]}" \
      "$SRC/" "loginhpc:$DEST/"

REMOTE_COUNT=$(ssh loginhpc "find '$DEST' -name '*.fastq.gz' -size +0 | wc -l" | tr -d ' ')
LOCAL_COUNT=$(find "$SRC" -name "*.fastq.gz" -size +0 | wc -l | tr -d ' ')
echo
echo "Local files:  $LOCAL_COUNT"
echo "Remote files: $REMOTE_COUNT"
