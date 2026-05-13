#!/bin/bash
# Aspera-based FASTQ download for the KpSC expansion subset.
#
# Runs on your local Mac (not HPC) because EBI rate-limits the HPC IP.
# Downloads R1 + R2 FASTQ.gz files to $DEST_DIR for every accession in the
# subset manifest.
#
# Prerequisites:
#   1. IBM Aspera Connect installed:
#        https://www.ibm.com/aspera/connect/
#      The script auto-detects `ascp` at the standard Mac install paths.
#   2. Subset manifest TSV (cols: accession, layout, r1_url, r2_url) from
#      data/setup/select_expansion_subset.py.
#
# Usage:
#   bash data/setup/aspera_download_subset.sh \
#        assets/kpsc_expansion_subset_5k.tsv \
#        /Volumes/Orange/cnvrock-dataset/fastq
#
# Optional env vars:
#   PARALLEL   max concurrent ascp processes (default 8)
#   BW_MBPS    per-connection bandwidth cap, Mbps (default 200)
#   ASCP_PATH  override ascp binary location
#   ASCP_KEY   override Aspera private-key file

set -euo pipefail

MANIFEST="${1:?usage: $0 <subset.tsv> <dest_dir>}"
DEST="${2:?usage: $0 <subset.tsv> <dest_dir>}"

PARALLEL="${PARALLEL:-8}"
BW_MBPS="${BW_MBPS:-200}"

# ── Find ascp binary ──────────────────────────────────────────────────────────
if [[ -n "${ASCP_PATH:-}" ]]; then
    ASCP="$ASCP_PATH"
elif [[ -x "/Applications/Aspera Connect.app/Contents/Resources/ascp" ]]; then
    ASCP="/Applications/Aspera Connect.app/Contents/Resources/ascp"
elif [[ -x "$HOME/Applications/Aspera Connect.app/Contents/Resources/ascp" ]]; then
    ASCP="$HOME/Applications/Aspera Connect.app/Contents/Resources/ascp"
elif command -v ascp > /dev/null; then
    ASCP="$(command -v ascp)"
else
    echo "ERROR: ascp not found. Install IBM Aspera Connect or set ASCP_PATH." >&2
    exit 1
fi

# ── Find Aspera SSH key ───────────────────────────────────────────────────────
if [[ -n "${ASCP_KEY:-}" ]]; then
    KEY="$ASCP_KEY"
else
    for candidate in \
        "/Applications/Aspera Connect.app/Contents/Resources/asperaweb_id_dsa.openssh" \
        "$HOME/Applications/Aspera Connect.app/Contents/Resources/asperaweb_id_dsa.openssh"; do
        if [[ -f "$candidate" ]]; then
            KEY="$candidate"
            break
        fi
    done
    if [[ -z "${KEY:-}" ]]; then
        echo "ERROR: Aspera SSH key not found. Set ASCP_KEY." >&2
        exit 1
    fi
fi

echo "ascp:    $ASCP"
echo "key:     $KEY"
echo "dest:    $DEST"
echo "parallel: $PARALLEL  (per-conn $BW_MBPS Mbps)"

mkdir -p "$DEST"

# ── Build job list: one line per file to fetch ────────────────────────────────
JOBS=$(mktemp)
trap 'rm -f "$JOBS"' EXIT

awk -F'\t' 'NR>1 {
    acc=$1; layout=$2; r1=$3; r2=$4
    # Strip "https://ftp.sra.ebi.ac.uk/" → leaves "vol1/fastq/..." which is the Aspera path
    sub(/^https:\/\/ftp\.sra\.ebi\.ac\.uk\//, "", r1)
    sub(/^https:\/\/ftp\.sra\.ebi\.ac\.uk\//, "", r2)
    print acc "\t" r1
    if (layout == "PAIRED" && r2 != "") print acc "\t" r2
}' "$MANIFEST" > "$JOBS"

TOTAL=$(wc -l < "$JOBS")
echo "Files queued: $TOTAL"

# ── Worker function: download one file via ascp ───────────────────────────────
download_one() {
    local acc="$1" aspath="$2"
    local fname="${aspath##*/}"
    local out="$DEST/$fname"
    if [[ -s "$out" ]]; then
        return 0   # already done
    fi
    "$ASCP" -P 33001 -O 33001 -QT -l "${BW_MBPS}m" \
            -i "$KEY" \
            "era-fasp@fasp.sra.ebi.ac.uk:${aspath}" "$out" \
            > /dev/null 2>&1
    if [[ -s "$out" ]]; then
        echo "ok   $acc  $fname"
    else
        echo "FAIL $acc  $fname" >&2
        rm -f "$out"
        return 1
    fi
}
export -f download_one
export ASCP KEY DEST BW_MBPS

# ── Run in parallel with xargs -P ─────────────────────────────────────────────
# Each line of $JOBS contains "<acc>\t<aspath>" → pass as two args to bash -c.
< "$JOBS" xargs -P "$PARALLEL" -I {} -d $'\n' bash -c '
    line="{}"
    acc="${line%%	*}"
    aspath="${line#*	}"
    download_one "$acc" "$aspath"
'

# ── Final report ──────────────────────────────────────────────────────────────
GOT=$(find "$DEST" -name "*.fastq.gz" -size +0 | wc -l | tr -d ' ')
echo
echo "Done. Files on disk: $GOT / $TOTAL"
