#!/bin/bash
#SBATCH --job-name=kpsc_kleb_exp
#SBATCH --output=/home/lshlt19/CNVRock/logs/kleb_exp_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/kleb_exp_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --array=1-5000%100   # overridden at submission time

# Run Kleborate v3 on Phase 2 expansion assemblies (one task per assembly).
#
# Input:  assets/kpsc_expansion_assembly_urls.tsv (col 1 = sample_accession)
#         data/assemblies/{sample_accession}.fa.gz  (from download_assemblies_s3.sh)
# Output: data/raw/kleborate_expansion/{sample_accession}.tsv  (one row per sample)
#
# After all tasks complete, merge with:
#   python3 data/setup/merge_kleborate_expansion.py
#
# Submit (after Stream B assemblies are downloaded):
#   N=$(( $(wc -l < assets/kpsc_expansion_assembly_urls.tsv) - 1 ))
#   sbatch --array=1-${N}%100 hpc/run_kleborate_expansion.sh
#
# Or via wave_submit.sh (handles MaxJobCount=10000 limit automatically).

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
export PATH="/home/lshlt19/miniconda3/envs/cnvrock/bin:/home/lshlt19/miniconda3/bin:$PATH"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cnvrock 2>/dev/null || true

URL_TSV="$REPO_DIR/assets/kpsc_expansion_assembly_urls.tsv"
ASM_DIR="$REPO_DIR/data/assemblies"
OUT_DIR="$REPO_DIR/data/raw/kleborate_expansion"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

# Support BATCH_OFFSET for chunked/wave submissions
LINE=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID + 1 ))  # +1 to skip header
ROW=$(sed -n "${LINE}p" "$URL_TSV")

if [[ -z "$ROW" ]]; then
    echo "No row for task $SLURM_ARRAY_TASK_ID (line $LINE) — skipping."
    exit 0
fi

SAMPLE_ACC=$(echo "$ROW" | awk -F'\t' '{print $1}')
ASM_GZ="$ASM_DIR/${SAMPLE_ACC}.fa.gz"
OUT="$OUT_DIR/${SAMPLE_ACC}.tsv"

# Skip if already done
if [[ -f "$OUT" ]]; then
    echo "$SAMPLE_ACC already done — skipping."
    exit 0
fi

# Check assembly exists
if [[ ! -f "$ASM_GZ" ]]; then
    echo "ERROR: assembly not found for $SAMPLE_ACC: $ASM_GZ" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ACC"

# Decompress to tmp
TMP_FA="$(mktemp /tmp/${SAMPLE_ACC}_XXXXXX.fa)"
trap "rm -f $TMP_FA" EXIT
gunzip -c "$ASM_GZ" > "$TMP_FA"

# Detect Kleborate version and run accordingly
KLEB_VER=$(kleborate --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+' | head -1 | cut -d. -f1)
if [[ "$KLEB_VER" == "3" ]]; then
    kleborate \
        --assemblies "$TMP_FA" \
        --preset     kpsc \
        --output     "${OUT}.tmp" \
        --threads    "$SLURM_CPUS_PER_TASK" 2>/dev/null
else
    # v2 syntax
    kleborate \
        -a           "$TMP_FA" \
        --resistance \
        -o           "${OUT}.tmp"
fi

# Annotate with sample_accession (Kleborate uses the filename as strain name)
# Replace the filename-derived strain name with the canonical sample accession
python3 - "$SAMPLE_ACC" "${OUT}.tmp" "$OUT" << 'PYEOF'
import sys, csv
sample_acc, inp, outp = sys.argv[1], sys.argv[2], sys.argv[3]
with open(inp) as fin, open(outp, 'w', newline='') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    fields = list(reader.fieldnames)
    writer = csv.DictWriter(fout, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for row in reader:
        row[fields[0]] = sample_acc  # first column = strain/sample name
        writer.writerow(row)
PYEOF

rm -f "${OUT}.tmp"
echo "Done: $SAMPLE_ACC → $OUT"
