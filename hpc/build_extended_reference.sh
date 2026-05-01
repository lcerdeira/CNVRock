#!/bin/bash
#SBATCH --job-name=kpsc_extend_ref
#SBATCH --output=/home/lshlt19/CNVRock/logs/extend_ref_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/extend_ref_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Build extended reference and BWA index.
#
# Run AFTER get_plasmid_references.py has produced:
#   assets/HS11286_extended.fasta
#   assets/plasmid_refs/plasmid_gene_coords.tsv
#
# Usage:
#   sbatch hpc/build_extended_reference.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718
module load samtools/1.20

EXTENDED="$REPO_DIR/assets/HS11286_extended.fasta"

if [[ ! -f "$EXTENDED" ]]; then
    echo "ERROR: $EXTENDED not found. Run get_plasmid_references.py first."
    exit 1
fi

# ── samtools faidx ──────────────────────────────────────────────────────────
echo "Indexing FASTA …"
samtools faidx "$EXTENDED"

# ── BWA index ───────────────────────────────────────────────────────────────
echo "Building BWA index (this takes ~5 min for a ~6 Mb genome) …"
bwa index "$EXTENDED"

# ── GATK sequence dictionary ────────────────────────────────────────────────
module load gatk/4.6.0.0 java/20.0.1
DICT="${EXTENDED%.fasta}.dict"
if [[ ! -f "$DICT" ]]; then
    echo "Building sequence dictionary …"
    gatk CreateSequenceDictionary \
        --java-options "-Xmx6g" \
        -R "$EXTENDED" \
        -O "$DICT"
fi

# ── Build plasmid interval list for GATK CollectReadCounts ──────────────────
# Reads plasmid_gene_coords.tsv, adds ±FLANK bp around each gene,
# and writes assets/kpsc-plasmid-genes.interval_list.

COORDS="$REPO_DIR/assets/plasmid_refs/plasmid_gene_coords.tsv"
INTERVAL_LIST="$REPO_DIR/assets/kpsc-plasmid-genes.interval_list"
FLANK=10000   # 10 kb each side of the gene

echo "Building plasmid interval list …"

# Write SAM-format header from the extended reference dict
grep "^@" "$DICT" | grep -v "^@HD" > "$INTERVAL_LIST.tmp" || true
echo -e "@HD\tVN:1.6" | cat - "$INTERVAL_LIST.tmp" > "$INTERVAL_LIST"
rm -f "$INTERVAL_LIST.tmp"

# Append one interval per plasmid gene (gene ± FLANK, clamped to contig length)
python3 - <<'PYEOF'
import csv, sys
from pathlib import Path

repo     = Path("/home/lshlt19/CNVRock")
coords   = repo / "assets/plasmid_refs/plasmid_gene_coords.tsv"
fai      = repo / "assets/HS11286_extended.fasta.fai"
out_path = repo / "assets/kpsc-plasmid-genes.interval_list"
flank    = 10000

# Parse contig lengths from .fai
contig_len = {}
with open(fai) as f:
    for line in f:
        parts = line.split("\t")
        contig_len[parts[0]] = int(parts[1])

with open(coords) as f, open(out_path, "a") as out:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        contig = row["contig"]
        if contig not in contig_len:
            print(f"WARNING: {contig} not in extended FASTA — skipping {row['gene']}", flush=True)
            continue
        start = max(1, int(row["start"]) - flank)
        end   = min(contig_len[contig], int(row["end"]) + flank)
        # GATK interval_list format: contig start end strand name
        out.write(f"{contig}\t{start}\t{end}\t+\t{row['gene']}\n")
        print(f"  {row['gene']:15s} {contig}:{start}-{end}")

print(f"Interval list → {out_path}", flush=True)
PYEOF

echo ""
echo "Done. Files produced:"
echo "  $EXTENDED"
echo "  $EXTENDED.fai"
echo "  $EXTENDED.{amb,ann,bwt,pac,sa}  (BWA index)"
echo "  $DICT"
echo "  $REPO_DIR/assets/kpsc-plasmid-genes.interval_list"
echo ""
echo "Next steps:"
echo "  sbatch hpc/collect_plasmid_readcounts.sh"
echo "  Update REFERENCE in hpc/download_sra.sh → assets/HS11286_extended.fasta"
