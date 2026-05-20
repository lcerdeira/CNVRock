#!/bin/bash
#SBATCH --job-name=cnv_scan
#SBATCH --output=/home/lshlt19/CNVRock/logs/cnv_scan_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/cnv_scan_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Genome-wide CNV–phenotype scan (Phase E exploratory analysis)
# Output: data/results/cnv_scan_phase_e/{scan_full,scan_significant,scan_top50_per_drug}.tsv
set -euo pipefail

REPO=/home/lshlt19/CNVRock
cd "$REPO"
export PATH="$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"

# Fetch HS11286 GFF3 (chromosome annotation) if not present
GFF="$REPO/assets/HS11286.gff3"
if [[ ! -s "$GFF" ]]; then
    echo "Fetching HS11286 GFF3 from NCBI…"
    python3 - <<'PY'
from Bio import Entrez
Entrez.email = "louise.cerdeira@gmail.com"
import requests, sys
url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       "?db=nuccore&id=NC_016845.1&rettype=gff3&retmode=text")
r = requests.get(url, timeout=120); r.raise_for_status()
open("assets/HS11286.gff3", "w").write(r.text)
print("ok", len(r.text), "bytes")
PY
fi

python3 analysis/genome_wide_cnv_scan.py \
    --store-dir   "$REPO/data/inputs/KpSC-expansion-10k-mq20-1000bp-npy" \
    --results-dir "$REPO/data/results/33_kpsc_expansion_10k" \
    --meta        "$REPO/assets/kpsc_expansion_metadata_runlevel.tsv" \
    --cabbage     "$REPO/assets/cabbage_kpsc_phenotypes.tsv" \
    --gff         "$GFF" \
    --out-dir     "$REPO/data/results/cnv_scan_phase_e"
