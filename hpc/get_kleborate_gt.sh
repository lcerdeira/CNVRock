#!/bin/bash
#SBATCH --job-name=kpsc_kleborate_gt
#SBATCH --output=/home/lshlt19/CNVRock/logs/kleborate_gt_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/kleborate_gt_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Build Kleborate-based ground truth for ompK35, ompK36, ramR.
#
# Prerequisites:
#   1. Assemblies available in $ASSEMBLIES_DIR (one FASTA per sample, named
#      <sample_id>.fasta).  AllTheBacteria assembles can be downloaded with:
#        wget https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/<accession>.fasta.gz
#   2. kleborate installed in the cnvrock conda env:
#        conda install -n cnvrock -c conda-forge -c bioconda kleborate -y
#      OR: pip install kleborate  (v3 only)
#   3. BLAST+ installed (for ramR):
#        conda install -n cnvrock -c bioconda blast -y
#
# Usage:
#   sbatch hpc/get_kleborate_gt.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cnvrock

module load samtools/1.20

ASSEMBLIES_DIR="$REPO_DIR/data/assemblies"
KLEBORATE_OUT="$REPO_DIR/assets/kleborate_out.tsv"
AMR_GT="$REPO_DIR/assets/kpsc_amrfinder_ground_truth.tsv"
REFERENCE="$REPO_DIR/assets/HS11286.fasta"
OUT="$REPO_DIR/assets/kpsc_kleborate_ground_truth.tsv"

# ── Step 1: Run Kleborate on all assemblies ───────────────────────────────
echo "Running Kleborate on assemblies in $ASSEMBLIES_DIR …"
N_FASTA=$(find "$ASSEMBLIES_DIR" -maxdepth 1 -name "*.fasta" -o -name "*.fa" -o -name "*.fna" | wc -l)
echo "  Found $N_FASTA assemblies."

if [[ "$N_FASTA" -eq 0 ]]; then
    echo "ERROR: no assemblies found in $ASSEMBLIES_DIR"
    echo "  Download with:"
    echo "    while IFS= read -r acc; do"
    echo "      wget -q https://ftp.ebi.ac.uk/.../\${acc}.fasta.gz -O $ASSEMBLIES_DIR/\${acc}.fasta.gz && gunzip \${acc}.fasta.gz"
    echo "    done < assets/kpsc_sra_accessions.txt"
    exit 1
fi

# Detect Kleborate version and run accordingly
if kleborate --version 2>&1 | grep -q "^3\|Kleborate 3"; then
    echo "  Using Kleborate v3 (--preset kpsc)"
    kleborate \
        --assemblies "$ASSEMBLIES_DIR"/*.fasta "$ASSEMBLIES_DIR"/*.fa "$ASSEMBLIES_DIR"/*.fna \
        --preset kpsc \
        --output "$KLEBORATE_OUT" \
        --threads "$SLURM_CPUS_PER_TASK" 2>/dev/null || \
    kleborate \
        -a "$ASSEMBLIES_DIR"/*.fasta \
        --preset kpsc \
        -o "$KLEBORATE_OUT" \
        --threads "$SLURM_CPUS_PER_TASK"
else
    echo "  Using Kleborate v2 (--resistance)"
    kleborate \
        -a "$ASSEMBLIES_DIR"/*.fasta "$ASSEMBLIES_DIR"/*.fa "$ASSEMBLIES_DIR"/*.fna \
        --resistance \
        -o "$KLEBORATE_OUT" 2>/dev/null || \
    kleborate \
        -a "$ASSEMBLIES_DIR"/*.fasta \
        --resistance \
        -o "$KLEBORATE_OUT"
fi

echo "  Kleborate done → $KLEBORATE_OUT"

# ── Step 2: Parse Kleborate output + BLAST ramR ──────────────────────────
echo ""
echo "Parsing Kleborate TSV and calling ramR via BLAST …"

BLAST_ARGS=""
if command -v blastn &>/dev/null; then
    BLAST_ARGS="--reference $REFERENCE --assemblies-dir $ASSEMBLIES_DIR"
else
    echo "  WARNING: blastn not found — ramR will be -1 (uncallable). Install BLAST+:"
    echo "    conda install -n cnvrock -c bioconda blast -y"
fi

AMR_ARG=""
if [[ -f "$AMR_GT" ]]; then
    AMR_ARG="--amrfinder-gt $AMR_GT"
fi

python3 "$REPO_DIR/data/setup/get_kleborate_gt.py" \
    --kleborate-tsv "$KLEBORATE_OUT" \
    $AMR_ARG \
    $BLAST_ARGS \
    --out "$OUT"

echo ""
echo "Done. Output → $OUT"
echo ""
echo "Next steps:"
echo "  Update kpsc_kleborate_gt_path in models/experiments/22/config.yaml"
echo "  Re-run wrap-up on experiment 22 to evaluate with the new GT:"
echo "    cd ~/CNVRock && python -m models.training.wrap_up models/experiments/22/config.yaml"
