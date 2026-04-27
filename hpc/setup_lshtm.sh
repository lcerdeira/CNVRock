#!/bin/bash
# HPC setup script — run once after first SSH to loginhpc.lshtm.ac.uk
# Usage: bash hpc/setup_lshtm.sh [scratch_root]
#
# What this does:
#   1. Clones the CNVRock repo to your HPC scratch space
#   2. Creates the Python venv with all dependencies
#   3. Checks for required HPC modules (Nextflow, Singularity, bwa)
#   4. Prints next steps
#
# Prerequisites (run these first):
#   ssh pryor              # gateway node
#   ssh loginhpc           # HPC login node
#   bash hpc/setup_lshtm.sh /scratch/lshlt19   # or your project dir
set -euo pipefail

REPO_URL="git@github.com:lcerdeira/CNVRock.git"
SCRATCH="${1:-/scratch/$USER}"
REPO_DIR="$SCRATCH/CNVRock"

echo "=== CNVRock HPC Setup ==="
echo "User:        $USER"
echo "Host:        $(hostname)"
echo "Scratch:     $SCRATCH"
echo "Repo target: $REPO_DIR"
echo ""

# ── 1. Check required modules ──────────────────────────────────────────────
echo "--- Checking HPC modules ---"
for mod in nextflow singularity bwa samtools python; do
    if command -v "$mod" &>/dev/null; then
        echo "  [OK]  $mod: $(command -v $mod)"
    else
        echo "  [!!]  $mod not found — you may need: module load $mod"
    fi
done
echo ""
echo "If modules are missing, check available modules with:"
echo "  module avail nextflow"
echo "  module avail singularity"
echo ""

# ── 2. Clone repo ─────────────────────────────────────────────────────────
if [[ -d "$REPO_DIR/.git" ]]; then
    echo "--- Repo already exists — pulling latest ---"
    git -C "$REPO_DIR" pull
else
    echo "--- Cloning CNVRock to $REPO_DIR ---"
    mkdir -p "$SCRATCH"
    git clone "$REPO_URL" "$REPO_DIR"
fi
echo ""

# ── 3. Create Python venv ─────────────────────────────────────────────────
echo "--- Setting up Python venv ---"
cd "$REPO_DIR"
if [[ ! -d ".venv" ]]; then
    python3 -m venv .venv
    .venv/bin/pip install --upgrade pip
    if [[ -f requirements.txt ]]; then
        .venv/bin/pip install -r requirements.txt
        echo "  Installed from requirements.txt"
    else
        echo "  WARNING: no requirements.txt found — install packages manually."
        echo "  Typical: pip install torch numpy pandas scikit-learn hmmlearn tqdm pyyaml pyarrow"
    fi
else
    echo "  Venv already exists — skipping."
fi
echo ""

# ── 4. Create data directories ────────────────────────────────────────────
echo "--- Creating data directories ---"
mkdir -p "$SCRATCH/CNVRock/data/inputs"
mkdir -p "$SCRATCH/CNVRock/data/results"
mkdir -p "$SCRATCH/CNVRock/assets"
echo "  $SCRATCH/CNVRock/data/inputs    ← read count NPYs go here"
echo "  $SCRATCH/CNVRock/data/results   ← experiment outputs go here"
echo "  $SCRATCH/CNVRock/assets         ← reference files, BEDs, ground truth"
echo ""

# ── 5. Print next steps ───────────────────────────────────────────────────
echo "=== Next steps ==="
echo ""
echo "1. Update experiment 21 reference paths in nextflow.config:"
echo "   $REPO_DIR/models/experiments/21/nextflow.config"
echo "   Set reference_fasta/fai/dict to paths under $SCRATCH/CNVRock/assets/"
echo ""
echo "2. Confirm SLURM partition name:"
echo "   sinfo -o '%P %a %l %F'   # list partitions"
echo "   Then update 'queue' in models/experiments/21/nextflow.config"
echo ""
echo "3. Download HS11286 reference genome:"
echo "   cd $SCRATCH/CNVRock/assets"
echo "   # Option A — NCBI datasets CLI (if available):"
echo "   datasets download genome accession GCF_000240185.1 --include genome,gff3"
echo "   # Option B — direct wget:"
echo "   wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz'"
echo "   gunzip *.fna.gz && mv *.fna HS11286.fasta"
echo "   bwa index HS11286.fasta"
echo "   samtools faidx HS11286.fasta"
echo "   gatk CreateSequenceDictionary -R HS11286.fasta"
echo ""
echo "4. Select KpSC samples (on local machine or HPC):"
echo "   .venv/bin/python data/setup/atb_sample_selection.py select \\"
echo "     --atb-metadata atb_metadata.tsv \\"
echo "     --out-accessions assets/kpsc_sra_accessions.txt \\"
echo "     --out-metadata assets/kpsc_sample_metadata.tsv \\"
echo "     --max-samples 5000"
echo ""
echo "5. Download raw reads (can be parallelised with a SLURM array job):"
echo "   sbatch hpc/download_sra.sh   # created by this repo — see hpc/ directory"
echo ""
echo "Setup complete. See models/experiments/21/README.md for full pipeline docs."
