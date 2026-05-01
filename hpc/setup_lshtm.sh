#!/bin/bash
# HPC setup script — run once on loginhpc.lshtm.ac.uk
# Usage: bash hpc/setup_lshtm.sh
#
#   ssh pryor && ssh loginhpc
#   git clone git@github.com:lcerdeira/CNVRock.git ~/CNVRock
#   bash ~/CNVRock/hpc/setup_lshtm.sh
set -euo pipefail

REPO_DIR="$HOME/CNVRock"
CONDA_DIR="$HOME/miniconda3"
CONDA_ENV="cnvrock"

echo "=== CNVRock HPC Setup ==="
echo "User:      $USER"
echo "Host:      $(hostname)"
echo "Home:      $HOME  ($(df -h "$HOME" | awk 'NR==2{print $4}') free)"
echo "Repo:      $REPO_DIR"
echo "Conda:     $CONDA_DIR"
echo ""

# ── 1. Activate conda ─────────────────────────────────────────────────────
echo "--- Conda ---"
source "$CONDA_DIR/etc/profile.d/conda.sh"

if conda env list | grep -q "^$CONDA_ENV "; then
    echo "  Env '$CONDA_ENV' already exists."
else
    echo "  Creating env '$CONDA_ENV' with Python 3.11 ..."
    conda create -n "$CONDA_ENV" python=3.11 -c conda-forge --override-channels -y
fi

conda activate "$CONDA_ENV"
echo "  [OK]  $(python --version) at $(which python)"
echo ""

# ── 2. PyTorch with CUDA 12.x ─────────────────────────────────────────────
echo "--- PyTorch ---"
if python -c "import torch" 2>/dev/null; then
    echo "  Already installed: $(python -c 'import torch; print(torch.__version__)')"
else
    module load cuda/12.4.0 2>/dev/null || true
    pip install -q torch --index-url https://download.pytorch.org/whl/cu121
    echo "  [OK]  $(python -c 'import torch; print(torch.__version__)')"
fi
echo ""

# ── 3. Remaining requirements ─────────────────────────────────────────────
echo "--- Python packages ---"
grep -v "^torch" "$REPO_DIR/requirements.txt" | pip install -q -r /dev/stdin
echo "  [OK]  All packages installed."
echo ""

# ── 4. Data directories ───────────────────────────────────────────────────
echo "--- Data directories ---"
mkdir -p "$REPO_DIR/data/inputs" \
         "$REPO_DIR/data/results" \
         "$REPO_DIR/data/raw/fastq" \
         "$REPO_DIR/data/raw/bam" \
         "$REPO_DIR/assets" \
         "$REPO_DIR/logs"
echo "  Created under $REPO_DIR/"
echo ""

# ── 5. Activation helper ──────────────────────────────────────────────────
cat > "$HOME/activate_cnvrock.sh" <<ACTIVATE
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate $CONDA_ENV
module load gatk/4.6.0.0 bwa/0.718 samtools/1.20 nextflow/24.0.4 singularity/4.4.1 2>/dev/null || true
cd $REPO_DIR
ACTIVATE
echo "  Activation helper: source ~/activate_cnvrock.sh"
echo ""

# ── 6. Verify ─────────────────────────────────────────────────────────────
echo "--- Verification ---"
python -c "
import torch, numpy, scipy, pandas, sklearn, hmmlearn, yaml, tqdm, pyarrow
print(f'  torch     {torch.__version__}  | CUDA available on login node: {torch.cuda.is_available()} (False is normal)')
print(f'  numpy     {numpy.__version__}')
print(f'  pandas    {pandas.__version__}')
print(f'  pyarrow   {pyarrow.__version__}')
"
echo ""

# ── 7. Next steps ─────────────────────────────────────────────────────────
echo "=== Next steps ==="
echo ""
echo "1. Find SLURM partition names:"
echo "   sinfo -o '%P %a %G %l'"
echo "   Then update --partition in hpc/download_sra.sh and hpc/train_gpu.sh"
echo ""
echo "2. Download HS11286 reference:"
echo "   cd $REPO_DIR/assets"
echo "   wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz'"
echo "   gunzip *.fna.gz && mv *.fna HS11286.fasta"
echo "   module load bwa/0.718 samtools/1.20 gatk/4.6.0.0"
echo "   bwa index HS11286.fasta && samtools faidx HS11286.fasta"
echo "   gatk CreateSequenceDictionary -R HS11286.fasta"
echo ""
echo "3. Select KpSC samples:"
echo "   python data/setup/atb_sample_selection.py select \\"
echo "     --atb-metadata atb_metadata.tsv \\"
echo "     --out-accessions assets/kpsc_sra_accessions.txt \\"
echo "     --out-metadata assets/kpsc_sample_metadata.tsv"
echo ""
echo "4. Submit SRA download array job:"
echo "   sbatch hpc/download_sra.sh"
echo ""
echo "5. After BAMs ready, run GATK read counts:"
echo "   cd $REPO_DIR/models/experiments/21"
echo "   module load nextflow/24.0.4 singularity/4.4.1"
echo "   nextflow run main.nf -profile lshtm_slurm -resume"
echo ""
echo "Setup complete!"
