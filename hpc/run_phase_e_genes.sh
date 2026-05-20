#!/bin/bash
#SBATCH --job-name=phase_e_genes
#SBATCH --output=/home/lshlt19/CNVRock/logs/phase_e_genes_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/phase_e_genes_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Phase E: add 8 new resistance genes to the extended reference, build the
# reference + indices, and produce the per-sample GATK counts for the new
# contigs. The downstream NPY rebuild + retrain runs as a separate job.
set -eo pipefail   # no -u: conda hooks reference unset vars

REPO=/home/lshlt19/CNVRock
cd "$REPO"
export PATH="$HOME/miniconda3/envs/blast_env/bin:$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"

# cnvrock env has Biopython; blast_env has blastn binary. We need both: keep
# cnvrock's python3 FIRST on PATH (so `import Bio` works) and APPEND blast_env
# at the end so bare `blastn` is still discoverable by the subprocess call.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cnvrock
export PATH="$PATH:$HOME/miniconda3/envs/blast_env/bin"

echo "PATH check — python3: $(command -v python3)"
echo "PATH check — blastn:  $(command -v blastn || echo 'NOT FOUND')"
python3 -c "import Bio; print('Biopython', Bio.__version__)"

echo "[1/3] Adding Phase E gene rows + plasmid contigs…"
python3 data/setup/add_phase_e_genes.py

echo "[2/3] Rebuilding BWA index for HS11286_extended.fasta…"
module load bwa samtools gatk/4.6.0.0 java/20.0.1 2>/dev/null || true
bwa index assets/HS11286_extended.fasta
samtools faidx assets/HS11286_extended.fasta
gatk CreateSequenceDictionary -R assets/HS11286_extended.fasta || true
gatk PreprocessIntervals \
    --reference assets/HS11286_extended.fasta \
    --bin-length 1000 --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output assets/HS11286_extended_1kb.interval_list

echo "[3/3] Submitting per-sample re-mapping (BWA-MEM) + counting:"
N=$(ls data/raw/bam_subset/*.bam | wc -l)
echo "  BAMs to re-process: $N"
# Submit re-alignment + counting array (writes to data/raw/readcounts_subset_mq20
# and data/raw/readcounts_subset_mq0). The shell script auto-skips samples whose
# count files already exist on disk — but for the new contigs we need a fresh
# pass, so we redirect to phase_e-specific output dirs.
sbatch --array=1-${N}%50 hpc/realign_for_phase_e.sh

echo "Phase E setup submitted."
