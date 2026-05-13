#!/bin/bash
#SBATCH --job-name=kpsc_aspera
#SBATCH --output=/home/lshlt19/CNVRock/logs/aspera_%A_%a.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/aspera_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# All-in-one Aspera download → BWA align → GATK CollectReadCounts pipeline
# for the KpSC expansion subset. Replaces the wget-based
# download_expansion_sra.sh, which kept failing under EBI's per-IP rate
# limiting.
#
# Per task:
#   1. Read one row from assets/kpsc_expansion_subset_5k.tsv
#   2. ascp-download R1 (and R2 if PAIRED) into /tmp scratch
#   3. BWA mem → samtools sort → BAM in scratch
#   4. samtools index
#   5. GATK CollectReadCounts at 1 kb across chromosome + plasmids
#   6. Move count file to data/raw/readcounts_subset/{ACC}.counts.tsv
#   7. Delete scratch (FASTQs + BAM)
#
# Submit:
#   N=$(( $(wc -l < assets/kpsc_expansion_subset_5k.tsv) - 1 ))
#   sbatch --array=1-${N}%50 hpc/aspera_subset_pipeline.sh

set -euo pipefail

REPO_DIR=/home/lshlt19/CNVRock
ASPERA_KEY=$HOME/.aspera/sdk/ebi_aspera_key.openssh
ASCP=$HOME/.aspera/sdk/ascp

export PATH="$HOME/miniconda3/envs/cnvrock/bin:$HOME/miniconda3/bin:$PATH"
module load bwa/0.718 samtools/1.20 gatk/4.6.0.0 java/20.0.1

MANIFEST="$REPO_DIR/assets/kpsc_expansion_subset_5k.tsv"
REFERENCE="$REPO_DIR/assets/HS11286_extended.fasta"
INTERVALS="$REPO_DIR/assets/HS11286_extended_1kb.interval_list"
COUNTS_DIR="$REPO_DIR/data/raw/readcounts_subset"

mkdir -p "$COUNTS_DIR" "$REPO_DIR/logs"

# Read one row of the subset manifest (cols: accession, layout, r1_url, r2_url)
LINE_NUM=$(( ${BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID + 1 ))   # +1 to skip header
ROW=$(sed -n "${LINE_NUM}p" "$MANIFEST")
if [[ -z "$ROW" ]]; then
    echo "No row for task $SLURM_ARRAY_TASK_ID (line $LINE_NUM) — skipping."
    exit 0
fi

ACC=$(echo    "$ROW" | cut -f1)
LAYOUT=$(echo "$ROW" | cut -f2)
R1_URL=$(echo "$ROW" | cut -f3)
R2_URL=$(echo "$ROW" | cut -f4)

OUT="$COUNTS_DIR/${ACC}.counts.tsv"
if [[ -f "$OUT" ]]; then
    echo "$ACC already done — skipping."
    exit 0
fi
if [[ -z "$R1_URL" ]]; then
    echo "ERROR: no R1 URL for $ACC" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID: $ACC (layout=$LAYOUT) on $(hostname)"

SCRATCH=/tmp/${ACC}_aspera_$$
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

R1="$SCRATCH/${ACC}_1.fastq.gz"
R2="$SCRATCH/${ACC}_2.fastq.gz"
BAM="$SCRATCH/${ACC}.bam"

# ── Build Aspera paths from https URLs ────────────────────────────────────────
# https://ftp.sra.ebi.ac.uk/<path>  →  era-fasp@fasp.sra.ebi.ac.uk:<path>
ascp_path() {
    echo "${1#https://ftp.sra.ebi.ac.uk/}"
}

ascp_fetch() {
    local URL="$1" DST="$2"
    local PATH_REL
    PATH_REL=$(ascp_path "$URL")
    local ATTEMPT=0
    while (( ATTEMPT < 3 )); do
        if "$ASCP" -P 33001 -QT -l 200m -i "$ASPERA_KEY" \
                  "era-fasp@fasp.sra.ebi.ac.uk:${PATH_REL}" "$DST" \
                  > "$SCRATCH/ascp_${PATH_REL##*/}.log" 2>&1; then
            [[ -s "$DST" ]] && return 0
        fi
        ATTEMPT=$(( ATTEMPT + 1 ))
        echo "  ascp attempt $ATTEMPT failed for $URL; retrying in 30s …"
        rm -f "$DST"
        sleep 30
    done
    echo "ERROR: ascp failed for $URL after 3 attempts" >&2
    tail -5 "$SCRATCH/ascp_${PATH_REL##*/}.log" >&2 || true
    return 1
}

# ── Download R1 (and R2 if paired) ────────────────────────────────────────────
ascp_fetch "$R1_URL" "$R1"
if [[ "$LAYOUT" == "PAIRED" && -n "$R2_URL" ]]; then
    ascp_fetch "$R2_URL" "$R2"
fi

# ── Align ─────────────────────────────────────────────────────────────────────
if [[ "$LAYOUT" == "PAIRED" && -s "$R2" ]]; then
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" "$R2" \
      | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
else
    echo "  Aligning single-end reads for $ACC"
    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" \
      | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM"
fi
samtools index "$BAM"
rm -f "$R1" "$R2"   # free FASTQ space ASAP

# ── GATK CollectReadCounts ────────────────────────────────────────────────────
gatk CollectReadCounts \
    --java-options            "-Xmx8g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality 40 \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --tmp-dir                 "$SCRATCH" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"
echo "Done: $ACC → $OUT"
