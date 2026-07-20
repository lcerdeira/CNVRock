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
#   1. Read one row from the manifest TSV (override via MANIFEST=... env var,
#      default is assets/kpsc_expansion_subset_5k.tsv)
#   2. ascp-download R1 (and R2 if PAIRED) into data/raw/fastq_subset/
#   3. BWA mem → samtools sort → BAM in data/raw/bam_subset/
#   4. samtools index
#   5. GATK CollectReadCounts at 1 kb across chromosome + plasmids,
#      with MIN_MQ=0 (keeps multi-mapped reads — required to detect plasmid
#      AMR genes whose nearly-identical alleles live on multiple reference
#      contigs; reads are then aggregated by gene family downstream)
#   6. Move count file to data/raw/readcounts_subset_mq0/{ACC}.counts.tsv
#
# Idempotent: skips download/align if FASTQs/BAM already exist on disk, and
# the script's early `[[ -f $OUT ]] && exit 0` check skips already-counted
# samples entirely. FASTQs + BAMs are KEPT (HPC has plenty of space per Alex
# 2026-05-15) so we can re-extract counts at different MQ thresholds without
# re-downloading. Only /tmp scratch (ascp logs, GATK tmp-dir) is removed at
# task exit.
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

MANIFEST="${MANIFEST:-$REPO_DIR/assets/kpsc_expansion_subset_5k.tsv}"
# Reference / output dirs are env-overridable so the same pipeline serves
# other organisms (e.g. A. baumannii). Defaults = KpSC.
REFERENCE="${REFERENCE:-$REPO_DIR/assets/HS11286_extended.fasta}"
INTERVALS="${INTERVALS:-$REPO_DIR/assets/HS11286_extended_1kb.interval_list}"
COUNTS_DIR="${COUNTS_DIR:-$REPO_DIR/data/raw/readcounts_subset_mq0}"
FASTQ_DIR="${FASTQ_DIR:-$REPO_DIR/data/raw/fastq_subset}"
BAM_DIR="${BAM_DIR:-$REPO_DIR/data/raw/bam_subset}"

# MQ filter: 0 (was 40 in Phase 1, then 20 in our first re-run).
#
# Diagnosis 2026-05-17 (GenBank + BLAST cross-check): the
# plasmid_gene_coords.tsv coordinates are correct. The true cause of zero
# plasmid-gene counts is *intentional multi-mapping*: nearly-identical AMR
# allele variants (NDM-1 vs NDM-5, OXA-48 vs OXA-181, CTX-M-15 vs
# CTX-M-65 vs CTX-M-27) live on different contigs in HS11286_extended.fasta,
# and reads from a sample carrying ANY one of them multi-map equally to all,
# get MQ=0 from BWA, and are filtered out at any MQ threshold above 0.
#
# At MQ=0 these multi-mapped reads are kept (BWA places them randomly at
# one of the homologous positions). Aggregating bin counts across ALL
# allele-variant positions per gene family (see plasmid_genes_to_npy_kpsc.py
# --families option and assets/plasmid_refs/plasmid_gene_families.tsv)
# recovers the full per-sample per-family signal.
#
# For chromosomal CNV detection (blaSHV): the KpSC chromosome is largely
# uniquely-mappable, so MQ=0 vs MQ=20 affects only repeat / transposon
# regions which the VAE treats as noise. Single pipeline, single count
# file per sample.
# Env-overridable: A. baumannii targets (AdeABC/AdeIJK, blaOXA) are
# chromosomal and uniquely-mappable, so that run sets MIN_MQ=20.
MIN_MQ="${MIN_MQ:-0}"

mkdir -p "$COUNTS_DIR" "$FASTQ_DIR" "$BAM_DIR" "$REPO_DIR/logs"

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
# Note: per user request (2026-05-17), we now KEEP FASTQs and BAMs persistently
# under data/raw/fastq_subset/ and data/raw/bam_subset/ — only /tmp scratch
# (used for ascp logs + GATK tmp-dir) is cleaned up at task exit.
trap "rm -rf $SCRATCH" EXIT

R1="$FASTQ_DIR/${ACC}_1.fastq.gz"
R2="$FASTQ_DIR/${ACC}_2.fastq.gz"
BAM="$BAM_DIR/${ACC}.bam"

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

# ── Download R1 (and R2 if paired) — skip if already present ─────────────────
if [[ ! -s "$R1" ]]; then
    ascp_fetch "$R1_URL" "$R1"
fi
if [[ "$LAYOUT" == "PAIRED" && -n "$R2_URL" && ! -s "$R2" ]]; then
    ascp_fetch "$R2_URL" "$R2"
fi

# ── Align — skip if BAM already present ──────────────────────────────────────
if [[ ! -s "$BAM" || ! -s "${BAM}.bai" ]]; then
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
fi

# ── GATK CollectReadCounts at MQ=10 ───────────────────────────────────────────
gatk CollectReadCounts \
    --java-options            "-Xmx8g" \
    --reference               "$REFERENCE" \
    --intervals               "$INTERVALS" \
    --input                   "$BAM" \
    --format                  TSV \
    --read-filter             MappingQualityReadFilter \
    --minimum-mapping-quality "$MIN_MQ" \
    --interval-merging-rule   OVERLAPPING_ONLY \
    --tmp-dir                 "$SCRATCH" \
    --output                  "${OUT}.tmp"

mv "${OUT}.tmp" "$OUT"

# ── Optional FASTQ cleanup ───────────────────────────────────────────────────
# By default FASTQs are kept (per user request 2026-05-17). For large
# multi-organism runs where storage is the binding constraint, set
# KEEP_FASTQ=0 to drop the raw reads once the BAM and counts exist. The BAM is
# retained (needed for the MQ>=20 chromosomal second pass), and the FASTQ is
# re-downloadable from ENA, so this is a safe, reversible space saving.
KEEP_FASTQ="${KEEP_FASTQ:-1}"
if [[ "$KEEP_FASTQ" == "0" && -s "$OUT" && -s "$BAM" ]]; then
    rm -f "$R1" "$R2" "$FASTQ_DIR/${ACC}.fastq.gz"
    echo "Done: $ACC → $OUT (BAM kept, FASTQ removed)"
else
    echo "Done: $ACC → $OUT (FASTQs + BAM kept)"
fi
