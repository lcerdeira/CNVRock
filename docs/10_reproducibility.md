# 10. Reproducibility

This page is the **one-shot recipe** to regenerate the manuscript from a
fresh clone.

## Software

- Python 3.11
- PyTorch 2.5.1 (CUDA 12.1)
- pandas, numpy, scikit-learn, hmmlearn, tqdm, PyYAML
- BWA 0.7.18, samtools 1.20, GATK 4.6.0.0 (HPC modules)
- IBM Aspera SDK 4.4.7+ (`ascp`)

`environment.yml` at the repo root pins exact versions.

## HPC environment

LSHTM compute cluster (Slurm), partition `normal` for CPU tasks and the GPU
partition for training. Conda environment name: `cnvrock`.

```bash
# One-time setup (HPC):
conda env create -f environment.yml
conda activate cnvrock
conda install -y -c bioconda aspera-cli
ascli config ascp install   # installs ~/.aspera/sdk/ascp

# EBI Aspera SSH key (one-time):
curl -sL https://www.ebi.ac.uk/bioimage-archive/help-download/ \
    | sed -n '/BEGIN OPENSSH PRIVATE KEY/,/END OPENSSH PRIVATE KEY/p' \
    > ~/.aspera/sdk/ebi_aspera_key.openssh
chmod 600 ~/.aspera/sdk/ebi_aspera_key.openssh
```

## End-to-end recipe

```bash
# ── 1. Subset manifests ───────────────────────────────────────────────────────
python3 data/setup/select_expansion_subset.py
# Writes assets/kpsc_expansion_subset_{5k,10k,20k,40k,80k}.tsv

# ── 2. Reference + 1 kb intervals ────────────────────────────────────────────
gatk CreateSequenceDictionary -R assets/HS11286_extended.fasta
gatk PreprocessIntervals \
    --reference assets/HS11286_extended.fasta \
    --bin-length 1000 \
    --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output assets/HS11286_extended_1kb.interval_list

# ── 3. Ground truth files (BioSample → run accession bridge) ─────────────────
python3 - <<'EOF'
import pandas as pd
kleb = pd.read_csv('assets/kpsc_expansion_kleborate_gt.tsv',
                   sep='\t', dtype=str, low_memory=False)
meta = pd.read_csv('assets/kpsc_expansion_metadata.tsv', sep='\t', dtype=str)
meta['run_acc'] = meta['sample_id'].str.split(',').str[0]

bridged = kleb.merge(meta[['sample_accession', 'run_acc']],
                     left_on='strain', right_on='sample_accession',
                     how='inner')
bridged = bridged.drop(columns=['strain', 'sample_accession']) \
                 .rename(columns={'run_acc': 'sample_id'})
bridged.to_csv('assets/kpsc_expansion_kleborate_gt_runlevel.tsv',
               sep='\t', index=False)

run_meta = meta.copy()
run_meta['sample_id'] = run_meta['sample_id'].str.split(',')
run_meta = run_meta.explode('sample_id').reset_index(drop=True)
run_meta = run_meta.rename(columns={'species': 'Species'})
run_meta = run_meta.merge(
    kleb[['strain', 'ST']].rename(columns={'strain': 'sample_accession'}),
    on='sample_accession', how='left'
)
run_meta[['sample_id', 'Species', 'ST', 'sample_accession']].to_csv(
    'assets/kpsc_expansion_metadata_runlevel.tsv', sep='\t', index=False
)
EOF

# ── 4. Download + align + count (per tier; idempotent) ───────────────────────
# Each tier reuses count files from the previous tier — only run the largest.
for chunk in 1 2 3 4; do
    case $chunk in
        1) OFFSET=0;     ARR=1-4999 ;;
        2) OFFSET=4999;  ARR=1-4999 ;;
        3) OFFSET=9998;  ARR=1-4999 ;;
        4) OFFSET=14997; ARR=1-5000 ;;
    esac
    sbatch --parsable \
        --export=ALL,MANIFEST=$PWD/assets/kpsc_expansion_subset_20k.tsv,BATCH_OFFSET=$OFFSET \
        --array=$ARR%10 hpc/aspera_subset_pipeline.sh
done
# For 40k: scale to 8 chunks of ≤ 5000 each.

# ── 5. NPY stores (per tier) ─────────────────────────────────────────────────
for N in 5k 10k 20k 40k; do
    python3 data/setup/readcounts_to_npy_kpsc.py \
        --counts-dir data/raw/readcounts_subset_mq0 \
        --manifest   assets/kpsc_expansion_subset_${N}.tsv \
        --out-dir    data/inputs/KpSC-expansion-${N}-mq0-1000bp-npy \
        --keep-contigs NC_016845.1

    python3 data/setup/plasmid_genes_to_npy_kpsc.py \
        --counts-dir  data/raw/readcounts_subset_mq0 \
        --manifest    assets/kpsc_expansion_subset_${N}.tsv \
        --gene-coords assets/plasmid_refs/plasmid_gene_coords.tsv \
        --families    assets/plasmid_refs/plasmid_gene_families.tsv \
        --out-dir     data/inputs/KpSC-expansion-${N}-mq0-plasmid-1000bp-npy
done

# ── 6. Train + evaluate ──────────────────────────────────────────────────────
sbatch hpc/train_gpu.sh experiments/32/config.yaml   # 5K
sbatch hpc/train_gpu.sh experiments/33/config.yaml   # 10K
sbatch hpc/train_gpu.sh experiments/34/config.yaml   # 20K
sbatch hpc/train_gpu.sh experiments/36/config.yaml   # 40K
```

## Commit hashes

Each scaling-study experiment's `evaluation.txt` is committed at a known
hash. To exactly reproduce a result, check out the repo at that hash and
run the recipe.

| Experiment | Final commit | Date |
|---|---|---|
| exp 32 (5K, MQ=20) | (pending) | (pending) |
| exp 33 (10K, MQ=20) | (pending) | (pending) |
| exp 34 (20K, MQ=20) | (pending) | (pending) |
| exp 36 (40K, MQ=20) | (pending) | (pending) |

## Random seeds

- Subset selection: `numpy.random.default_rng(42)` → identical manifests
  every time.
- VAE training: PyTorch default (non-deterministic on GPU). Run-to-run
  variability is small; for the manuscript we report a single run per tier
  and discuss variability in the supplementary.

## Storage

Per 20K samples:

| Asset | Size (≈) |
|---|---|
| FASTQs (paired, gzipped) | 6–10 TB |
| BAMs (sorted + indexed) | 18–24 TB |
| Count TSVs (7,365 bins) | 1.5 GB |
| NPY stores (chrom + plasmid) | 1.5 GB |
| Trained checkpoints + outputs | 200 MB / exp |

FASTQs and BAMs are kept on the HPC scratch until the 40K experiment is
complete; then archived or deleted depending on storage pressure. Count
files and NPY stores are kept indefinitely.
