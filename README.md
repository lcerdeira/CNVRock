[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478464.svg)](https://doi.org/10.5281/zenodo.19478464)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-blue.svg)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/cnvrock/badge/?version=latest)](https://cnvrock.readthedocs.io/en/latest/?badge=latest)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://cnvrock.streamlit.app/)

[![CNVRock](./assets/cnvrock-logo.png)](./assets/cnvrock-logo.png)

# CNVRock — AMR gene copy-number variation in *Klebsiella pneumoniae* using variational autoencoders

Antimicrobial resistance (AMR) is one of the leading global health threats. In *Klebsiella pneumoniae* species complex (KpSC), resistance is driven not just by gene **presence** but by **how many copies** of a resistance gene a bacterium carries — chromosomal amplification (e.g. `blaSHV`) and plasmid copy-number variation (PCN) of clinically critical genes (`blaKPC`, `blaNDM`, `blaCTX-M`, `blaOXA-48`-like). Both effects matter clinically; conventional presence-callers see neither.

CNVRock detects both from short-read sequencing using:

- a **1D convolutional VAE** that learns a 10-dim latent embedding of genome-wide 1 kb-binned read depth,
- a **6-state Gaussian HMM** that segments the per-sample reconstructions into copy-number states,
- separate **per-gene CNV callers** — CRR-based for chromosomal genes, PCN-based for plasmid-borne AMR genes (aggregated by gene family so multi-mapped reads between near-identical allele variants are correctly counted).

## Scaling study (manuscript headline)

We report performance across **five training-set sizes** holding architecture, HMM, and CNV-caller parameters constant — only n varies:

| Tier | n samples | Manifest |
|---|---|---|
| Phase 1 baseline | 545 | (legacy) |
| **exp 32 / 5K** | 5,000 | `assets/kpsc_expansion_subset_5k.tsv` |
| **exp 33 / 10K** | 10,000 | `assets/kpsc_expansion_subset_10k.tsv` |
| **exp 34 / 20K** | 20,000 | `assets/kpsc_expansion_subset_20k.tsv` |
| **exp 36 / 40K** | 40,000 | `assets/kpsc_expansion_subset_40k.tsv` |

Each manifest is a **strict superset** of the previous one (see `data/setup/select_expansion_subset.py`) — any observed change in metrics is attributable to *more training data*, not to a different sample mix.

Results populate as experiments complete. See [`docs/08_scaling_study.md`](docs/08_scaling_study.md) for the live table.

## Pipeline at a glance

```
ENA (~78,000 KpSC FASTQ pairs)
    │
    │  IBM Aspera (ascp via bioconda aspera-cli, concur=10, MQ=0 GATK output)
    ▼
data/raw/fastq_subset/     ── BWA-MEM ── data/raw/bam_subset/  (kept)
                                                │
                                                ▼ GATK CollectReadCounts
                              data/raw/readcounts_subset_mq0/  (1 kb bins, chrom + plasmids)
                                                │
                          ┌─────────────────────┴─────────────────────┐
                          ▼                                            ▼
       data/inputs/KpSC-expansion-{N}-mq0-1000bp-npy/    KpSC-expansion-{N}-mq0-plasmid-1000bp-npy/
       (NC_016845.1 only, 5,334 bins)                    (7 gene families: blaKPC, blaCTX-M,
                          │                                            blaNDM, blaOXA-48-like,
                          │                                            blaTEM, qnrB, aac6-Ib-cr)
                          ▼                                            ▼
                ┌──────────────────────────────────────────────────────────┐
                │   VAE → Gaussian HMM → CNV callers → evaluation.txt      │
                │                  (models/train.py)                       │
                └──────────────────────────────────────────────────────────┘
```

## Quickstart

```bash
# One-time setup (HPC):
conda env create -f environment.yml && conda activate cnvrock
conda install -y -c bioconda aspera-cli
ascli config ascp install
curl -sL https://www.ebi.ac.uk/bioimage-archive/help-download/ \
    | sed -n '/BEGIN OPENSSH PRIVATE KEY/,/END OPENSSH PRIVATE KEY/p' \
    > ~/.aspera/sdk/ebi_aspera_key.openssh
chmod 600 ~/.aspera/sdk/ebi_aspera_key.openssh

# Generate subset manifests (deterministic, seed=42):
python3 data/setup/select_expansion_subset.py

# Build 1 kb interval list across chrom + plasmids:
gatk CreateSequenceDictionary -R assets/HS11286_extended.fasta
gatk PreprocessIntervals \
    --reference assets/HS11286_extended.fasta \
    --bin-length 1000 --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output assets/HS11286_extended_1kb.interval_list

# Download + align + count one tier (idempotent; samples already done are skipped):
sbatch --export=ALL,MANIFEST=assets/kpsc_expansion_subset_20k.tsv \
       --array=1-4999%10 hpc/aspera_subset_pipeline.sh
# (Plus three more chunks of size ≤5000 — SLURM MaxArraySize limit.)

# Build NPY stores (per tier):
python3 data/setup/readcounts_to_npy_kpsc.py \
    --counts-dir data/raw/readcounts_subset_mq0 \
    --manifest   assets/kpsc_expansion_subset_5k.tsv \
    --out-dir    data/inputs/KpSC-expansion-5k-mq0-1000bp-npy \
    --keep-contigs NC_016845.1

python3 data/setup/plasmid_genes_to_npy_kpsc.py \
    --counts-dir  data/raw/readcounts_subset_mq0 \
    --manifest    assets/kpsc_expansion_subset_5k.tsv \
    --gene-coords assets/plasmid_refs/plasmid_gene_coords.tsv \
    --families    assets/plasmid_refs/plasmid_gene_families.tsv \
    --out-dir     data/inputs/KpSC-expansion-5k-mq0-plasmid-1000bp-npy

# Train one experiment:
sbatch hpc/train_gpu.sh experiments/32/config.yaml
```

Full per-step explanation lives in the **[ReadTheDocs site](https://cnvrock.readthedocs.io/en/latest/)**.

## Repository layout

```
assets/                    sample manifests, ground truth, reference + intervals
  HS11286_extended.fasta              chromosome + 11 plasmid contigs
  HS11286_extended_1kb.interval_list  7,365 1 kb bins for GATK CollectReadCounts
  kpsc_expansion_subset_{5k,10k,20k,40k,80k}.tsv  nested manifests
  kpsc_expansion_kleborate_gt_runlevel.tsv         Kleborate v3 ground truth
  amrfinder_gt_expansion.tsv                       AMRFinder+ ground truth
  plasmid_refs/
    plasmid_gene_coords.tsv     per-gene CDS coordinates on the extended ref
    plasmid_gene_families.tsv   gene-family groupings (aggregation key)

data/
  raw/
    fastq_subset/                  persistent FASTQs (Aspera downloads)
    bam_subset/                    persistent BWA BAMs
    readcounts_subset_mq0/         GATK 1 kb count TSVs
  inputs/
    KpSC-expansion-{N}-mq0-1000bp-npy/                chromosome NPY store
    KpSC-expansion-{N}-mq0-plasmid-1000bp-npy/        plasmid-family NPY store
  results/
    {32,33,34,36}_kpsc_expansion_{5k,10k,20k,40k}/   per-tier CNVRock outputs

hpc/
  aspera_subset_pipeline.sh        ascp + BWA + GATK per SLURM array task
  train_gpu.sh                     submits one CNVRock training run on an A40

data/setup/
  select_expansion_subset.py       stratified nested 5K/10K/20K/40K/80K
  readcounts_to_npy_kpsc.py        TSV → chromosome NPY store
  plasmid_genes_to_npy_kpsc.py     TSV → plasmid-family NPY store
  reannotate_plasmid_genes_blast.py     BLAST-based coord verification
  reannotate_plasmid_genes_genbank.py   GenBank-based coord verification

models/
  train.py                         single entry point: VAE + HMM + CNV + eval
  architecture/06_conv_vae.py      1D Conv-VAE (10-dim latent)
  hmm/02_gaussian_hmm.py           6-state Gaussian HMM segmenter
  cnv/06_gene_cnv_caller.py        chromosomal CRR-based caller
  cnv/07_plasmid_cnv_caller.py     plasmid PCN-based caller
  evaluation/04_kpsc_evaluation.py per-gene MCC / FNR / PPV
  experiments/32-36/               scaling-study experiment configs

docs/                              Sphinx + MyST documentation (publishes to RTD)
archive/experiments_pf/            historical Phase 1 experiments (exp 01–21)
```

## Key parameter choices

These are documented at length in [`docs/09_methods.md`](docs/09_methods.md); the highlights:

- **MQ filter = 0** for GATK `CollectReadCounts`. Multi-mapping between nearly-identical AMR allele variants (NDM-1 / NDM-5; OXA-48 / OXA-181; CTX-M-15 / -14 / -65 / -27) is *the signal*, not noise — filtering it out makes those genes invisible. We aggregate by gene family downstream.
- **Aspera concurrency = 10** per array job (EBI's per-source-IP cap).
- **Stratification**: species × `Bla_Carb` × ST cap of 150 per stratum, with 1.5× carbapenemase-carrier oversample weighting.
- **Nested manifest construction**: one weighted shuffle of the eligible pool, then prefixes of length 5K / 10K / 20K / 40K. Strictly nested by construction.

## Documentation

The full pipeline narrative — overview, data acquisition, reference + intervals, subset selection, NPY stores, training, evaluation, scaling-study results, methods, reproducibility — is at **[cnvrock.readthedocs.io](https://cnvrock.readthedocs.io/en/latest/)** and builds automatically on every push to `main` from the [`docs/`](docs/) folder.

## Interactive diagnostics

A live Streamlit app showing the monitor / sample-viewer / latent-coverage pages on a 200-sample subsample of exp 32 is hosted at **<https://cnvrock.streamlit.app/>**. No install required for reviewers — just click and explore. The full pipeline (training + all 5 scaling tiers) reproduces locally via the [Quickstart](#quickstart) above; the app then auto-detects the full `data/results/` and switches out of demo mode.

## License

MIT. See [`LICENSE`](LICENSE).
