# CNVRock — AMR gene copy-number variation in *Klebsiella pneumoniae* using variational autoencoders

Antimicrobial resistance (AMR) is one of the leading global health threats. In *Klebsiella pneumoniae* Species Complex (KpSC), resistance is driven not just by gene presence but by **how many copies** of a resistance gene a bacterium carries — amplification of chromosomal genes and variable plasmid copy number (PCN) both matter clinically.

CNVRock adapts the [autoresearch](https://github.com/karpathy/autoresearch) strategy — where an AI agent continuously proposes and runs ML experiments overnight — to detect AMR-related copy-number variation in KpSC whole-genome sequencing data. A convolutional VAE learns a low-dimensional representation of genome-wide read depth; a Gaussian HMM segments the latent trajectories into copy-number states; a gene caller converts those states into per-gene calls. Claude proposes the next experiment, emails a summary, and a background daemon runs it after authorisation.

**Results — full cohort (exp 29, n=545) and hold-out validation (exp 30, 20% stratified split, n=109):**

| Gene | Type | MCC (full) | MCC (hold-out) | FNR (hold-out) | PPV (hold-out) | Notes |
|------|------|-----------|---------------|---------------|---------------|-------|
| blaSHV | chrom amp | — | — | — | — | GT unreliable: assembly collapse hides tandem dups |
| blaKPC-2 | plasmid | 1.00 | 1.00 | 0.00 | 1.00 | |
| blaNDM-1 | plasmid | 0.99 | 1.00 | 0.00 | 1.00 | |
| blaCTX-M-15 | plasmid | 0.82 | 0.88 | 0.13 | 1.00 | Irreducible FNs: reads cross-map to chromosomal blaSHV |
| qnrB1 | plasmid | 0.98 | 0.92 | 0.12 | 1.00 | |
| blaOXA-48 | plasmid | 0.98 | 0.96 | 0.06 | 1.00 | GT covers OXA-48-like family (OXA-181/232/244) |
| aac(6')-Ib-cr | plasmid | 0.86 | 0.93 | 0.04 | 0.93 | |

Hold-out split: 20% stratified by blaKPC/blaCTX-M/blaNDM presence (seed=42). Model trained on 437 samples, evaluated on 109 unseen samples. Performance is consistent with the full-cohort results, confirming generalisation.

## Layout

```
data/
  inputs/       read-count NPY stores (chromosomal + plasmid)
  results/      per-experiment outputs (one folder per experiment)
  setup/        scripts to prepare reference, extract counts, build stores

assets/         sample manifests, BAM accession lists, reference files,
                AMRFinder+ and Kleborate ground-truth TSVs,
                plasmid reference FASTAs and gene coordinate TSV

hpc/            SLURM scripts for LSHTM HPC
  build_extended_reference.sh   BWA-index the extended reference
  remap_unmapped_to_plasmids.sh remap unmapped reads to plasmid contigs

models/
  train.py          entry point — runs a full experiment from a config
  architectures/    versioned VAE definitions  (06_conv_vae.py, …)
  hmm/              versioned HMM segmenters   (02_gaussian_hmm.py, …)
  cnv/              versioned CNV callers       (06_gene_cnv_caller.py,
                                                 07_plasmid_cnv_caller.py, …)
  evaluation/       versioned evaluators        (04_kpsc_evaluation.py, …)
  training/         dataset loader, trainer, inference (non-versioned)
  experiments/      one self-contained folder per experiment

diagnostics/    Streamlit app for interactive sample inspection
```

## Running an experiment

```bash
cd models/experiments/27
bash run.sh
```

Or directly (after `conda activate cnvrock`):
```bash
python models/train.py models/experiments/30/config.yaml
```

## Adding a new experiment

```bash
cp -r models/experiments/26 models/experiments/27
# edit models/experiments/27/config.yaml
```

A new experiment can reuse any existing versioned component — just point `architecture`, `hmm`, `cnv`, and `evaluation` in `config.yaml` at the same numbered files and adjust parameters. Only create a new versioned file (e.g. `08_plasmid_cnv_caller.py`) when the algorithm itself changes, not just the parameters.

Outputs are written to the `out_dir` defined in the config: `checkpoint.pth`, `latents.npy`, `reconstructions.npy`, `sample_ids.npy`, `segments.parquet`, `gene_calls.tsv`, `plasmid_gene_calls.tsv`, `evaluation.txt`.

## How it works

```mermaid
flowchart TD
    subgraph exp["models/experiments/N/"]
        cfg["config.yaml\narchitecture · hmm · cnv · evaluation\nhyperparameters · data paths"]
        runsh["run.sh"]
    end

    subgraph vc["Versioned components"]
        arch["architectures/\n06_conv_vae.py …"]
        hmm["hmm/\n02_gaussian_hmm.py …"]
        cnv["cnv/\n06_gene_cnv_caller.py\n07_plasmid_cnv_caller.py …"]
        ev["evaluation/\n04_kpsc_evaluation.py …"]
    end

    subgraph out["data/results/N/"]
        o1["latents.npy\nreconstructions.npy"]
        o2["segments.parquet"]
        o3["gene_calls.tsv\nplasmid_gene_calls.tsv"]
        o4["evaluation.txt"]
    end

    subgraph loop["Autonomous proposal loop"]
        claude["Claude Code\n/propose-experiment"]
        mail["📧 Proposal email"]
        you(["You"])
        daemon["Daemon\nlaunchd · 60 s"]
    end

    subgraph diag["Diagnostics"]
        app["Streamlit app\ndiagnostics/app.py"]
    end

    cfg -->|"selects component\nversions & params"| vc
    cfg --> runsh
    runsh -->|"wrap_up.py\ninference → HMM → CNV → eval"| out

    out -->|"evaluation.txt"| claude
    claude -->|"creates N+1 folder\n& config.yaml"| exp
    claude --> mail
    mail --> you
    you -->|"AUTHORISE"| daemon
    you -->|"feedback"| daemon
    daemon -->|"runs experiment"| runsh
    daemon -->|"on feedback: flags\nfor Claude to revise"| claude

    cfg -->|"resolves all paths\n& versions"| app
    out --> app
```

## Experiment proposal workflow

Claude analyses the latest `evaluation.txt`, proposes the next experiment, creates the folder, and emails a summary. Reply "AUTHORISE" to run it on the Mac mini; reply with feedback to get a revised proposal.

**First-time setup** (install the background polling daemon):
```bash
bash tools/install_daemon.sh
```

**To propose the next experiment** (invoke from Claude Code):
```
/propose-experiment
```
Claude sets up the experiment folder, writes a README, and emails a ≤100-line summary.

**The daemon** (`tools/check_and_run.sh`, running via launchd every 60 s) checks for a reply:
- `AUTHORISE` → runs the experiment automatically
- Anything else → flags that feedback is waiting; open Claude Code and run `/check-reply`

**Privacy:** `check_reply.py` searches only by the exact Message-ID of the proposal email. It never lists or reads any other email. Reply body content is never written to disk or logs.

## Diagnostics

```bash
cd diagnostics
streamlit run app.py
```

Select an experiment from the dropdown — the app loads that experiment's config to resolve data paths, component versions, and all calling parameters automatically. Displays the VAE latent space, per-sample read-depth tracks, HMM segmentation, and gene calls side by side.

## Data and reference

- **Cohort:** KpSC samples from the [AllTheBacteria](https://www.allthebacteria.org/) catalog (quality-filtered; SRA accessions in `assets/kpsc_bam_accessions.txt`)
- **Reference:** *K. pneumoniae* HS11286, NC_016845.1 (~5.3 Mb; GCF_000240185.1)
- **Extended reference:** HS11286 + plasmid contigs carrying resistance genes (`assets/HS11286_extended.fasta`)
- **Ground truth:** AMRFinder+ on AllTheBacteria assemblies + Kleborate v3 for porin/efflux status

## Plasmid gene expansion

New plasmid genes are added without any code changes:

1. Run `data/setup/add_phase_c_genes.py` — downloads representative plasmid from NCBI, BLASTs the gene, appends contig to `HS11286_extended.fasta`, adds a row to `assets/plasmid_refs/plasmid_gene_coords.tsv`
2. `sbatch hpc/build_extended_reference.sh` — re-index BWA
3. `sbatch --array=1-N%50 hpc/remap_unmapped_to_plasmids.sh` — remap per-sample
4. `python3 data/setup/merge_plasmid_counts.py ...` — update NPY store
5. Create next experiment config; all new genes are called automatically

## Installation

### Prerequisites

| Tool | Version | Purpose |
|------|---------|---------|
| conda / mamba | ≥23 | environment management |
| CUDA | 12.x | GPU training (CPU works but is ~10× slower) |
| BWA | ≥0.7.17 | read alignment |
| SAMtools | ≥1.20 | BAM processing |
| BLAST+ | ≥2.14 | plasmid gene discovery |

### 1. Clone and create the environment

```bash
git clone https://github.com/lcerdeira/CNVRock.git
cd CNVRock
conda env create -f environment.yml   # installs PyTorch + all dependencies
conda activate cnvrock
```

For CPU-only (no GPU), edit `environment.yml` before creating: replace `pytorch-cuda=12.4` with `cpuonly`.

### 2. Prepare read-count stores

The pipeline expects per-sample BAMs aligned to HS11286 (NC_016845.1). Starting from BAM files listed in `assets/kpsc_bam_accessions.txt`:

```bash
# Extract 1 kb bin read counts from each BAM
sbatch --array=1-N%50 hpc/collect_readcounts.sh

# Convert per-sample count files to a single NPY store
python data/setup/readcounts_to_npy.py \
    --counts-dir data/raw/ \
    --store-path data/inputs/KpSC-HS11286-1000bp-core-npy/
```

### 3. Prepare plasmid read-count store

Unmapped reads from step 2 are remapped to plasmid gene contigs:

```bash
sbatch hpc/build_extended_reference.sh          # BWA-index extended reference
sbatch --array=1-N%50 hpc/remap_unmapped_to_plasmids.sh
python data/setup/merge_plasmid_counts.py \
    --counts-dir data/inputs/plasmid_remap_counts/ \
    --store-path data/inputs/KpSC-plasmid-1000bp-npy/
```

### 4. Prepare ground truth

```bash
# AMRFinder+ copy counts (requires assemblies in data/assemblies/)
sbatch hpc/get_amrfinder_gt.sh

# Kleborate MLST / resistance (writes ST to assets/kpsc_kleborate_ground_truth.tsv)
sbatch hpc/get_kleborate_gt.sh
```

### 5. Run an experiment

```bash
# Full pipeline: train VAE → HMM segmentation → gene calls → evaluation
python models/train.py models/experiments/29/config.yaml

# Or submit to SLURM GPU queue
sbatch hpc/train_gpu.sh models/experiments/29/config.yaml
```

### 6. Launch the diagnostics app

```bash
cd diagnostics
streamlit run app.py
```

Select an experiment from the dropdown — the app auto-resolves all paths from `config.yaml`.

---

## Key findings

### blaCTX-M-15 false negatives: variant-specific, not lineage-specific

33 samples (6% of CTX-M-positive cases) are consistently missed across all experiments. Root cause: their CTX-M reads map to chromosomal blaSHV in the original BWA alignment and never appear as unmapped reads available for plasmid remapping.

ST11 has the highest FNR (0.46) because it is enriched for **blaCTX-M-65**, a variant not represented in the plasmid reference panel. Per-sample AMRFinder+ analysis confirms: 9/12 ST11 FNs carry CTX-M-65 exclusively (PCN ≈ 0.000), whereas ST11 CTX-M-15 carriers are detected normally (PCN 0.52–4.30). Adding a blaCTX-M-65 reference contig is expected to recover these cases.

### blaSHV amplification: assembly-based GT underestimates copy number

AMRFinder+ reports ≤ 1 copy for all 545 samples — consistent with assembly collapse of tandem duplications in short-read de novo assembly. CNVRock calls 34 samples as amplified (CRR 1.75–10.5×), with copy-ratio signal strongly above the chromosomal median. These are likely true tandem duplications invisible to assembly-based tools; long-read sequencing would be required to confirm.

---

## Setup

See [data/setup/](data/setup/) for all data-preparation scripts. The full workflow from BAMs to a trained model takes approximately 4–6 hours on an LSHTM HPC node (1× NVIDIA A100, 32 GB RAM, 8 CPUs).
