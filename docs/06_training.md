# 6. Training

CNVRock training is **one entry point**: `models/train.py`. It loads a YAML
config, builds the dataset, trains the VAE, runs inference + HMM segmentation,
calls per-gene CNVs (chromosomal + plasmid), and writes evaluation outputs.

## Entry point

```bash
python models/train.py models/experiments/32/config.yaml
```

## SLURM wrapper

`hpc/train_gpu.sh` requests a GPU node, activates the conda env, and `cd`s
into `models/` before invoking `train.py`. Submit with:

```bash
sbatch hpc/train_gpu.sh experiments/32/config.yaml
```

Run-time on an A40: **~4 min for 5K samples, ~6 min for 10K samples**
(150 epochs × ~40 batches × 128 samples/batch).

## Config schema

Every experiment lives at `models/experiments/{N}/config.yaml`. The configs
for exp 32–36 share **identical architecture, HMM, CNV-caller and threshold
parameters** — only `store_path`, `plasmid_store_path`, and `out_dir` vary
across the scaling tiers.

```yaml
# Modules (resolved via importlib at runtime)
architecture: "06_conv_vae"
hmm:          "02_gaussian_hmm"
cnv:          "06_gene_cnv_caller"
evaluation:   "04_kpsc_evaluation"

# Data (per-tier varies)
store_path:           "../../../data/inputs/KpSC-expansion-5k-mq20-1000bp-npy"
plasmid_store_path:   "../../../data/inputs/KpSC-expansion-5k-mq20-plasmid-1000bp-npy"
out_dir:              "../../../data/results/32_kpsc_expansion_5k"

# Ground truth
kpsc_gt_path:           "../../../assets/amrfinder_gt_expansion.tsv"
kpsc_kleborate_gt_path: "../../../assets/kpsc_expansion_kleborate_gt_runlevel.tsv"
kpsc_meta_path:         "../../../assets/kpsc_expansion_metadata_runlevel.tsv"

# Plasmid genes
plasmid_gene_coords_path: "../../../assets/plasmid_refs/plasmid_gene_coords.tsv"
pcn_absent_threshold:     0.20
pcn_amp_threshold:        1.50

# VAE
latent_dim:    10
batch_size:    128
epochs:        150
lr:            1.0e-4
weight_decay:  1.0e-5
max_beta:      1.0
warmup_epochs: 20
patience:      20

# HMM
hmm_n_states:        6
hmm_self_transition: 0.80
hmm_low_cov_threshold: 10

# Chromosomal CNV caller
cnv_min_cn1_proportion:         0.55
cnv_min_confidence:             0.50
cnv_flank_padding:              100000
cnv_crr_amp_threshold:          1.75
cnv_crr_gate_threshold:         1.75
cnv_crr_min_bins_fallback:      3
cnv_min_gene_coverage_fraction: 0.50

eval_min_group_n: 10
```

## Architecture: 1D Conv-VAE

`models/architecture/06_conv_vae.py`. Encoder takes the per-sample 5,334-bin
vector through three 1D convolutions + a dense head, producing a 10-dim
latent. Decoder mirrors the encoder with transposed convs back to bin space.

Training minimises a weighted ELBO with **β-warmup** (β=0 → β=1 over the
first 20 epochs) plus a **CNV-pattern alignment auxiliary loss** at weight
1.0 (warmup 30 epochs). The auxiliary loss pulls the latent toward
biologically-meaningful structure (preventing the VAE from collapsing to a
global-depth-only representation).

## Segmenter: Gaussian HMM

`models/hmm/02_gaussian_hmm.py`. Per-sample inference reconstructions are
re-normalised and segmented with a **6-state Gaussian HMM** (state means
initialised at CN ∈ {0, 0.5, 1, 1.5, 2, 3}, self-transition 0.8). Low-coverage
bins (<10 reads) are masked and re-imputed from neighbours.

## Output artefacts

After `train.py` completes, `out_dir/` contains:

```
checkpoint.pth                model + optimiser state at best epoch
training_log.tsv              per-epoch loss curves
reconstructions.npy           (n_samples, n_bins) imputed depth
latents.npy                   (n_samples, 10)
segments.parquet              per-sample CN segments
gene_calls.tsv                per-sample chromosomal-gene CN calls
plasmid_gene_calls.tsv        per-sample plasmid-gene CN calls
evaluation.txt                MCC/FNR/PPV by gene and by ST
```

`evaluation.txt` is the headline artefact for the manuscript — see
{doc}`07_evaluation`.
