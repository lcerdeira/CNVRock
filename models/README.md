# models/

## Big picture — how the VAE signal works

The VAE produces a reconstruction of what "normal" read-count coverage should look like
for each sample, based on the overall coverage profile. At inference, CRR = input /
reconstruction: CNVs stick out because the reconstruction doesn't model them.

Two failure modes to always keep in mind:
- **Over-training**: the VAE memorises CNV profiles as normal → CNVs stop sticking out
  → HMM misses amplifications (high FNR).
- **Under-training**: reconstructions are too different from input everywhere → CRR is
  noisy → HMM produces many borderline calls (high FP).

CRR is a rough signal (gene/flank ratio over a 100 kb window). Do not over-optimise
HMM parameters to CRR percentiles — that is optimising a noisy proxy. When FP and FN
CRR distributions overlap, the problem is upstream (VAE reconstruction quality), not
a calling-layer problem. Address reconstruction quality first; tune the calling layer
after the signal is clean.

Architecture improvement directions (in rough priority order):
1. **Dropout** in encoder — prevents per-sample CNV memorisation (simplest, tried in exp 09)
2. **Masked training** — exclude CNV gene loci from reconstruction loss during training
3. **U-Net skip connections** in decoder — finer reconstruction in non-CNV regions
4. **Self-attention block** — captures long-range genomic dependencies

## Entry point

```bash
.venv/bin/python train.py path/to/experiment/config.yaml
```

Loads the config, dynamically imports the versioned components named in it, trains the model, then runs inference → HMM segmentation → CNV calling → evaluation in sequence.

## Modules

| Path | What it does |
|---|---|
| `train.py` | Reads `config.yaml`, loads versioned components, runs full pipeline |
| `architectures/01_conv_vae.py` | `ConvVAE` — 1-D convolutional encoder/decoder |
| `hmm/01_gaussian_hmm.py` | `run_hmm_all_samples`, `fit_hmm_sample` — Gaussian HMM segmentation |
| `cnv/01_gene_cnv_caller.py` | `run_cnv_calls`, `call_all_genes` — gene-level CNV calls from segments |
| `evaluation/01_pf9_evaluation.py` | `run_evaluation` — metrics vs Pf8 GATK ground truth |
| `training/dataset.py` | `ReadCountDataset` — reads NPY store, log2-normalises |
| `training/trainer.py` | `train_vae` — training loop with KL annealing + early stopping |
| `training/wrap_up.py` | `run_inference` — encodes all samples, saves latents + reconstructions |
| `experiments/` | One folder per experiment, each with a `config.yaml` and `run.sh` |

## Versioning

Each component family has its own folder with numbered implementation files (`01_*.py`, `02_*.py`, …). The `config.yaml` for each experiment declares which version to use:

```yaml
architecture: "01_conv_vae"
hmm:          "01_gaussian_hmm"
cnv:          "01_gene_cnv_caller"
evaluation:   "01_pf9_evaluation"
```

`train.py` loads these via `importlib` at runtime, so a new experiment can mix and match any existing versions just by editing the config. The `__init__.py` in each folder re-exports the latest version as a convenience for interactive use.

**Rule:** reuse an existing versioned file and change only the config parameters when only the parameters differ. Create a new numbered file only when the algorithm or interface itself changes.
