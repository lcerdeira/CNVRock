# deep-gw-cnv

Variational autoencoder-based genome-wide copy number variation detection

## Quick start

```bash
# 1. Install the package in editable mode (creates the deep_gw_cnv importable package)
pip install -e ".[dev]"

# 2. Ensure data is in place (see data/README.md for acquisition instructions)
ls data/raw/

# 3. Preprocess raw data once before first training run
python scripts/preprocess.py --raw data/raw/ --out data/processed/

# 4. Launch a training experiment
python scripts/train.py configs/experiments/exp_001_baseline.yaml
```

## Repository layout

```
src/deep_gw_cnv/     installable Python package
  models/            architecture definitions (nn.Module subclasses)
  datasets/          Dataset / DataLoader factory code
  augmentation/      signal-space transforms used during training
  utils/             shared helpers (logging, metrics, reproducibility)

configs/             YAML experiment configs (base defaults + per-experiment overrides)
scripts/             thin entry-point scripts: train, evaluate, predict, preprocess
experiments/         per-run artefacts: frozen config snapshot, results, notes
notebooks/           exploratory analysis, demos, post-hoc plots
data/                raw and processed data files  [gitignored]
weights/             model checkpoints             [gitignored]
assets/              images and static files used in documentation
```

## Running a new experiment

1. Copy an existing config:
   ```bash
   cp configs/experiments/exp_001_baseline.yaml configs/experiments/exp_002_mychange.yaml
   ```
2. Edit the fields you want to vary (model, learning rate, augmentation, …).
3. Run training:
   ```bash
   python scripts/train.py configs/experiments/exp_002_mychange.yaml
   ```
4. Results land in `experiments/exp_002_mychange/`. Add a `README.md` there
   summarising the hypothesis and outcome.

## Conventions

- **Config over code** — all hyperparameters live in YAML. No magic numbers in
  scripts or library code.
- **Registered names** — models and datasets are referenced by string keys in
  configs; see `src/deep_gw_cnv/models/__init__.py` and `datasets/__init__.py`.
- **Frozen configs** — `scripts/train.py` writes a resolved config snapshot into
  each experiment output directory so any run is fully reproducible from it.
- **Experiment READMEs** — every `experiments/<name>/README.md` records the
  hypothesis, key result, and conclusion. Commit these even when results are
  negative.
