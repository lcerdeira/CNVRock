# configs

YAML configuration files. Every value that influences a result — model
architecture, optimiser settings, data splits, augmentation probabilities —
must live here, **not** hard-coded in scripts or library code.

## Layout

```
configs/
├── base.yaml                       ← shared defaults inherited by all experiments
└── experiments/
    └── exp_001_baseline.yaml       ← per-experiment overrides on top of base
```

## Config schema (key top-level fields)

| Field | Type | Description |
|-------|------|-------------|
| `model` | str | Key into `deep_gw_cnv.models.REGISTRY` |
| `model_cfg` | dict | Keyword arguments forwarded to the model constructor |
| `dataset` | str | Key into `deep_gw_cnv.datasets.REGISTRY` |
| `data_dir` | path | Root of `data/processed/` |
| `train` | dict | `epochs`, `batch_size`, `lr`, `scheduler` |
| `augmentation` | dict | Per-transform enable flags and parameters |
| `output_dir` | path | Where results and checkpoints are written |
| `seed` | int | Global RNG seed for reproducibility |

## Agent instructions

To run a new experiment:
1. Copy an existing file: `cp configs/experiments/exp_001_baseline.yaml configs/experiments/exp_002_mychange.yaml`
2. Edit the fields you want to vary.
3. Run: `python scripts/train.py configs/experiments/exp_002_mychange.yaml`

The training script resolves `base.yaml` defaults and writes a frozen copy of
the merged config into the experiment output directory — results are always
reproducible from that snapshot.
