# scripts

Thin entry-point scripts. No business logic lives here — each script parses a
config file and delegates entirely to library code in `src/deep_gw_cnv/`.

## Scripts

| Script | Purpose |
|--------|---------|
| `train.py` | Full training loop with checkpointing |
| `evaluate.py` | Evaluation on a held-out split; writes `metrics.json` |
| `predict.py` | Inference on new data; writes predictions to disk |
| `preprocess.py` | `data/raw/` → `data/processed/` pipeline |

## Usage

```bash
# Train
python scripts/train.py configs/experiments/exp_001_baseline.yaml

# Evaluate the best checkpoint from that experiment
python scripts/evaluate.py configs/experiments/exp_001_baseline.yaml --split test

# Preprocess raw data (run once before first training)
python scripts/preprocess.py --raw data/raw/ --out data/processed/
```

## Design rules

- Scripts accept exactly **one positional argument**: the config file path.
- All other options are in the config. Avoid free-standing `--flags` for
  anything that affects results (those belong in the config).
- Scripts exit with code 0 on success and non-zero on failure — safe to use
  in job-scheduler pipelines.
