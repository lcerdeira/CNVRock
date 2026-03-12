# experiments

One subdirectory per training run, named `exp_<NNN>_<short_description>`
(zero-padded to three digits).

## Contents of each experiment directory

```
exp_001_baseline/
├── README.md          ← hypothesis, result, conclusion (commit this)
├── config.yaml        ← frozen snapshot of the config used (do not edit)
└── results/
    ├── metrics.json   ← train/val loss and task metrics at each epoch
    └── plots/         ← loss curves, confusion matrices, etc.
```

Checkpoints are stored under `weights/exp_001_baseline/` (see `weights/README.md`).

## README template for each experiment

```markdown
# exp_NNN — <short description>

**Date:** YYYY-MM-DD
**Config:** configs/experiments/exp_NNN_....yaml
**Weights:** weights/exp_NNN_.../best.ckpt

## Hypothesis
What change is being tested and why it was expected to help.

## Result
Key metric(s) vs. baseline (e.g. "val AUC 0.91 vs 0.88 baseline").

## Conclusion
Whether the change should be promoted to the next baseline, and any
follow-up actions.
```

## Workflow

- Results directories (`results/`) are gitignored — only narratives and
  frozen configs are committed.
- Never edit a `config.yaml` inside an experiment folder after the run;
  create a new experiment instead.
