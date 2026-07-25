# 8. Scaling study

The scaling study is the manuscript's main empirical result. **Same
architecture, same hyperparameters, same evaluation** — only the training-set
size varies.

## Tiers

| Experiment | Tier | n samples | Manifest |
|---|---|---|---|
| (Phase 1 baseline) | 545 | 545 | (legacy) |
| exp 32 | 5K | 5,000 | `assets/kpsc_expansion_subset_5k.tsv` |
| exp 33 | 10K | 10,000 | `assets/kpsc_expansion_subset_10k.tsv` |
| exp 34 | 20K | 20,000 | `assets/kpsc_expansion_subset_20k.tsv` |
| exp 36 | 40K | 40,000 | `assets/kpsc_expansion_subset_40k.tsv` |

Each manifest is a **strict superset** of the previous one — see
{doc}`04_subset_selection`.

## Headline results

```{note}
Results table populates as experiments complete. Currently scheduled:
exp 32–34 will retrain at MQ=20 once the 20K download finishes
(~3 days from 2026-05-17). Exp 36 follows after.
```

### Chromosomal blaSHV (extra-copy)

| Tier | MCC | FNR | PPV | call_rate | n_eval |
|---|---|---|---|---|---|
| 545 (Phase 1) | — | — | — | — | — |
| 5K | — | — | — | — | — |
| 10K | — | — | — | — | — |
| 20K | — | — | — | — | — |
| 40K | — | — | — | — | — |

### Plasmid genes (presence)

For each tier we report per-gene MCC across {`blaKPC`, `blaCTX-M`, `blaNDM`,
`blaOXA-48`, `qnrB1`, `aac6-Ib-cr`}.

| Tier | blaKPC | blaCTX-M | blaNDM | blaOXA-48 | qnrB1 | aac6-Ib-cr |
|---|---|---|---|---|---|---|
| 5K | — | — | — | — | — | — |
| 10K | — | — | — | — | — | — |
| 20K | — | — | — | — | — | — |
| 40K | — | — | — | — | — | — |

## Reading the curve

We expect:

- **Monotonic improvement in MCC** for under-represented genes (rare STs,
  rare plasmid carriages) as the training set grows.
- **Plateau** at some n* between 10K and 40K — that's the per-gene
  saturation point.
- For genes already at MCC ≈ 0.9+ at 545 samples (e.g. blaKPC in the Phase 1
  cohort), expect **little to no headroom** — diminishing returns from
  scaling. Scaling matters most for the long tail.

## Early observation (MQ=40 baseline, since superseded)

The initial exp 32 / 33 runs were done at MQ=40 — too strict for the
multi-mapping plasmid reads in the extended reference (see
{doc}`09_methods`). With those broken-plasmid results we already saw a
**striking jump in chromosomal blaSHV `call_rate`**:

| Tier | blaSHV call_rate (MQ=40 baseline) |
|---|---|
| 5K | 0.44 |
| 10K | 0.91 |

The signal is there — the rerun at MQ=20 should produce comparable
chromosomal numbers plus working plasmid detection.

## Reproducibility note

Every experiment's `evaluation.txt` is committed to `data/results/{exp}/`.
The commit hash that produced each result is recorded in the experiment
log. See {doc}`10_reproducibility` for the exact recipe to regenerate.
