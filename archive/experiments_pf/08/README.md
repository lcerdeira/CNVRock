# Experiment 08 — Baseline evaluation with segment diagnostics (03_pf9_evaluation)

**Status:** Complete 2026-04-13

**Evaluation-only — no training, HMM, or CNV calling required (~5 min).**
Reuses exp 07's `gene_calls.tsv` and `segments.parquet`.

## Hypothesis

Before modifying the VAE architecture in future experiments, we need a baseline
that includes the new `03_pf9_evaluation` metrics (callability, mean_transitions).
Without this baseline, we cannot tell whether any change in those diagnostics
is due to the new architecture or is just an artefact of changing the evaluation
at the same time as changing the model. This experiment fixes that.

## Changes from experiment 07

| Component    | Exp 07               | Exp 08               | Rationale                                            |
|--------------|----------------------|----------------------|------------------------------------------------------|
| `evaluation` | `02_pf9_evaluation`  | `03_pf9_evaluation`  | Adds callability (CN1 baseline) + transition percentiles (p10–p90) |
| Everything else | unchanged         | —                    | Same gene_calls.tsv and segments.parquet from exp 07 |

No model changes. Gene call metrics (FNR, PPV, MCC) are identical to exp 07 by construction.

## Expected outcome

- **FNR / PPV / MCC**: identical to exp 07 (same gene_calls.tsv).
- **callability**: establishes the baseline fraction of samples where the HMM
  ever leaves CN1 (haploid baseline), using the exp 06/07 VAE checkpoint.
  Future experiments with VAE changes can be compared against this number
  directly.
- **transitions (p10–p90)**: percentile distribution of within-chromosome CN
  state changes per sample. Near-zero p50 indicates the HMM is too sticky;
  high p90 means too jumpy. Useful as a sanity check on future HMM parameter
  changes.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  | Notes                                               |
|---------|------|------|------|-----------------------------------------------------|
| CRT     | 0.00 | 0.31 | 0.55 | Identical to exp 07 by construction                 |
| GCH1    | 0.17 | 0.83 | 0.81 | Identical to exp 07 by construction                 |
| MDR1    | 0.01 | 0.92 | 0.95 | Identical to exp 07 by construction                 |
| PM2_PM3 | 0.26 | 0.36 | 0.51 | Identical to exp 07 by construction                 |

**Segment diagnostics (exp 08 baseline):**
- Callability (samples with any non-CN1 segment): **0.995**
- Within-chrom CN transitions per sample: p10=11, p25=23, p50=123, p75=600, p90=1221
- n_samples analysed: 50,455

**Prediction vs reality:**
- FNR/PPV/MCC are exactly identical to exp 07, as expected — same gene_calls.tsv, different evaluation only.
- Callability at 0.995 confirms the HMM is not too conservative overall; almost every sample has at least one non-CN1 segment.
- The transition distribution (p50=123, p90=1221) indicates the HMM is quite active/jumpy. This high transition count is consistent with the HMM oscillating between CN states rather than maintaining sustained segments, which is characteristic of weak VAE reconstruction signal at the boundaries of amplified vs. normal coverage.
- These numbers now serve as the baseline for comparing future architecture changes.

→ See experiment 09 (VAE encoder dropout training with 02_conv_vae)

## What we could do instead

Skip the baseline and proceed directly to VAE architecture changes. Risk: if the
new evaluation metrics in future experiments look surprising, we won't know
whether that's a model effect or an evaluation artefact.

## Proposal history

**Original proposal (2026-04-13):** Raise `cnv_min_confidence` 0.50 → 0.65 to filter
borderline HMM calls at CRR ~1.20. No retraining required.

**Feedback received (2026-04-13):** Asked whether flash attention or dropout could
boost signal further.

**First revision (2026-04-13):** Proposal unchanged — flash attention not applicable
(Conv1d architecture); dropout would require retraining and root cause is not
representation quality at that margin. Confidence gate stood.

**Feedback received (2026-04-13):** Root problem is VAE reconstruction quality.
Suggested VAE architecture improvements and lightweight HMM segment diagnostics.

**Second revision (2026-04-13):** Pivoted to VAE encoder dropout (p=0.20) with new
`02_conv_vae` and added `03_pf9_evaluation` with segment diagnostics.

**Feedback received (2026-04-13):** Do not modify the VAE yet — establish a baseline
evaluation first so we have a point of comparison before architecture changes.

**Third revision (2026-04-13, current):** Evaluation-only. Runs `03_pf9_evaluation`
on exp 07's existing outputs. VAE dropout deferred to exp 09 once baseline is in hand.

**Feedback received (2026-04-13):** Two corrections: (1) the organism is haploid
(*P. falciparum*), so CN1 is the baseline — calling it CN2 is wrong throughout.
(2) Rather than reporting mean transitions, report percentiles (p10–p90) to
capture the distribution shape.

**Fourth revision (2026-04-13, current):** Fixed `03_pf9_evaluation.py`:
- `_segment_diagnostics` now checks `cn != 1` (not `cn != 2`) for callability.
- Replaced `mean_transitions` with `transitions_percentiles` (p10–p90).
- Updated all guidance text and report output to match.
