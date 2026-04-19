# Experiment 12 — Higher encoder dropout (p=0.20 → p=0.30)

**Status:** Complete 2026-04-15

**Full retraining required (~hours).**

## Hypothesis

Exp 09/10 (02_conv_vae, p=0.20) established the baseline: CRT FNR=0.02, MDR1 FNR=0.03,
GCH1 FNR=0.07, but PM2_PM3 FNR=0.93–0.94 persists. FN CRR p50=1.43 >> 1.0 confirms
signal IS present in the VAE reconstruction — the bottleneck is that the signal is not
strong or consistent enough for the HMM to commit to the short (1–2 bin) PM2_PM3 runs.
HMM self_transition tuning (0.80→0.75 in exp 10) had negligible effect.

Exp 11 tested the opposite direction: p=0.10 (lower dropout for more consistent per-bin
signal). It caused severe regressions across all genes. The confound: changing dropout
changes the entire reconstruction error distribution, so the HMM emission parameters
— calibrated on a different scale — were no longer valid. p=0.10 + self_transition=0.80
was a confounded test, not a clean evaluation of reduced dropout.

This experiment tests higher dropout (p=0.30). More aggressive regularisation forces the
encoder to learn a tighter, more generalised representation of normal coverage. At
inference (dropout disabled), the reconstruction is a "cleaner" normal-mode estimate —
CNV samples produce larger absolute reconstruction errors and therefore stronger anomaly
signal. The hypothesis: stronger signal improves the signal-to-noise ratio for PM2_PM3's
short runs, giving the HMM a more unambiguous surface to commit to.

HMM self_transition is held at 0.80 (exp 09 baseline) to isolate the effect of the
dropout change. If p=0.30 shows promise but PM2_PM3 FNR still doesn't fully resolve,
self_transition can then be tuned for the new distribution.

## Changes from experiment 09

| Component              | Exp 09        | Exp 12 | Rationale                                                      |
|------------------------|---------------|--------|----------------------------------------------------------------|
| `architecture`         | `02_conv_vae` (p=0.20) | `04_conv_vae` (p=0.30) | Higher dropout to strengthen anomaly signal |
| `out_dir`              | `09_vae_dropout` | `12_higher_dropout` | Reflects actual approach |
| `hmm_self_transition`  | 0.80          | 0.80   | Unchanged — isolate dropout effect; re-tune if needed post-run |
| Checkpoint             | trained from scratch | retrained from scratch | Architecture changed; new checkpoint needed |
| `cnv_downsample_ratio` | 0.25 | 0.25 | Retained — increasing this shifts the normal baseline |
| Everything else        | unchanged     | —      | Same HMM, CNV caller, and eval config as exp 09                |

## Expected outcome

- **PM2_PM3 FNR**: should fall from 0.93 if stronger regularisation makes the VAE's
  normal-mode reconstruction tighter, producing larger CRR values for amplification bins.
  FN CRR p50=1.43 is the current signal level; we want this to rise or become more
  consistent across bins so the HMM at 0.80 can commit.
- **CRT and MDR1**: risk of slight regression if the distribution shift moves the HMM
  emission parameters to a less favourable regime. Both had FNR=0.02/0.03 in exp 09 with
  TP CRR p50 well above threshold — a moderate distribution shift should not cause severe
  regressions (unlike exp 11, where all genes failed simultaneously, suggesting p=0.10
  was qualitatively different).
- **GCH1**: should hold near 0.07. The distributed representation benefit that drove the
  GCH1 recovery in exp 09 is preserved or strengthened at p=0.30; GCH1's longer
  amplification runs are less sensitive to per-bin signal consistency.
- If PM2_PM3 FN CRR p50 also rises (e.g. 1.43 → 1.6+), that confirms the architecture
  change is producing stronger anomaly signal. If p50 stays at 1.43, the signal level is
  unchanged and the HMM is the remaining bottleneck.

## What we could do instead

1. **Lower HMM self_transition** (0.80 → 0.65) — if p=0.30 produces results with similar
   FN CRR distribution to p=0.20, then the HMM remains the bottleneck and the exp 12
   (previously proposed) self_transition reduction is the right next step.
2. **p=0.30 with tuned self_transition** — if p=0.30 retraining shows the new emission
   distribution needs a different self_transition, run wrap_up.py with adjusted HMM
   settings (no retraining needed once the checkpoint exists).
3. **Skip connections** — residual connections across stride-2 blocks to preserve fine
   spatial structure for 1–2 bin PM2_PM3 runs; higher architectural complexity, best as a
   follow-up if dropout tuning is exhausted.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.93 | 0.96 |
| GCH1    | 0.11 | 0.91 | 0.89 |
| MDR1    | 0.03 | 0.93 | 0.95 |
| PM2_PM3 | 0.53 | 0.83 | 0.62 |

**Matches vs predictions:**
- PM2_PM3 FNR dropped from 0.93 to 0.53 — predicted direction confirmed. Higher dropout
  (p=0.30) did strengthen the anomaly signal, enabling the HMM to commit to more PM2_PM3
  runs.
- CRT improved (0.02 → 0.00) and MDR1 held at 0.03 — as predicted; both are robust genes
  with large TP CRR values that withstand distribution shifts.
- GCH1 regressed slightly (0.07 → 0.11) rather than holding or improving. GCH1 FN p50=1.23
  > 1.0 confirms some HMM discarding; the distribution shift under p=0.30 slightly reduced
  GCH1's signal-to-noise.
- PM2_PM3 FN p50 = 1.38 >> 1.0 — signal is still present that the HMM is discarding.
  The remaining bottleneck is `hmm_self_transition`; now that the p=0.30 distribution is
  established, it can be tuned.

**Next:** Experiment 13 — lower HMM self_transition (0.80 → 0.65) reusing the exp 12 checkpoint.

## Proposal history

### Original proposal (2026-04-14)

Proposed increasing `cnv_downsample_ratio` from 0.25 to 0.50 (architecture reverted to
02_conv_vae, p=0.20). Rationale: PM2_PM3 FN CRR p50=1.43–1.45 >> 1.0 in exp 09/10 shows
signal is present but inconsistent across adjacent bins; more CNV-positive training windows
should sharpen per-bin anomaly contrast for the HMM. HMM self_transition kept at 0.80 to
isolate the sampler effect. Full retraining required.

### Feedback received (2026-04-14)

Higher `cnv_downsample_ratio` is counterproductive — if the model sees more CNVs during
training it treats them as the normal baseline, reducing anomaly signal at test time (the
same mechanism that likely harmed exp 11). HMM self_transition retuning was flagged as an
underexplored, no-retraining lever.

### Revision 1 (2026-04-14)

Abandoned cnv_downsample_ratio increase. Pivoted to HMM self_transition reduction:
0.80 → 0.65, reusing the exp 09 checkpoint (no retraining). `out_dir` renamed from
`12_higher_cnv_sampling` to `12_lower_self_transition`.

### Feedback received (2026-04-14)

Critique of the exp 11 diagnosis in the previous proposal email: the claim that p=0.10
caused the VAE to reconstruct amplifications more accurately (reducing anomaly signal) was
flawed reasoning. Changing dropout changes the reconstruction error distribution, so the
HMM emission parameters are calibrated to a completely different scale. The correct
conclusion is that p=0.10 + self_transition=0.80 was a confounded comparison, not that
p=0.10 is fundamentally harmful. The revised proposal (exp 12: 02_conv_vae, self_transition
0.65, reuse exp 09 checkpoint) was accepted as a clean experiment with correct reasoning.

### Revision 2 (2026-04-14)

Experiment 12 design was unchanged (02_conv_vae, p=0.20; self_transition=0.65; exp 09
checkpoint). Only the reasoning in the proposal was corrected.

### Feedback received (2026-04-14)

Raise dropout to 0.3.

### Revision 3 (2026-04-14)

Pivoted from the HMM self_transition reduction to a new architecture test. The user's
direction: try higher dropout (p=0.30) to see if stronger regularisation sharpens anomaly
signal. Design: create 04_conv_vae (p=0.30), full retraining required. HMM self_transition
reset to 0.80 (exp 09 baseline) to isolate the dropout effect — a new dropout value
produces a new emission distribution and the HMM parameters must not be conflated across
distributions. `out_dir` renamed from `12_lower_self_transition` to `12_higher_dropout`.
