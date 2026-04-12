# Experiment 03 — Aggressive HMM + threshold tuning

**Status:** Proposed 2026-04-12

**No retraining.** Reuses the experiment 01 VAE checkpoint.

## Hypothesis

Experiment 02 showed that reducing `hmm_self_transition` from 0.99 → 0.95 improved recall
on CRT and MDR1, but all genes still have FN p50 >> 1.0, meaning the HMM continues to
discard real amplification signal. PM2_PM3 remains catastrophically broken (FNR=1.00, only
2 TPs) due to two compounding issues: (1) residual HMM stickiness, and (2) the CN=1
sanity check rejecting samples with a depressed chr14 baseline (pred_normal CRR p50=0.90).

Since PPV has been 1.00 across all experiments, we have full headroom to push stickiness
lower and relax thresholds further.

## Changes from experiment 02

| Parameter                  | Exp 02 | Exp 03 | Rationale |
|----------------------------|--------|--------|-----------|
| `hmm_self_transition`      | 0.95   | 0.90   | FN p50 still >> 1.0 for all genes → still too sticky |
| `cnv_min_cn1_proportion`   | 0.70   | 0.55   | PM2_PM3 pred_normal CRR p50=0.90; 0.70 still rejects most chr14 samples |
| `cnv_min_confidence`       | 0.60   | 0.50   | Relaxed proportionally with the CN=1 threshold |

## Expected outcome

- GCH1, MDR1 FNR should drop further (signal is clear; HMM just needs to transition more)
- PM2_PM3 FNR should recover substantially — less sticky HMM plus a sanity check that
  actually passes chr14 samples
- PPV may decline slightly but should remain high given four experiments of zero FPs

## If this doesn't work

- If FN p50 is still >> 1.0 after self_transition=0.90: consider changing the number of
  HMM states, or a sliding-window approach for PM2_PM3
- If PPV drops severely: tighten thresholds selectively per-gene rather than globally
