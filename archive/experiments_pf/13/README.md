# Experiment 13 — Lower HMM self_transition (0.80 → 0.65)

**Status:** Complete 2026-04-15

**No retraining required — reuses exp 12 checkpoint (~minutes).**

## Hypothesis

Exp 12 (04_conv_vae, p=0.30) confirmed that higher dropout strengthens PM2_PM3 anomaly
signal: FNR dropped from 0.93 to 0.53. However, PM2_PM3 FN p50 = 1.38 >> 1.0 and GCH1
FN p50 = 1.23 > 1.0 both confirm real signal is still being discarded by the HMM.

With the p=0.30 emission distribution now established, the HMM `self_transition` parameter
can be tuned without confounding the architecture. `self_transition` controls how sticky
the HMM is in its current state — at 0.80, the model requires a long run of elevated bins
before committing to a CN>1 state. Reducing it to 0.65 lowers this inertia, allowing the
HMM to commit to shorter, weaker (but real) amplification runs.

This is a no-retraining test: one parameter, same checkpoint, fast to run.

## Changes from experiment 12

| Component              | Exp 12        | Exp 13 | Rationale                                               |
|------------------------|---------------|--------|---------------------------------------------------------|
| `architecture`         | `04_conv_vae` | same   | Unchanged — reusing checkpoint                         |
| `hmm_self_transition`  | 0.80          | 0.65   | Lower inertia so HMM commits to PM2_PM3/GCH1 FN signal |
| `out_dir`              | `12_higher_dropout` | `13_lower_self_transition` | Reflects actual change |
| Checkpoint             | trained from scratch | reused from exp 12 | Architecture unchanged |
| Everything else        | unchanged     | —      | Same CNV caller and eval config                         |

## Expected outcome

- **PM2_PM3 FNR**: should fall from 0.53. FN p50=1.38 >> 1.0 confirms the signal is
  there — the HMM just needs less persuasion to commit. A 0.65 self_transition should
  allow more of the 1–2 bin PM2_PM3 runs to be captured.
- **GCH1 FNR**: should improve from 0.11. FN p50=1.23 > 1.0 indicates HMM discarding.
  Lower self_transition helps GCH1 as well.
- **CRT and MDR1**: should hold. Both have large TP CRR values (p50 ~2.0 and ~1.9)
  well above threshold; even a more trigger-happy HMM should not introduce false positives.
  PPV may decrease slightly if very noisy bins get misclassified, but FNR should stay low.
- If PM2_PM3 FNR drops substantially (e.g. toward 0.20–0.30) and PPV holds above 0.75,
  this confirms the HMM was the remaining bottleneck.

## What we could do instead

1. **Intermediate self_transition** (0.80 → 0.72) — if 0.65 overshoots and CRT/MDR1
   PPV collapses, a smaller step might thread the needle.
2. **Skip connections** — residual connections across stride-2 blocks to preserve fine
   spatial structure for 1–2 bin PM2_PM3 runs; best if HMM tuning is exhausted.
3. **Dropout p=0.40** — if signal strength (not HMM inertia) is still the limiting
   factor after seeing 0.65 results.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.86 | 0.93 |
| GCH1    | 0.11 | 0.91 | 0.89 |
| MDR1    | 0.02 | 0.86 | 0.91 |
| PM2_PM3 | 0.34 | 0.73 | 0.69 |

**Matches vs predictions:**
- PM2_PM3 FNR dropped from 0.53 to 0.34 — predicted direction confirmed. Lower
  self_transition (0.65) allowed the HMM to commit to shorter PM2_PM3 runs.
  PM2_PM3 FN p50 = 1.26 (was 1.38): signal improved, some still discarded.
- MDR1 FNR improved to 0.02 and CRT held at 0.00 — as predicted.
- GCH1 FNR unchanged at 0.11 despite the self_transition reduction. GCH1 FN p50 = 1.23
  is identical to exp 12 — the HMM change had zero effect on GCH1. This rules out
  HMM inertia as the GCH1 bottleneck; the issue is upstream (VAE signal quality).
- HMM transitions p50 rose from 73 (exp 12) to 101 — the HMM is already quite jumpy.
  Further self_transition reduction risks over-segmentation and PPV collapse.

**Next:** Experiment 14 — skip connections in VAE encoder (05_conv_vae) to preserve
fine spatial resolution for 1–2 bin PM2_PM3 runs and improve GCH1 signal quality.
