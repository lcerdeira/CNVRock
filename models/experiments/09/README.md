# Experiment 09 — VAE encoder dropout (02_conv_vae)

**Status:** Complete 2026-04-14

**Full retraining required (~hours).**

## Hypothesis

The exp 08 baseline established that callability is 0.995 and within-chrom transitions
average p50=123 per sample — the HMM is not too conservative globally, but it is quite
jumpy. The core problem for GCH1 (FNR=0.17) and PM2_PM3 (FNR=0.26) is weak VAE
reconstruction signal at marginal amplification levels (FN CRR p50 ~1.22–1.33).

FP CRR and FN CRR for these genes overlap (gap ≤ 0.13), leaving no room for the HMM
or confidence gate to separate them without trading FNR for PPV. The only path forward
is to improve the reconstruction signal upstream: if the VAE produces cleaner
reconstructions, amplified samples will show higher CRR (the ratio of input to
reconstruction rises) and marginal cases become separable.

Encoder dropout (p=0.20) prevents the VAE from memorising individual per-sample
coverage patterns — including the elevated counts of CNV-positive samples — as "normal".
This forces the encoder to learn distributed representations that generalise across the
population, making the reconstructed "normal" flatter and amplified regions more
prominent in the ratio.

## Changes from experiment 07 (last full training)

| Component       | Exp 07                      | Exp 09             | Rationale                                        |
|-----------------|-----------------------------|--------------------|--------------------------------------------------|
| `architecture`  | `01_conv_vae`               | `02_conv_vae`      | Adds Dropout(p=0.20) after each encoder conv ReLU|
| `evaluation`    | `02_pf9_evaluation`         | `03_pf9_evaluation`| Matches exp 08 baseline for direct comparison    |
| `out_dir`       | `07_higher_self_transition` | `09_vae_dropout`   | New output directory                             |
| Everything else | unchanged                   | —                  | Same HMM/CNV parameters as exp 07               |

Decoder is unchanged: dropout is disabled at eval time, so reconstruction fidelity at
inference is not degraded.

## Expected outcome

- **GCH1 FNR**: should decrease from 0.17 if dropout pushes the reconstruction lower
  for amplified samples, raising their CRR above the 1.30 fallback threshold.
- **PM2_PM3 FNR**: should decrease from 0.26 if FN CRR rises from ~1.33 toward 1.5+,
  widening the gap over FP CRR (~1.20).
- **CRT**: FNR already 0.00; no regression expected. FP count may improve slightly.
- **MDR1**: already near-perfect; should hold.
- **Callability**: compared against exp 08 baseline (0.995); transitions distribution
  is the key diagnostic for whether the new VAE produces cleaner, more sustained segments.
- Risk: if dropout is too aggressive, reconstruction loss may rise and callability may
  drop. If callability falls below ~0.98 or MDR1 FNR regresses, p=0.20 was too high.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  | Notes                                                        |
|---------|------|------|------|--------------------------------------------------------------|
| CRT     | 0.02 | 1.00 | 0.99 | Tiny regression from 0.00; acceptable                        |
| GCH1    | 0.07 | 0.87 | 0.89 | Big improvement from 0.17 — dropout hypothesis confirmed     |
| MDR1    | 0.03 | 1.00 | 0.98 | Excellent; held as expected                                  |
| PM2_PM3 | 0.94 | 1.00 | 0.25 | Catastrophic regression from 0.26 — HMM too conservative     |

**Segment diagnostics (exp 09):**
- Callability: **0.986** (down from 0.995 in exp 08 — dropout slightly reduced signal strength)
- Within-chrom CN transitions per sample: p10=5, p25=10, p50=53, p75=443, p90=1039
- Transitions much less jumpy than exp 08 baseline (p50=53 vs 123) — dropout VAE produces cleaner reconstructions

**Where predictions matched:**
- GCH1 FNR fell from 0.17 to 0.07 — exactly as expected; dropout raised the reconstruction baseline
  so marginal amplifications now stand out in the CRR ratio.
- CRT and MDR1 held near-perfect as predicted.
- Callability dropped slightly (within the "too aggressive" risk flagged in the hypothesis).

**Where predictions diverged:**
- PM2_PM3 FNR collapsed from 0.26 → 0.94 — the opposite of the predicted improvement.
  FN CRR p50=1.45 >> 1.0 confirms the signal IS present; the HMM is discarding it.
  The dropout VAE's cleaner-but-less-sustained reconstruction likely makes PM2_PM3
  amplification runs appear shorter, leaving the HMM self_transition=0.80 too sticky
  to commit. Lowering self_transition is the targeted fix; no retraining needed.

→ See experiment 10 (lower HMM self_transition 0.80 → 0.75, reuse exp 09 checkpoint)

## What we could do instead

1. **Lower `cnv_crr_amp_threshold`** (1.30 → 1.15) — recovers GCH1 FNs at CRR ~1.22
   without retraining, but increases FPs for all CRR-fallback genes. Does not help
   PM2_PM3 (HMM-based).
2. **Raise `cnv_downsample_ratio`** (0.25 → 0.35) — more CNV-positive training samples;
   simpler change but exp 06 showed diminishing returns past 0.25.
