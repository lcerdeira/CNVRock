# Experiment 20 — Lower GCH1 fallback threshold to 1.25

**Status:** Proposed 2026-04-19

**Reuses exp 19 VAE checkpoint + HMM segments + reconstructions — ultra-fast (CNV calling + evaluation only).**

## Hypothesis

Exp 19 fixed PM2_PM3 almost completely (FNR 0.74→0.01), but GCH1 remains at FNR=0.07
with 139 FNs. All FNs have CRR just below the fallback threshold:
p10=1.15, p25=1.21, p50=1.24, p75=1.28, p90=1.29.

GCH1 is a short gene (~2 bins at 1000bp resolution), so `fallback_eligible=True` in
`05_gene_cnv_caller`. It is called amplified purely via CRR threshold
(`cnv_crr_amp_threshold`), not via the fractional coverage / spanning logic. The current
threshold of 1.30 places virtually all FNs just below the call line.

Lowering `cnv_crr_amp_threshold` from 1.30 to 1.25 should recover FNs in the [1.25, 1.30]
range — roughly p50 to p90 of the FN distribution, i.e. ~50% of the 139 FNs. The risk of
additional FPs is low because pred_normal GCH1 p90=1.14, so very few true-normal samples
have CRR in [1.25, 1.30].

Note: `cnv_crr_amp_threshold` affects only GCH1 (the sole fallback-eligible gene).
PM2_PM3, MDR1, and CRT use the fractional coverage + gate path — they are unaffected.

## Changes from experiment 19

| Component               | Exp 19 | Exp 20 | Rationale                                              |
|-------------------------|--------|--------|--------------------------------------------------------|
| `cnv_crr_amp_threshold` | 1.30   | 1.25   | Recover GCH1 FNs with CRR in [1.25, 1.30]             |
| VAE checkpoint          | exp 19 symlink | exp 19 symlink | No retraining                          |
| HMM segments            | exp 19 symlink | exp 19 symlink | No re-segmentation                     |
| Reconstructions         | exp 19 symlink | exp 19 symlink | No re-inference                        |
| Everything else         | unchanged | — | Same fractional coverage caller, gate=1.35, HMM params |

## Expected outcome

- **GCH1 FNR**: recover FNs in [1.25, 1.30] — roughly the top 50% of the FN distribution
  (p50=1.24 to p90=1.29). Predicted FNR ≈ 0.03–0.04.
- **GCH1 PPV**: pred_normal p90=1.14 leaves a large gap below 1.25; minimal new FPs expected.
  PPV should stay ≥ 0.89.
- **PM2_PM3 / CRT / MDR1**: unchanged — `cnv_crr_amp_threshold` is not used in their
  calling path.

If GCH1 FNR barely improves, the FNs in [1.25, 1.30] are not a simple threshold problem —
they may reflect VAE signal weakness for AF-E historical samples (pre-2010), and the next
step would be investigating population-stratified training.

## What we could do instead

1. **Gate at 1.20** — if most FNs are in [1.20, 1.25], pushing lower recovers more but
   risks PPV degradation if some pred_normal samples reach that range.
2. **Population-specific thresholds** — AF-E (FNR=0.32) drives most of the GCH1 FNR.
   A gene×population-specific gate would concentrate the recovery where needed without
   loosening the global threshold.
3. **self_transition=0.75** — still untested. Longer HMM segments might improve signal
   aggregation for borderline GCH1 bins.
