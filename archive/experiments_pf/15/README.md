# Experiment 15 — Skip connections + HMM self_transition 0.80

**Status:** Complete 2026-04-18

**Checkpoint reuse (no retraining) — fast run via wrap-up only.**

## Hypothesis

Exp 14 (05_conv_vae, self_transition=0.65) produced two contradictory results: a huge
PM2_PM3 FNR win (0.34→0.08) alongside a complete PPV collapse (CRT 0.86→0.11, PM2_PM3
0.73→0.05). HMM transitions p50 jumped from 101 (exp 13) to 195 (exp 14), and callability
hit 0.997 — nearly every sample triggered the amplification state.

The cause: skip connections in 05_conv_vae roughly doubled the reconstruction anomaly signal
amplitude. self_transition=0.65, calibrated for the weaker 04_conv_vae signal, is far too
permissive for the new architecture. The HMM fires continuously on noise that was previously
below its activation threshold.

The hypothesis: raising self_transition from 0.65 to 0.80 (the exp 12 reference value) will
bring transitions back to the exp 13 level (~100), recovering PPV, while 05_conv_vae retains
the higher sensitivity that drove the PM2_PM3 FNR improvement. The 05_conv_vae inference
outputs are identical regardless of self_transition — no retraining needed; only HMM
segmentation, CNV calling, and evaluation are re-run.

## Changes from experiment 14

| Component             | Exp 14                    | Exp 15                         | Rationale                                                   |
|-----------------------|---------------------------|--------------------------------|-------------------------------------------------------------|
| `hmm_self_transition` | 0.65                      | 0.80                           | Dampen over-segmentation; match exp 12 reference level      |
| `out_dir`             | `14_skip_connections`     | `15_skip_connections_st0.80`   | Reflects the HMM change                                     |
| Architecture          | 05_conv_vae (retrained)   | 05_conv_vae (checkpoint reused)| No retraining — inference outputs copied from exp 14        |
| Everything else       | unchanged                 | —                              | Same CNV caller, eval config, dropout p=0.30                |

## Expected outcome

- **HMM transitions**: should drop from p50=195 toward ~100, matching exp 13. If transitions
  remain above ~150, a further increase to 0.85 may be needed.
- **PM2_PM3 FNR**: should remain well below exp 13's 0.34. Some increase from exp 14's 0.08
  is acceptable if PPV recovers — the architecture sensitivity gain should survive a stricter
  prior. FN p50=1.22 confirms real signal remains; 0.80 should still allow the HMM to commit.
- **PM2_PM3 PPV**: should recover substantially from 0.05. A tighter HMM prior eliminates
  borderline segments flagged by the stronger skip-connection signal.
- **CRT PPV**: should recover from 0.11 toward exp 13's 0.86. CRT FNR=0.00 in both exp 13
  and 14 — the gain here is purely PPV cleanup.
- **GCH1 FNR**: should improve from exp 13's 0.11. GCH1 FN p50=1.20 in exp 14 confirms
  signal is present; the architecture is more sensitive to short runs even at 0.80.
- **MDR1**: should hold and recover. FNR=0.01 is excellent; PPV should recover from 0.66.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.13 | 0.35 |
| GCH1    | 0.18 | 0.90 | 0.84 |
| MDR1    | 0.01 | 0.76 | 0.86 |
| PM2_PM3 | 0.24 | 0.06 | 0.20 |

HMM segment diagnostics: callability=0.996, transitions p50=144 (exp 14: 195, exp 13: 101).

**Where predictions matched:** PM2_PM3 FNR rose from 0.08 (exp 14) to 0.24, which stayed below
exp 13's 0.34 — the HMM becoming more selective was expected. MDR1 and CRT FNR held steady.
MDR1 PPV partially recovered from 0.66 to 0.76.

**Where predictions diverged:** CRT and PM2_PM3 PPV did not recover — CRT stayed at 0.13
(expected exp 13's 0.86) and PM2_PM3 stayed at 0.06. HMM transitions fell only from 195 to
144, not toward the hoped-for ~100. GCH1 FNR worsened slightly to 0.18 rather than improving.

**Diagnosis:** Raising self_transition from 0.65 to 0.80 was insufficient. The 05_conv_vae
architecture produces a broader distribution of anomaly signals that pass the cnv_crr_amp_threshold
= 1.30 as noise — the same threshold that held PM2_PM3 PPV at 0.73 for 04_conv_vae. The PPV
collapse is a CNV-caller-level problem, not purely an HMM-level one. The FP vs TP evaluation
CRR separation is clear (CRT FP p50=1.12 vs TP p50=1.77; PM2_PM3 FP p50=1.16 vs TP p50=1.51),
suggesting that raising cnv_crr_amp_threshold is the right lever.

**Next:** Exp 16 — raise cnv_crr_amp_threshold from 1.30 to 1.40; HMM and checkpoint unchanged.

## What we could do instead

1. **self_transition=0.75** — if 0.80 is too strict and PM2_PM3 FNR climbs back toward 0.34,
   try the intermediate 0.75. The signal from 05_conv_vae is roughly 2× that of 04_conv_vae,
   so the sweet spot may sit between the 0.65 that was too loose and 0.80.
2. **Raise cnv_crr_amp_threshold (1.30→1.40)** — filter borderline-CRR calls at the caller
   level rather than the HMM level; combines well with any self_transition outcome.
