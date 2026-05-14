# Experiment 02 — HMM stickiness + CNV threshold tuning

**No retraining.** Reuses the experiment 01 VAE checkpoint. Only HMM segmentation
and CNV calling parameters change. Run via `wrap_up.py` (inference → HMM → CNV → eval).

## Hypothesis

Experiment 01 showed FN p50 well above 1.0 for every gene:

| Gene    | FNR  | FN p50 CRR | Interpretation               |
|---------|------|------------|------------------------------|
| PM2_PM3 | 1.00 | 1.53       | Signal present, HMM rejecting |
| GCH1    | 0.81 | 1.35       | Signal present, HMM rejecting |
| MDR1    | 0.46 | 1.60       | Signal present, HMM rejecting |
| CRT     | 0.23 | 1.67       | Signal present, HMM rejecting |

Per the evaluation guidance: FN p50 >> 1.0 means the HMM is seeing the elevated
copy ratio but labelling it CN=1 anyway — the transition matrix is too sticky
(`self_transition=0.99` makes the model strongly prefer staying in its current state).

Additionally, PM2_PM3 `pred_normal` CRR p50 = 0.90 (slightly below 1.0) suggests
chromosome 14 has a mildly depressed baseline in many samples, causing the 80%
CN=1 sanity check to disqualify otherwise-callable samples.

## Changes from experiment 01

| Parameter | Exp 01 | Exp 02 | Reason |
|---|---|---|---|
| `hmm_self_transition` | 0.99 | 0.95 | Reduce stickiness; let HMM transition into amplified states more readily |
| `cnv_min_cn1_proportion` | 0.80 | 0.70 | Accommodate slightly depressed chr14 baseline |
| `cnv_min_confidence` | 0.70 | 0.60 | Proportionally relaxed alongside proportion threshold |

PPV was 1.00 in exp 01 with zero FP, so we have room to trade stickiness for recall.

## Expected outcome

- FNR should decrease, particularly for PM2_PM3 and GCH1
- PPV may decrease slightly (lower self_transition → more transitions → possible FP)
- If PPV drops badly, the next step is to tune `cnv_min_confidence` back up

## Status

Proposed — awaiting AUTHORISE.
