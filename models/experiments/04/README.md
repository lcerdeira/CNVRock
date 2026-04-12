# Experiment 04 -- HMM state boundary fix + CRR fallback

**Status:** Proposed 2026-04-12

**No retraining.** Reuses the experiment 01 VAE checkpoint.

## Root cause diagnosis

Three experiments of reducing self_transition (0.99 -> 0.95 -> 0.90) failed
to move PM2_PM3 recall.  FN p50 CRR stayed at 1.51 in both exp 02 and exp 03.

The reason is mathematical: HMM state 2 was at mean=2.0, putting the
CN=1/CN=2 Viterbi boundary at CRR=1.5 -- exactly where the FN signal lives.
Emission ratio at CRR=1.51 is only 1.04x; the Viterbi self_transition penalty
is 90x.  1.04^7 * (0.01/0.90) = 0.015 << 1, so the decoder always stays in
state 1 regardless of signal length.

For GCH1 (2 bins, FN p50=1.33): even with perfect boundary placement,
Viterbi cannot call a 2-bin gene when self_transition > ~0.2.

## Changes from experiment 03

| Component             | Exp 03 | Exp 04 | Rationale |
|-----------------------|--------|--------|-----------|
| hmm version           | 01     | 02     | State 2 mean 2.0 -> 1.5 |
| cnv version           | 01     | 02     | CRR fallback added |
| hmm_self_transition   | 0.90   | 0.75   | Needed for 7-bin path to flip at CRR~1.5 |
| cnv_crr_amp_threshold | (none) | 1.30   | GCH1 2-bin fallback; TN p90=1.09-1.17 |

HMM means change: [0,1,2,3,4,5] -> [0,1,1.5,2.5,3.5,5]
  CN=1/CN=2 boundary: 1.5 -> 1.25
  Emission ratio at CRR=1.51: 1.04x -> 1.68x
  With self_transition=0.75: path ratio 1.68^7 * (0.05/0.75) = 1.49 > 1

## Expected outcome

- PM2_PM3: pred_amp jumps from 12 to hundreds; FNR should drop from 1.00
- MDR1: same mechanism (FN p50=1.57, 7 bins)
- GCH1: CRR fallback catches ~50-60% of FNs (FN p25=1.20, p50=1.33 vs 1.30)
- CRT: FNR=0.21 should hold (FN p50=1.70, well above boundary)
- PPV: should remain 1.00 (TN p90 is 1.09-1.17, well below 1.30)

## If this does not work

- If PM2_PM3 FNR remains ~1.00: VAE absorbs chr14 signal; retrain with
  chr14 gene region masked during loss computation.
- If PPV drops: apply CRR fallback only to GCH1/PM2_PM3, raise to 1.35.
