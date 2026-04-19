# Experiment 17 — Skip connections + CRR gate on long genes (04_gene_cnv_caller)

**Status:** Complete 2026-04-18

**Checkpoint and HMM segments reused from exp 16 — ultra-fast (CNV calling + evaluation only).**

## Hypothesis

Exp 16 revealed the root cause of the PPV collapse: `cnv_crr_amp_threshold` in
`03_gene_cnv_caller` only operates via the short-gene fallback path, which is gated on
`crr_min_bins_fallback=3`. GCH1 has 2 bins and is fallback-eligible; CRT (4 bins), MDR1
(7 bins), and PM2_PM3 (7 bins) are not. Raising the threshold to 1.40 was therefore a
no-op for the three problem genes — all their FPs come from the HMM calling CN≥2
unconditionally and the caller accepting every HMM amplification without a CRR check.

The fix is a new caller (`04_gene_cnv_caller`) that adds a CRR gate for long genes: after
the HMM calls CN≥2, reject the call (set to CN=1) if the gene/flank CRR is below
`cnv_crr_gate_threshold=1.40`. This directly targets the noise-level HMM calls (CRT FP
p50=1.12, PM2_PM3 FP p50=1.16) while preserving confirmed amplifications (CRT TP p10=1.57,
PM2_PM3 TP p25=1.42). Short genes (GCH1) keep the existing fallback logic unchanged, and
`cnv_crr_amp_threshold` is restored to 1.30 to recover the GCH1 FNR lost in exp 16.

## Changes from experiment 16

| Component                  | Exp 16                              | Exp 17                          | Rationale                                                        |
|----------------------------|-------------------------------------|---------------------------------|------------------------------------------------------------------|
| `cnv`                      | `03_gene_cnv_caller`                | `04_gene_cnv_caller`            | New caller adds CRR gate for long genes                          |
| `cnv_crr_gate_threshold`   | —                                   | 1.40                            | Veto HMM CN≥2 calls with CRR below this for long genes           |
| `cnv_crr_amp_threshold`    | 1.40                                | 1.30                            | Restored; only affects GCH1 fallback (1.40 hurt GCH1 FNR)       |
| `out_dir`                  | `16_skip_st0.80_crrthresh1.40`      | `17_skip_st0.80_crr_gate`       | Reflects the architectural change in calling                     |
| HMM / checkpoint           | segments from exp 15, copied to 16  | segments copied from exp 16     | Identical params; skip HMM step                                  |
| Everything else             | unchanged                           | —                               | Same architecture, self_transition=0.80                          |

## Expected outcome

- **CRT PPV**: should recover dramatically. Gate at 1.40 vs FP p90=1.29 means essentially
  all 375 FPs are cut. TP p10=1.57 means no TPs are lost. FNR stays 0.00.
- **MDR1 PPV**: should recover toward exp 13's 0.86. TP p10=1.69 means no TPs cut; most
  of the 182 FPs (FP p75=1.33) are below the gate.
- **GCH1**: should restore to exp 15 baseline (FNR≈0.18, PPV≈0.90). The fallback threshold
  is back to 1.30; GCH1 logic is unchanged from exp 15.
- **PM2_PM3**: mixed. Gate at 1.40 vs FP p90=1.31 cuts most of the 591 FPs. TP p25=1.42
  means ~15% of TPs (CRR 1.30–1.42) may be cut, pushing FNR from 0.24 toward ~0.28–0.32.
  PM2_PM3 PPV should improve from 0.06 to 0.30–0.50 even with some TP loss.

The key question remains PM2_PM3: if FNR stays below 0.33 and PPV reaches 0.30+,
the 05_conv_vae architecture is viable with the CRR gate in place.

## What we could do instead

1. **cnv_crr_gate_threshold=1.35** — if exp 17 shows PM2_PM3 FNR climbing above 0.33,
   a less aggressive gate preserves more borderline TPs (CRR 1.30–1.40) while still
   cutting most FPs (FP p90=1.31).
2. **self_transition=0.75 + gate=1.40** — if the gate suppresses too many TPs, a looser
   HMM (0.75) would generate more CN≥2 segments, recovering borderline TPs that the HMM
   with 0.80 self_transition is discarding (FN p50=1.27 confirms real signal remains).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.02 | 0.91 | 0.95 |
| GCH1    | 0.18 | 0.90 | 0.84 |
| MDR1    | 0.01 | 0.97 | 0.98 |
| PM2_PM3 | 0.41 | 0.69 | 0.64 |

**CRT and MDR1**: Exceeded expectations. CRT PPV recovered from 0.13 to 0.91 (only 5 FPs
remain, all with CRR ≥ 1.43 — just above the gate). MDR1 PPV hit 0.97 vs the predicted
0.86. FNR near-zero for both. Gate at 1.40 worked exactly as designed.

**GCH1**: Baseline precisely restored (FNR=0.18, PPV=0.90), confirming the fallback path
with crr_amp_threshold=1.30 is stable.

**PM2_PM3**: FNR=0.41 — worse than the predicted 0.28–0.32. The gate at 1.40 cut too many
TPs: FN CRR distribution is p50=1.28, p75=1.35, p90=1.38, meaning the bulk of the 20 FNs
have CRR below 1.40. These are genuine borderline amplifications that fall below the gate.
PPV=0.69 is better than predicted (0.30–0.50), showing the gate is highly precise but
over-aggressive for PM2_PM3.

**Next**: Experiment 18 — lower `cnv_crr_gate_threshold` from 1.40 to 1.35 to recover
PM2_PM3 FNs in the CRR 1.35–1.40 range while keeping CRT/MDR1 precision intact.
