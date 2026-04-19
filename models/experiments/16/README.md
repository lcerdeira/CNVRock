# Experiment 16 — Skip connections + raise cnv_crr_amp_threshold (1.30 → 1.40)

**Status:** Complete 2026-04-18

**Checkpoint and HMM segments reused from exp 15 — ultra-fast (CNV calling + evaluation only).**

## Hypothesis

Exp 15 (05_conv_vae, self_transition=0.80) showed that raising the HMM prior alone is
insufficient to recover PPV. CRT PPV stayed at 0.13 and PM2_PM3 PPV stayed at 0.06 despite
transitions falling from 195 to 144. The 04_conv_vae architecture achieved PM2_PM3 PPV=0.73
with the same threshold (1.30), confirming the problem is not the threshold per se but the
05_conv_vae's broader anomaly distribution that places more noise calls above the threshold.

Evaluation CRR diagnostics show a clear separation between FPs and TPs across all genes:
- CRT: FP p50=1.12 vs TP p50=1.77
- MDR1: FP p50=1.25 vs TP p50=2.02
- GCH1: FP p50=1.23 vs TP p50=1.81
- PM2_PM3: FP p50=1.16 vs TP p50=1.51

Raising cnv_crr_amp_threshold from 1.30 to 1.40 should cut the bulk of the FP calls, which
have evaluation CRR clustered around 1.10-1.25, while preserving true positives whose CRR
is substantially higher. This is a pure CNV-calling-stage change — the HMM segments and
VAE checkpoint are identical to exp 15.

## Changes from experiment 15

| Component                | Exp 15                          | Exp 16                              | Rationale                                           |
|--------------------------|---------------------------------|-------------------------------------|-----------------------------------------------------|
| `cnv_crr_amp_threshold`  | 1.30                            | 1.40                                | Filter 05_conv_vae noise FPs that passed 1.30       |
| `out_dir`                | `15_skip_connections_st0.80`    | `16_skip_st0.80_crrthresh1.40`      | Reflects the threshold change                       |
| HMM                      | run from exp 14 inference       | segments.parquet copied from exp 15 | Identical params → identical segments; skip HMM step|
| Everything else          | unchanged                       | —                                   | Same architecture, self_transition=0.80, dropout    |

## Expected outcome

- **CRT PPV**: should recover dramatically. FP eval CRR p90=1.29 means most FPs will be
  filtered at the higher threshold; TP CRR p10=1.57 means no TPs are cut. FNR stays 0.00.
- **MDR1 PPV**: should recover toward exp 13's 0.86. TP p10=1.69 means no TPs are cut;
  the threshold will filter most of the 182 FPs.
- **GCH1 PPV**: should improve from 0.90 (already good). FP p75=1.33, TP p25=1.48 — most
  FPs cut, a small fraction of borderline TPs may be lost.
- **PM2_PM3**: mixed picture. FP p90=1.31 so most FPs will be cut, but TP p25=1.42 means
  some borderline TPs (CRR 1.30-1.40) may be cut → FNR could worsen from 0.24 to ~0.30.
  Even so, PM2_PM3 PPV should improve substantially from 0.06.
- **HMM diagnostics**: unchanged (callability=0.996, transitions p50=144) since HMM is not re-run.

The critical question is PM2_PM3: if FNR stays below 0.30 and PPV reaches 0.20+, the
skip-connection architecture is viable with threshold tuning. If FNR climbs above 0.35,
the threshold lever is insufficient and we need to adjust self_transition simultaneously.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.13 | 0.35 |
| GCH1    | 0.28 | 0.91 | 0.79 |
| MDR1    | 0.01 | 0.76 | 0.86 |
| PM2_PM3 | 0.24 | 0.06 | 0.20 |

HMM diagnostics: callability=0.996, transitions p50=144 (identical to exp 15 — segments were copied).

**Where predictions matched:** CRT FNR stayed at 0.00 and MDR1 FNR at 0.01. GCH1 PPV
improved marginally (0.90 → 0.91). PM2_PM3 FNR stayed at 0.24, not worsening past 0.30.

**Where predictions diverged:** CRT and PM2_PM3 PPV did not change at all (CRT=0.13,
PM2_PM3=0.06). The threshold increase had essentially no effect on the genes where FP
recovery was expected. GCH1 FNR worsened from 0.18 to 0.28 — the 1.40 threshold cut
borderline TPs (CRR 1.30–1.40 range) in the fallback path.

**Root cause:** `cnv_crr_amp_threshold` is only applied via the fallback path, which is
gated on `crr_min_bins_fallback=3`. Only GCH1 (2 bins) is fallback-eligible; CRT (4 bins),
MDR1 (7 bins), and PM2_PM3 (7 bins) rely entirely on the HMM call. The threshold change
was therefore a no-op for the three genes with collapsed PPV. All FPs for CRT, MDR1, and
PM2_PM3 originate from HMM segments calling CN≥2 over those genes — and the caller accepts
any HMM amplification call unconditionally for long genes.

**Next:** Exp 17 — new `04_gene_cnv_caller` adds a CRR gate for long genes: if HMM calls
CN≥2 but gene CRR < `cnv_crr_gate_threshold`, the call is rejected. Short genes keep the
existing fallback. Restore `cnv_crr_amp_threshold=1.30` for GCH1 to recover exp 15 FNR.

## What we could do instead

1. **self_transition=0.75 + cnv_crr_amp_threshold=1.40** — if exp 16 shows PM2_PM3 FNR
   climbing above 0.30, a slightly more permissive HMM (0.75) may recover borderline
   FNs whose CRR is near 1.27 (the current FN p50), while the 1.40 threshold still
   suppresses the noise FPs.
2. **cnv_crr_amp_threshold=1.35** — a less aggressive step if 1.40 proves too strict for
   PM2_PM3. PM2_PM3 TP p10=1.31 vs FP p90=1.31 suggests 1.35 cuts most FPs while
   keeping more borderline TPs.
