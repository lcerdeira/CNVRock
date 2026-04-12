# Experiment 05 — Per-gene CRR fallback thresholds

**Status:** Proposed 2026-04-12

**No retraining. No re-segmentation.** Reuses exp 01 checkpoint and exp 04 HMM segments.
Only the CNV caller changes (03_gene_cnv_caller.py).

## Hypothesis

Exp 04 achieved near-perfect recall (CRT FNR=0.00, MDR1 FNR=0.02, PM2_PM3 FNR=0.06)
but at the cost of precision: CRT PPV=0.23 (77% of calls are false positives). The
global CRR fallback threshold of 1.30 is too permissive — it catches genuine signal
but also accepts noisy CN=1 samples that happen to have slightly elevated CRR.

Per-gene CRR distributions from exp 04 show gene-specific gaps between FPs and TPs:

| Gene    | FP p75 | FP p90 | TP p10 | Gap  | Action         |
|---------|--------|--------|--------|------|----------------|
| CRT     | 1.58   | 2.35   | 1.64   | 0.06 | Raise to 1.60  |
| MDR1    | 1.41   | 1.48   | 1.51   | 0.03 | Raise to 1.50  |
| GCH1    | 1.37   | 1.40   | 1.37   | 0.00 | Keep at 1.30   |
| PM2_PM3 | 1.43   | 1.46   | 1.40   | neg  | Keep at 1.30   |

## Changes from experiment 04

| Parameter                          | Exp 04 | Exp 05 | Rationale                                |
|------------------------------------|--------|--------|------------------------------------------|
| `cnv` version                      | 02     | 03     | Per-gene threshold support               |
| `cnv_crr_amp_threshold_crt`        | 1.30   | 1.60   | FP p75=1.58, TP p10=1.64 — clean gap    |
| `cnv_crr_amp_threshold_mdr1`       | 1.30   | 1.50   | FP p90=1.48, TP p10=1.51 — clean gap    |
| `cnv_crr_amp_threshold_gch1`       | 1.30   | 1.30   | TP/FP overlap — can't raise safely       |
| `cnv_crr_amp_threshold_pm2_pm3`    | 1.30   | 1.30   | TP/FP overlap — can't raise safely       |

## Expected outcome

- CRT: PPV 0.23 → ~0.50+ (eliminating FPs with CRR 1.30–1.58); FNR stays 0.00
- MDR1: PPV 0.84 → ~0.95+ (most FPs have CRR < 1.50); FNR stays ~0.02
- GCH1: unchanged (~0.29 FNR, 0.95 PPV)
- PM2_PM3: unchanged (~0.06 FNR, 0.95 PPV)

## If this doesn't work

- If CRT PPV remains poor: the residual FPs at CRR 1.58–2.35 may be genuine
  amplifications the Pf8 GATK ground truth missed — consider reviewing them manually.
- If GCH1/PM2_PM3 FPs are still problematic: consider a population-stratified
  threshold or a read-depth-weighted CRR to reduce noise in low-coverage samples.
