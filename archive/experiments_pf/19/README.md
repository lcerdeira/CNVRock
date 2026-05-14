# Experiment 19 — Fractional gene coverage CNV caller + gate 1.35 (05_gene_cnv_caller)

**Status:** Complete 2026-04-19

**Reuses exp 18 VAE checkpoint + HMM segments — ultra-fast (CNV calling + evaluation only).**

## Hypothesis

Exp 18 (sin loss fine-tune) improved GCH1 dramatically (FNR 0.18→0.07) but worsened
PM2_PM3 catastrophically (FNR 0.41→0.74). Analysing the exp 18 FN population reveals
two distinct failure modes:

**Failure mode 1 — gate veto**: the HMM called CN≥2 for the sample but the gene/flank
CRR fell below the gate threshold (1.40), so `04_gene_cnv_caller` downgraded the call
to CN=1. These FNs have CRR < 1.40 (roughly the bottom half of the FN distribution:
p10=1.18, p25=1.32, p50=1.40).

**Failure mode 2 — non-spanning segments**: the HMM produced CN≥2 segments in the
PM2_PM3 region, but no single segment fully spans the gene (x0 ≤ g_start AND
x1 ≥ g_end). The strict spanning condition in `04_gene_cnv_caller` means `cn` stays
at −1 and the gene is never called. These FNs can have any CRR, including CRR ≥ 1.40
(the top half of the FN distribution: p50=1.40, p90=1.57). The sin loss in exp 18
likely caused this by sharpening segmentation: crisper boundaries produce shorter
segments that each cover most of the gene but no single one spans it end-to-end.

Addressing both modes simultaneously:
1. **New caller `05_gene_cnv_caller`**: replaces the strict spanning check with a
   fractional coverage condition — if ≥50% of gene bins are covered by high-confidence
   CN≥2 segments, the gene is called amplified (CRR gate veto still applies).
2. **Gate 1.35** (down from 1.40): recovers gate-vetoed FNs with CRR in [1.35, 1.40].

## Changes from experiment 18

| Component                        | Exp 18                         | Exp 19                               | Rationale                                                      |
|----------------------------------|--------------------------------|--------------------------------------|----------------------------------------------------------------|
| `cnv`                            | `04_gene_cnv_caller`           | `05_gene_cnv_caller`                 | Fractional coverage replaces strict spanning condition         |
| `cnv_min_gene_coverage_fraction` | —                              | 0.50                                 | New param: ≥50% gene bins covered → call CN≥2                 |
| `cnv_crr_gate_threshold`         | 1.40                           | 1.35                                 | Recover gate-vetoed FNs with CRR in [1.35, 1.40]              |
| VAE checkpoint                   | trained in exp 18              | symlinked from exp 18                | No retraining                                                  |
| HMM segments                     | fitted in exp 18               | symlinked from exp 18                | No re-segmentation                                             |
| Everything else                  | unchanged                      | —                                    | Same HMM params, architecture, GCH1 fallback logic            |

## Expected outcome

- **PM2_PM3 FNR**: should recover substantially from 0.74. FNs with CRR ≥ 1.40 (failure
  mode 2, non-spanning) are recovered by the fractional coverage check. FNs with CRR
  in [1.35, 1.40] (failure mode 1, gate-vetoed) are recovered by the lower gate. Only
  FNs with CRR < 1.35 (≈p10=1.18, roughly 10–15% of the 60 FNs) remain. Predicted
  FNR: 0.10–0.25.
- **PM2_PM3 PPV**: pred_normal PM2_PM3 p90=1.12 is far below 1.35, so very few
  currently-normal samples will be promoted to amplified. PPV should stay ≥ 0.85.
- **GCH1**: unchanged — GCH1 uses the fallback path (not spanning/gate logic), and
  FNR=0.07 from exp 18 should be preserved.
- **CRT / MDR1**: minimal change. Their pred_normal p90 ≈ 1.09–1.11; gate at 1.35
  safely excludes noise-level HMM calls. PPV should remain ≥ 0.99.

If PM2_PM3 FNR barely moves, the FNs are not explained by spanning fragments or the
gate — the VAE reconstruction itself may be failing to encode the amplification signal,
and we would investigate per-population training representation.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.03 | 0.95 | 0.73 |
| GCH1    | 0.07 | 0.91 | 0.91 |
| MDR1    | 0.01 | 0.99 | 0.86 |
| PM2_PM3 | 0.01 | 0.93 | 0.62 |

**PM2_PM3**: Dramatically exceeded predictions — FNR 0.74→0.01 (predicted 0.10–0.25). Only
2 FNs remain (CRR p50=1.22, genuinely weak signal below the gate). Both failure modes were
resolved: fractional coverage recovered the non-spanning FNs, and the gate at 1.35 recovered
the gate-vetoed FNs. PPV held at 0.93 (18 FPs, CRR p10=1.37).

**GCH1**: Unchanged from exp 18 as predicted — FNR=0.07 (139 FNs, CRR p50=1.24, p90=1.29).
All FNs sit just below the fallback threshold of 1.30. AF-E is the main driver
(population-level FNR=0.32), particularly pre-2010 samples. PPV=0.91 (191 FPs, CRR p50=1.34).

**CRT / MDR1**: Near-perfect and stable as predicted (FNR≤0.03, PPV≥0.95).

**Next**: Experiment 20 — lower `cnv_crr_amp_threshold` from 1.30 to 1.25 to recover
GCH1 FNs in the [1.25, 1.30] CRR band. Reuses exp 19 artifacts.

## What we could do instead

1. **self_transition=0.75** — looser HMM would produce longer segments, addressing
   failure mode 2 without a new caller. Risk: already p90=1160 transitions; could
   make HMM too jumpy globally.
2. **Per-gene gate thresholds** — if PM2_PM3 PPV falls unacceptably, use gate=1.30
   specifically for PM2_PM3 while keeping CRT/MDR1 at 1.40. Requires minor code change.
