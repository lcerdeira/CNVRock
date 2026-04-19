# Experiment 07 — Raise HMM self-transition (0.75 → 0.80)

**Status:** Complete 2026-04-13

**No retraining.** HMM self_transition changed, so re-segmentation is required.
Reuses exp 06 VAE checkpoint.

## Hypothesis

Exp 06 fixed recall with CNV-sample downsampling (PM2_PM3 FNR: 0.98→0.20, CRT: 0.04→0.00)
but introduced systematic false positives: 142 CRT FPs, 137 PM2_PM3 FPs, 333 GCH1 FPs, 60
MDR1 FPs. All cluster at CRR p50 ~1.21 — borderline samples the improved VAE now emits a
detectable amplification signal for.

The root cause: `hmm_self_transition=0.75` was calibrated for the weak exp 05 VAE, which
needed permissive transitions to catch anything. With the exp 06 VAE providing a much stronger
signal, the HMM can now afford to be more selective. Raising self_transition to 0.80 makes
state transitions slightly more expensive, requiring a stronger sustained signal. Clear TPs have
CRR p50 ~1.63–1.88; FPs cluster at ~1.21. The gap is wide enough that a modest increase should
preferentially cut FPs with minimal recall loss.

## Changes from experiment 06

| Parameter             | Exp 06 | Exp 07 | Rationale                                         |
|-----------------------|--------|--------|---------------------------------------------------|
| `hmm_self_transition` | 0.75   | 0.80   | Filter borderline FPs introduced by stronger VAE  |

Everything else (VAE checkpoint, downsampling ratio, CNV caller) is unchanged.

## Expected outcome

- **CRT**: FP count should drop sharply from 142; FNR should stay at 0.00 since CRT TPs
  have CRR p50=1.88, well above the borderline FP range (~1.21).
- **PM2_PM3**: FP count should drop from 137; FNR may rise slightly from 0.20 since FNs
  have CRR p50=1.30, but TPs are at p50=1.63 — a 0.33 gap gives room to filter FPs first.
- **MDR1**: small FP count (60) should largely clear; near-zero FNR should hold.
- **GCH1**: uses CRR fallback (2-bin gene), not HMM — self_transition has no effect here.
  Results should be unchanged from exp 06.
- Risk: if 0.80 overshoots and cuts genuine TPs, MCC on PM2_PM3 will drop. FN p50 CRR=1.30
  is the warning signal — if FNR climbs back significantly, step down to 0.77 in exp 08.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  | FN p50 CRR | FP count | Notes                                              |
|---------|------|------|------|------------|----------|----------------------------------------------------|
| CRT     | 0.00 | 0.31 | 0.55 | N/A (0 FN) | 118      | FPs reduced 142→118; PPV still low                 |
| GCH1    | 0.17 | 0.83 | 0.81 | 1.22       | 318      | Unchanged — uses CRR fallback, unaffected by HMM   |
| MDR1    | 0.01 | 0.92 | 0.95 | 1.38       | 46       | FPs reduced 60→46; near-perfect overall            |
| PM2_PM3 | 0.26 | 0.36 | 0.51 | 1.33       | 118      | FPs down 137→118, but FNR rose 0.20→0.26           |

**Prediction vs reality:**
- FP reduction was directionally correct but modest. CRT: 142→118, PM2_PM3: 137→118. CRT PPV
  improved 0.27→0.31 but remains low. The borderline FPs (CRR p50 ~1.20–1.21) survived; they
  sit just far enough below the TP range (CRT TPs p10=1.50, PM2_PM3 TPs p10=1.31) that the
  transition penalty at 0.80 does not cleanly exclude them.
- PM2_PM3 FNR rose 0.20→0.26, as flagged as a risk. FN p50 CRR=1.33 — right at the edge of
  the FP range (p90=1.38). Further tightening self_transition would push more FNs below the
  call threshold.
- GCH1 was correctly predicted to be unaffected (CRR fallback, no HMM involvement).
- MDR1 improved cleanly: both FP count and PPV moved in the right direction.

**Population patterns (from by-population tables):**
- CRT FPs are concentrated in AF-E (PPV=0.03 — essentially all positive calls are "false") and
  AS-SE-E (PPV=0.20). These may be genuine amplifications absent from the GATK ground truth.
- PM2_PM3 AF-W has PPV=0.04 (6,342 samples). This population's PM2/PM3 biology is likely
  under-captured by the Pf8 GATK pipeline — most calls may be real signal.

**Root cause of remaining FPs:** The self_transition lever is nearly exhausted. FP CRR p50=1.20
and FN CRR p50=1.33 are too close (0.13 gap) for a stricter transition to discriminate them
cleanly. The next mechanism is confidence-based post-processing: `cnv_min_confidence`.

→ See experiment 08 (baseline evaluation with 03_pf9_evaluation segment diagnostics)

## What we could do instead

1. **Step up further** (0.80 → 0.85) if PPV is still low — more conservative, lower FP
   at risk of recall regression on PM2_PM3.
2. **Raise cnv_min_confidence** (0.50 → 0.65) — gates out uncertain HMM calls rather than
   changing the transition penalty; independent of self_transition, could combine both.

## Proposal history

**Original proposal (exp 07 v1, 2026-04-13):** Raise `hmm_self_transition` 0.75 → 0.80 to
filter borderline FPs introduced by the stronger exp 06 VAE. Evaluation module unchanged
(01_pf9_evaluation).

**Feedback received:** HMM self_transition change approved. Two additions requested:
(1) don't over-weight PPV — ground truth is imperfect and some apparent FPs may be real
signal; update guidance text accordingly. (2) Add per-(Year, Population) stratified metrics
to the evaluation report, gated on minimum group size to avoid blowing up the report.

**Revised proposal (exp 07 v2, 2026-04-13):** HMM change unchanged. Evaluation upgraded to
`02_pf9_evaluation`: adds BY YEAR × POPULATION tables (one per gene, groups with n_eval ≥ 10
only), updates guidance text to note ground-truth imperfection. Config gains `eval_min_group_n: 10`.
