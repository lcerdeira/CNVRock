# Experiment 05 — Bin-count-gated CRR fallback

**Status:** Proposed 2026-04-12

**No retraining. No re-segmentation.** Reuses exp 01 checkpoint and exp 04 HMM segments.
Only the CNV caller changes (03_gene_cnv_caller.py).

## Hypothesis

Exp 04's CRR fallback (override CN=1 → CN=2 when CRR ≥ 1.30) rescued GCH1 recall
by handling its 2-bin span, but it introduced ~190 CRT and ~128 MDR1 false positives
by overriding correct HMM CN=1 calls on genes that have enough bins for Viterbi to
work properly. The HMM boundary fix in exp 04 already gives CRT FNR=0.00 and
MDR1 FNR=0.02 without any CRR help — the fallback is redundant for these genes and
only adds noise.

The principled fix: **apply the CRR fallback only to genes with fewer than 3 bins**,
where the HMM is physically unable to transition regardless of self_transition.
This applies uniformly across the genome without requiring per-gene threshold tuning
against labelled ground truth.

Gene bin counts at 1 kbp resolution:
  GCH1    ~2 bins  → fallback applies   (2 < 3)
  CRT     ~4 bins  → HMM only           (4 ≥ 3)
  PM2_PM3 ~7 bins  → HMM only           (7 ≥ 3)
  MDR1    ~7 bins  → HMM only           (7 ≥ 3)

## Changes from experiment 04

| Parameter                  | Exp 04  | Exp 05  | Rationale                                     |
|----------------------------|---------|---------|-----------------------------------------------|
| `cnv` version              | 02      | 03      | Bin-count gate added                          |
| `cnv_crr_min_bins_fallback`| (none)  | 3       | Fallback only for genes < 3 bins              |
| `cnv_crr_amp_threshold`    | 1.30    | 1.30    | Unchanged — not tuned to labelled data        |

## Expected outcome

- CRT: FPs eliminated (HMM was already calling correctly; fallback no longer overrides)
- MDR1: FPs eliminated (same reason); FNR stays ~0.02
- GCH1: unchanged (2-bin fallback still active)
- PM2_PM3: slight FP reduction; small FNR risk if HMM misses some without the fallback

## Long-term: genome-wide CNV calling

The four labelled genes are a proxy for evaluation only. The design of this
caller — HMM for well-covered regions, CRR fallback for sub-resolution loci —
is intended to generalise. For unlabelled regions of the genome:

- The HMM will call CNVs in any region with ≥ 3 bins without needing thresholds
  tuned to labels.
- The Pf8 GATK ground truth is itself imperfect; apparent FPs may be real
  amplifications the pipeline missed. Threshold tuning against imperfect labels
  risks miscalibrating the caller for the rest of the genome.
- Future calibration approaches: population-level CRR distributions as a null,
  holdout regions with known copy number, or cross-validation across populations.

## Proposal history

**Original proposal (before feedback):** Per-gene CRR thresholds — CRT raised to 1.60,
MDR1 raised to 1.50, GCH1 and PM2_PM3 kept at 1.30.

**Feedback received:** (1) Pf8 GATK ground truth isn't perfect; some apparent FPs may
be real amplifications it missed. (2) The goal is genome-wide CNV detection, not 4-gene
tuning — per-gene thresholds overfit to labelled data and won't generalise.

**Revised approach:** Instead of per-gene thresholds, gate the CRR fallback on bin count.
The threshold (1.30) is unchanged. Only the *applicability* of the fallback changes —
it now requires fewer than 3 bins, which is a principled physical criterion that applies
uniformly across the whole genome.
