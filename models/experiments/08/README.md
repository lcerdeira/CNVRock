# Experiment 08 — Raise CNV confidence threshold (0.50 → 0.65)

**Status:** Proposed 2026-04-13

**No retraining. No HMM re-segmentation.** `cnv_min_confidence` is applied during
CNV calling (post-segmentation), so only the CNV caller and evaluation re-run.
Reuses exp 06 VAE checkpoint and exp 07 HMM segments.

## Hypothesis

Exp 07 raised `hmm_self_transition` 0.75 → 0.80, reducing FP counts modestly but
leaving 118 CRT FPs and 118 PM2_PM3 FPs, all clustering at CRR p50 ~1.20–1.21.
PM2_PM3 FNR rose from 0.20 to 0.26 in the process. The self_transition lever is
nearly exhausted: FP CRR p50=1.20 and FN CRR p50=1.33 are only 0.13 apart, so
any further tightening would cut TPs alongside FPs.

The orthogonal approach is `cnv_min_confidence`: this gates on the HMM's probability
mass in the amplified state at call time. Even if the HMM enters the amplified state
at CRR ~1.20, it likely does so with lower confidence than at CRR ~1.64–1.88. Raising
the threshold from 0.50 to 0.65 should suppress these uncertain borderline calls without
altering the HMM transition dynamics that affect FN boundaries.

## Changes from experiment 07

| Parameter             | Exp 07 | Exp 08 | Rationale                                            |
|-----------------------|--------|--------|------------------------------------------------------|
| `cnv_min_confidence`  | 0.50   | 0.65   | Filter borderline calls by HMM confidence, not penalty |

Everything else (VAE checkpoint, HMM segments, self_transition, downsampling) is unchanged.

## Expected outcome

- **CRT**: 118 FPs should reduce; FNR should remain 0.00 since TPs have CRR p50=1.88
  — they should be called with high confidence. CRT FPs at p50=1.21, p90=1.38 — much
  lower CRR, likely lower confidence.
- **PM2_PM3**: 118 FPs should reduce; FNR may edge up slightly from 0.26. TPs are at
  CRR p50=1.64, FNs at 1.33. If FNs are being called with confidence just above 0.50,
  they may fall below the new 0.65 gate — main risk.
- **MDR1**: 46 FPs at CRR p50=1.20 should reduce further; FNR should stay near 0.01.
- **GCH1**: uses CRR fallback (2-bin gene), not HMM confidence — unaffected.
- Population-specific note: CRT in AF-E (PPV=0.03) and PM2_PM3 in AF-W (PPV=0.04)
  are likely genuine amplifications the GATK ground truth doesn't capture; apparent
  PPV improvement from filtering those calls may be spurious.

## What we could do instead

1. **Combine both levers** — leave self_transition at 0.80 and raise confidence to 0.70
   for a more aggressive filter; riskier for PM2_PM3 FNR.
2. **Revert self_transition to 0.77** and rely on confidence filtering instead — if
   PM2_PM3 FNR at 0.26 is already too high, a softer transition with tighter confidence
   may recover recall while still suppressing FPs.
