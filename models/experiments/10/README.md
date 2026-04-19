# Experiment 10 — Lower HMM self_transition (0.80 → 0.75)

**Status:** Complete 2026-04-14

**No retraining — reuses exp 09 checkpoint (~minutes).**

## Hypothesis

Exp 09 (dropout VAE) achieved a major GCH1 improvement (FNR 0.17 → 0.07) but caused
a catastrophic PM2_PM3 regression (FNR 0.26 → 0.94). The segment diagnostics pinpoint
the cause: PM2_PM3 FN CRR p50=1.45 >> 1.0, which by the guidance means signal IS present
and the HMM is discarding it — the fix is not upstream reconstruction but the HMM's
willingness to commit.

With the dropout VAE in place, the reconstruction landscape changed: transitions are
much less jumpy (p50=53 vs 123 in exp 08), which is healthy. But the more conservative
HMM (self_transition=0.80) is now too sticky for PM2_PM3's amplification runs.
Lowering self_transition from 0.80 to 0.75 makes the HMM transition more readily,
targeting the PM2_PM3 FNs at CRR~1.45 without touching the VAE or CNV caller.

The 0.80 setting was needed in exp 07 specifically to suppress FPs introduced by the
non-dropout VAE at CRR~1.21. With the dropout VAE's cleaner baseline, those borderline
FPs at ~1.21 are already gone (PM2_PM3 FP=0, GCH1 FP=281 at CRR p50=1.35). Reverting
to 0.75 should let the HMM catch the PM2_PM3 FNs (CRR~1.45) without re-introducing
the ~1.21 FPs that motivated the 0.80 setting.

## Changes from experiment 09

| Component           | Exp 09 | Exp 10 | Rationale                                              |
|---------------------|--------|--------|--------------------------------------------------------|
| `hmm_self_transition` | 0.80 | 0.75   | HMM too sticky after dropout VAE; PM2_PM3 FNR=0.94    |
| `out_dir`           | `09_vae_dropout` | `10_lower_self_transition` | New output directory |
| Checkpoint          | trained from scratch | reuses exp 09 `checkpoint.pth` | HMM-only change |
| Everything else     | unchanged | —   | Same dropout VAE, CNV caller, and eval as exp 09       |

## Expected outcome

- **PM2_PM3 FNR**: should fall substantially from 0.94 — FN CRR p50=1.45 is well above
  the TP CRR p50=1.80, and 0.75 gives the HMM enough freedom to commit to those runs.
- **GCH1 FNR**: should hold near 0.07. GCH1 uses the CRR fallback path (< 3 bins), so
  it is not affected by self_transition. Any change in GCH1 FNR would be unexpected.
- **CRT and MDR1**: should hold. Both had FNR near 0.00 with TP CRR well above 1.5;
  the HMM change won't affect them materially.
- **GCH1 FPs**: may increase slightly (from 281) since 0.75 is more permissive. Watch
  whether FP CRR p50 stays near 1.35 — if they stay there it means those are borderline
  real signals, not noise artefacts.
- **Transitions**: p50 expected to rise slightly from 53 (more willingness to transition)
  but should remain well below the exp 08 oscillation level of 123.
- Risk: if PM2_PM3 or other gene FPs spike significantly, 0.75 is too permissive for
  the dropout VAE's signal level.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  | Notes                                                              |
|---------|------|------|------|---------------------------------------------------------------------|
| CRT     | 0.02 | 1.00 | 0.99 | Unchanged from exp 09; excellent                                    |
| GCH1    | 0.07 | 0.87 | 0.89 | Unchanged from exp 09; GCH1 unaffected by self_transition (fallback path) |
| MDR1    | 0.03 | 1.00 | 0.98 | Unchanged from exp 09; excellent                                    |
| PM2_PM3 | 0.93 | 1.00 | 0.27 | Barely changed from 0.94 — self_transition 0.80→0.75 almost no effect |

**Segment diagnostics (exp 10):**
- Callability: **0.989** (up from 0.986 in exp 09 — very slight loosening confirmed)
- Within-chrom CN transitions per sample: p10=6, p25=11, p50=59, p75=473, p90=1102
- Transitions p50 rose slightly from 53 → 59, confirming the HMM is marginally more willing to commit

**Where predictions matched:**
- GCH1 FNR held at 0.07 — uses CRR fallback path; unaffected by self_transition as predicted.
- CRT and MDR1 unchanged as predicted.
- PM2_PM3 FPs remain 0 — confirming no FP bleed-over from the lower threshold.

**Where predictions diverged:**
- PM2_PM3 FNR barely changed (0.94 → 0.93). FN CRR p50=1.43 >> 1.0 still confirms HMM
  is discarding clear signal. The 0.75 setting was expected to produce a substantial
  recovery; instead it had negligible effect. This suggests the PM2_PM3 amplification
  runs may be very short (1–2 bins) and the HMM needs a significantly more aggressive
  self_transition reduction to capture them.

→ See experiment 11 (lower HMM self_transition further: 0.75 → 0.70, reuse exp 09 checkpoint)

## What we could do instead

1. **Lower self_transition further** (0.75 → 0.70) — more aggressive unlock; warranted
   given 0.75 produced negligible improvement. PM2_PM3 FP=0 so there is headroom.
2. **Retrain with lower dropout** (p=0.20 → p=0.10) — might make PM2_PM3 amplification
   runs appear longer/more sustained; requires full retraining.
   Keep this as fallback if 0.70 also fails.
