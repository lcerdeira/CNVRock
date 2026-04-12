# Experiment 06 — CNV-sample downsampling (ablation-clean)

**Status:** Proposed 2026-04-12 (revised x3)

**Full retrain required.** cnv_downsample_ratio is a training-time change; the exp 01/05
checkpoint was trained without it and cannot be reused.
HMM self_transition held at 0.75 (exp 05 value) to isolate the VAE change.

## Hypothesis

Exp 05 results suggest the VAE latent space does not cleanly separate amplified from normal
coverage. Two signals point to this:

1. **PM2_PM3 FNR=0.98, FN p50 CRR=1.47** — the copy-ratio signal is elevated (Viterbi
   discards it), but the underlying cause may be that the VAE learned amplified coverage as
   the default state. If ~50% of Pf9 samples carry PM2 amplification, the VAE sees amplified
   profiles far more often than normal ones during training.

2. **GCH1 FNR=0.29, FN p50 CRR=1.18 ≈ 1.0** — the latent representation for GCH1-amplified
   samples is barely distinguishable from normal, consistent with a model that encodes common
   CNV profiles as "normal."

Root cause: training-distribution bias. The fix is `cnv_downsample_ratio: 0.25` — CNV-positive
samples are drawn at 1/4 the rate of normal samples during training via WeightedRandomSampler.
This directly rebalances the effective training distribution so amplification is the exception,
not the baseline. All other parameters remain at exp 05 values for a clean ablation.

## Changes from experiment 05

| Parameter               | Exp 05 | Exp 06 | Rationale                                                         |
|-------------------------|--------|--------|-------------------------------------------------------------------|
| `cnv_downsample_ratio`  | —      | 0.25   | Rebalances training distribution — CNV-positive samples 4× rarer |

Everything else (architecture, HMM, CNV caller, thresholds, max_beta, warmup_epochs) is
unchanged from exp 05.

## Expected outcome

- **PM2_PM3**: primary target — PM2 is the most common CNV in Pf9 (~50% of samples), making
  it the dominant source of distribution bias. Downsampling corrects this directly. Expect FNR
  to improve from 0.98, though some signal may still be lost at the HMM step (FN p50 CRR=1.47
  indicates the signal is present but Viterbi discards it). If FNR improves only partially, a
  follow-up experiment can tune hmm_self_transition.
- **GCH1**: should also improve as the latent space is less contaminated by amplified-is-normal
  representations. Expect FNR to improve from 0.29.
- **CRT, MDR1**: should remain stable; changes are upstream of the calling logic.
- Risk: WeightedRandomSampler means CNV-positive samples appear less often per epoch. Monitor
  PPV alongside FNR — if reconstruction quality degrades on CNV-positive profiles, PPV could drop.

## What we could do instead

1. **Lower hmm_self_transition** (0.75 → 0.65) — directly addresses PM2_PM3 Viterbi bottleneck
   (FN p50 CRR=1.47), no retraining. Deferred to isolate the distribution fix first.
2. **Downsampling + higher beta** (max_beta 1.0 → 4.0) — the previous proposal; adds stronger
   KL regularisation on top of distribution fix. Could try if downsampling alone is insufficient.

## Proposal history

**Original proposal (exp 06 v1):** Lower `hmm_self_transition` 0.75 → 0.65 — exp 05 showed
PM2_PM3 FNR=0.98 with FN p50 CRR=1.47, pointing to Viterbi discarding real signal. No
retraining, ~25 min run time.

**First revision feedback:** Improve VAE training instead.

**First revision (exp 06 v2):** Retrain with `max_beta` 1.0 → 4.0 and `warmup_epochs` 50 → 100.
Stronger KL regularisation to force a more structured latent space; longer warmup to prevent
latent collapse.

**Second revision feedback:** Downsample CNV-positive samples during training.

**Second revision (exp 06 v3):** Added `cnv_downsample_ratio: 0.25` on top of higher-beta
changes — CNV-positive samples sampled at 1/4 the rate of normal via WeightedRandomSampler, to
directly address training-distribution bias from PM2 prevalence in Pf9.

**Third revision feedback:** Downsampling only — revert max_beta and warmup_epochs.

**Third revision (exp 06 v4, current):** Reverted `max_beta` to 1.0 and `warmup_epochs` to 50
(matching exp 05). Only `cnv_downsample_ratio: 0.25` remains as the single change. Clean
ablation that isolates the distribution fix without confounding it with beta or warmup changes.
