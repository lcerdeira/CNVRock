# Experiment 06 — Higher-beta VAE retrain (max_beta 1.0 → 4.0)

**Status:** Proposed 2026-04-12

**Full retrain required.** max_beta changed, so the checkpoint from exp 01 is not reusable.
HMM self_transition held at 0.75 (exp 05 value) to isolate the VAE change.

## Hypothesis

Exp 05 results suggest the VAE latent space does not cleanly separate amplified from normal
coverage. Two signals point to this:

1. **GCH1 FNR=0.29, FN p50 CRR=1.18 ≈ 1.0** — the copy-ratio signal is barely above 1.0,
   meaning the VAE latent representation for these samples is not clearly elevated. If the
   VAE has learned a latent space where low-amplitude amplifications overlap with high-normal
   variation, the HMM will never reliably call them.

2. **Evaluation note**: "High FNR may indicate the model learned on too many CNV-positive
   samples, causing the CNV read-count profile to be treated as the default (normal) state."
   PM2 amplification is extremely common in the Pf9 dataset (drug resistance). If the VAE
   encodes this as normal, it will compress amplified coverage into the same latent region as
   truly normal samples.

Root cause: `max_beta=1.0` imposes weak KL regularisation, allowing the latent space to
develop an unstructured, entangled representation. Raising beta to 4.0 forces stronger
regularisation — the model must pack meaningful variation (including copy-number state) into
fewer effective latent dimensions, producing a more structured space where amplified and normal
samples separate more cleanly. The HMM then receives clearer input signal.

## Changes from experiment 05

| Parameter        | Exp 05 | Exp 06 | Rationale                                                   |
|------------------|--------|--------|-------------------------------------------------------------|
| `max_beta`       | 1.0    | 4.0    | Stronger KL pressure → better-structured latent space      |
| `warmup_epochs`  | 50     | 100    | Longer warmup so reconstruction stabilises before KL grows |

Everything else (architecture, HMM, CNV caller, thresholds) is unchanged.

## Expected outcome

- **GCH1**: primary target — better latent separation should lift FNR from 0.29 toward 0.15.
  FN p50 CRR=1.18 is low, but if the latent space currently mixes amplified with normal, a
  structured beta-VAE should pull them apart and raise the HMM emission signal.
- **PM2_PM3**: uncertain — FN p50 CRR=1.47 >> 1.0 is clearly an HMM transition issue, not
  a latent signal issue. Better representation may help marginally but is not the bottleneck.
- **CRT, MDR1**: should remain stable; changes are in the VAE, not the calling logic.
- Risk: higher beta can reduce reconstruction quality, leading to noisier latent trajectories
  that hurt HMM segmentation. If PPV drops sharply, beta is too high.

## What we could do instead

1. **Lower hmm_self_transition** (0.75 → 0.65) — directly fixes PM2_PM3 HMM bottleneck
   (FN p50 CRR=1.47), no retraining. Original exp 06 proposal; deferred to isolate VAE first.
2. **Increase latent_dim** (10 → 20) — more representational capacity; complementary to
   higher beta but less targeted at the latent structure problem.

## Proposal history

**Original proposal:** Lower `hmm_self_transition` 0.75 → 0.65 — exp 05 showed PM2_PM3
FNR=0.98 with FN p50 CRR=1.47, pointing to Viterbi discarding real signal. No retraining,
~25 min run time.

**Feedback received:** Improve VAE training instead.

**Revised approach:** Retrain with `max_beta` 1.0 → 4.0. The HMM self-transition change
addresses a downstream bottleneck; the feedback asks us to look upstream at the VAE
representation quality first. A better latent space benefits all genes and all downstream
components, whereas HMM tuning is a narrower fix.
