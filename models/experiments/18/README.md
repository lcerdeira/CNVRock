# Experiment 18 — Sin loss fine-tuning to boost integer signal (05_conv_vae)

**Status:** Complete 2026-04-19

**Fine-tunes from exp 14 checkpoint — full pipeline (training + HMM + CNV + eval).**

## Hypothesis

Exp 17 PM2_PM3 FNR=0.41 — worse than expected. The 20 FNs have CRR p50=1.28, and
lowering the gate threshold (the original exp 18 proposal) would only recover ~25%
of those FNs. The deeper issue is that the VAE's reconstruction at genuine amplification
sites isn't sharp enough: blurry reconstructions understate the per-bin CRR, leaving
borderline amplifications just below whatever gate we set.

A sin loss term — `λ · mean(sin²(π · recon))` applied to the log2(count+1)
reconstruction — creates pressure for the decoder to commit to discrete amplitude
levels. When x is an integer in log2 space, sin²(π·x) = 0; between integers it
peaks at 0.5. This doesn't map directly to copy number integers, but it does sharpen
the reconstruction surface: the model is penalised for hedging at fractional values.
A crisper reconstruction produces a larger input-vs-reconstruction residual at
amplified bins, pushing CRR higher on true TPs and potentially allowing the existing
gate at 1.40 to work without also cutting genuine amplifications.

## Changes from experiment 17

| Component                | Exp 17 / orig. 18 proposal     | Exp 18 (revised)                        | Rationale                                                        |
|--------------------------|--------------------------------|-----------------------------------------|------------------------------------------------------------------|
| `pretrained_checkpoint`  | exp 14 (implicit)              | exp 14 checkpoint (explicit in config)  | Fine-tune from trained base                                      |
| `sin_loss_max_weight`    | —                              | 1.0                                     | New integer-regularisation term (≈0.14% of recon loss at peak)  |
| `sin_loss_warmup_epochs` | —                              | 30                                      | Ramp from epoch 0 to epoch 30 before full weight                 |
| `warmup_epochs`          | 50                             | 0                                       | KL already warmed in exp 14; apply full beta immediately         |
| `epochs`                 | 300                            | 100                                     | Fine-tune budget; early stopping at patience=20                  |
| `lr`                     | 1.0e-3                         | 1.0e-4                                  | Reduced LR for fine-tuning                                       |
| `cnv_crr_gate_threshold` | 1.40                           | 1.40 (unchanged)                        | Hold gate constant to isolate sin loss effect on CRR             |
| `out_dir`                | `18_skip_st0.80_crr_gate1.35`  | `18_skip_st0.80_sin_finetune`           | Reflects actual approach                                         |
| Everything else          | unchanged                      | —                                       | Same HMM params, caller, architecture                            |

## Expected outcome

- **PM2_PM3 FNR**: should improve from 0.41 if sin loss pushes amplified-bin CRR
  above 1.40. FN p50=1.28 means the signal is real but weak — a crisper reconstruction
  at those bins should raise CRR. Predicted FNR 0.25–0.35.
- **PM2_PM3 PPV**: should stay high (≥0.70) since FPs are already cut by the 1.40
  gate; the sin loss doesn't directly affect the gate threshold.
- **CRT / MDR1**: minimal change expected. Their remaining FPs (CRR ≥ 1.43) are
  safely above the gate; their TPs have very high CRR and should remain called.
- **GCH1**: likely stable. Fallback path is CRR-threshold based; sharpened
  reconstruction may modestly lower FNR from 0.18.

If PM2_PM3 FNR doesn't improve, the FNs are structurally under-represented in training
(not a reconstruction sharpness problem) and we'd need a different approach (e.g. a
gene-region-masked sin loss, or loosening self_transition).

## What we could do instead

1. **Gate 1.35 (original proposal)** — still a valid quick diagnostic. Recovers ~25%
   of PM2_PM3 FNs in the 1.35–1.40 CRR range without requiring retraining.
2. **Gene-region-masked sin loss** — apply `sin²(π·recon)` only to bins within CNV
   gene flanks. Forces sharpness exactly where the HMM needs signal, avoids any
   possible reconstruction degradation on genome-wide normal regions.
3. **self_transition=0.75 + gate=1.40** — if reconstruction signal is fine but the HMM
   is too conservative to call borderline CN≥2 segments.

## Proposal history

### Original proposal (2026-04-18)

Lower `cnv_crr_gate_threshold` from 1.40 to 1.35 (no retraining; CNV calling +
evaluation only from exp 17 segments). Expected to recover ~25% of PM2_PM3 FNs
in the CRR 1.35–1.40 range.

### Feedback received (2026-04-19)

Three points: (1) use symlinks instead of copying large files between experiment
directories; (2) when FN CRR is low, think about boosting signal upstream rather
than only tuning the gate; (3) consider adding a sin loss term to encourage
integer-valued reconstructions, with appropriate scheduling, or other masking-based
signal-boost approaches.

### What changed and why

The original proposal was a gate tweak: fast to run but only addresses ~25% of FNs
and doesn't change the underlying signal. The feedback directed attention to upstream
signal quality. The revised proposal fine-tunes the VAE with a sin loss that
incentivises integer-level reconstruction sharpness. The gate threshold is held at
1.40 (same as exp 17) to isolate the effect of the architectural change on CRR.

The symlink approach applies to future experiments reusing HMM segments or other
large artifacts without modification — noted in run.sh and saved to memory.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.02 | 1.00 | 0.99 |
| GCH1    | 0.07 | 0.91 | 0.91 |
| MDR1    | 0.03 | 1.00 | 0.98 |
| PM2_PM3 | 0.74 | 0.95 | 0.50 |

**GCH1**: Dramatically improved as predicted — FNR 0.18→0.07. Sin loss sharpened
the reconstruction at GCH1 amplified bins, pushing CRR above the fallback threshold.
PPV held at 0.91 (191 FPs, CRR p50=1.34). Better than the predicted "likely stable".

**CRT / MDR1**: Unchanged as predicted. PPV improved slightly (CRT 0.91→1.00, MDR1
0.97→1.00) because the sin loss apparently eliminated the 5 remaining CRT FPs and
the 1 remaining MDR1 FP from exp 17.

**PM2_PM3**: Catastrophically worse — FNR 0.41→0.74. The predicted FNR was 0.25–0.35;
actual is more than double. Two failure modes identified post-hoc:

1. *Gate veto*: HMM called CN≥2 but CRR < 1.40 → downgraded to CN=1. These FNs
   have CRR < 1.40 (bottom half of FN distribution: p10=1.18, p25=1.32, p50=1.40).
2. *Non-spanning segments*: The sin loss sharpened HMM segmentation, producing shorter
   CN≥2 segments that together cover most of the gene but no single segment fully
   satisfies the strict spanning condition (x0 ≤ g_start AND x1 ≥ g_end). These FNs
   have CRR ≥ 1.40 (top half: p50=1.40, p90=1.57) — real signal, but the caller
   never sees a spanning segment and leaves cn at −1.

The sin loss side effect on segmentation was not anticipated; this points to a
structural fragility in `04_gene_cnv_caller`'s spanning condition that only became
visible when the VAE started producing sharper, shorter segments.

**Next**: Experiment 19 — new `05_gene_cnv_caller` with fractional gene coverage
(≥50% of gene bins covered by CN≥2 segments) + lower gate to 1.35. Addresses both
failure modes without retraining.
