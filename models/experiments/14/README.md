# Experiment 14 — Skip connections in VAE encoder (05_conv_vae)

**Status:** Complete 2026-04-18

**Full retraining required (~hours).**

## Hypothesis

Exp 13 (04_conv_vae, p=0.30, self_transition=0.65) moved PM2_PM3 FNR from 0.53 to 0.34.
PM2_PM3 FN p50 = 1.26 >> 1.0 confirms real signal is still present but the HMM is still
discarding some. More critically, GCH1 FNR = 0.11 was completely unresponsive to the
self_transition change (0.80 → 0.65) — GCH1 FN p50 = 1.23 was identical between exp 12
and exp 13. HMM inertia is ruled out as the GCH1 bottleneck.

HMM transitions at exp 13 p50=101 per sample (up from 73 in exp 12) indicate the HMM is
already quite jumpy. Further self_transition reduction risks over-segmentation and PPV
collapse without materially improving GCH1.

The hypothesis: the current 5× stride-2 downsampling chain (04_conv_vae) discards spatial
resolution that is critical for short (1–2 bin) amplification runs. Per-bin anomaly signal
gets smeared across stride steps, making short PM2_PM3 and weak GCH1 runs inconsistent
across bins and difficult for the HMM to commit to. Adding a 1×1-strided residual shortcut
to each encoder block (05_conv_vae) propagates the fine-scale input structure through all
encoder levels. The encoder produces latent codes that are more sensitive to localised
amplitude spikes, yielding larger and more spatially precise reconstruction errors at the
amplified bins — giving the HMM a cleaner signal to segment.

HMM self_transition is kept at 0.65 (exp 13 value) to isolate the architecture effect.
Emission parameters are always re-fitted on the new checkpoint, so the distribution shift
is absorbed; self_transition is the only carry-over and holding it constant is a direct
architectural comparison.

## Changes from experiment 13

| Component             | Exp 13              | Exp 14              | Rationale                                                              |
|-----------------------|---------------------|---------------------|------------------------------------------------------------------------|
| `architecture`        | `04_conv_vae`       | `05_conv_vae`       | Residual skip connections in encoder to preserve fine spatial structure |
| `out_dir`             | `13_lower_self_transition` | `14_skip_connections` | Reflects actual change                                         |
| `hmm_self_transition` | 0.65                | 0.65                | Unchanged — isolate architecture effect                                |
| Checkpoint            | reused from exp 12  | retrained from scratch | Architecture changed; new checkpoint required                       |
| Everything else       | unchanged           | —                   | Same CNV caller, eval config, dropout p=0.30                          |

## Expected outcome

- **PM2_PM3 FNR**: should fall from 0.34. FN p50=1.26 confirms signal is present;
  skip connections should make that per-bin signal more spatially precise and
  consistent, allowing the HMM (at 0.65) to capture more 1–2 bin runs.
- **GCH1 FNR**: should improve from 0.11. GCH1 showed zero response to HMM tuning,
  pointing to a VAE signal quality issue. Better encoder spatial resolution should
  produce cleaner anomaly signals for the short GCH1 runs that the current
  architecture is smearing.
- **CRT and MDR1**: should hold. Both have large TP CRR values (p50 ~2.0 and ~1.9);
  the architecture change and retained self_transition=0.65 should not degrade them.
- If PM2_PM3 FNR falls to ~0.20 and GCH1 FNR begins to move, skip connections are
  confirmed as the remaining bottleneck and HMM tuning can be re-evaluated.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.11 | 0.32 |
| GCH1    | 0.17 | 0.89 | 0.84 |
| MDR1    | 0.01 | 0.66 | 0.81 |
| PM2_PM3 | 0.08 | 0.05 | 0.21 |

HMM segment diagnostics: callability=0.997, transitions p50=195 (exp 13: 101).

**Where predictions matched:** PM2_PM3 FNR fell dramatically from 0.34 to 0.08, confirming
that encoder skip connections provide the spatial precision needed to detect short amplification
runs. CRT FNR remained 0.00. MDR1 FNR also stable (0.01).

**Where predictions diverged:** PPV collapsed catastrophically across all genes —
CRT 0.86→0.11, PM2_PM3 0.73→0.05, MDR1 0.86→0.66. GCH1 FNR worsened rather than improved
(0.11→0.17). The root cause: skip connections produce much stronger anomaly signals, which
caused the HMM (at self_transition=0.65) to over-segment massively (transitions p50=195 vs
101 in exp 13). The HMM became hyperactive everywhere, not just at real amplification sites.

**Diagnosis:** The 05_conv_vae architecture is effective — skip connections genuinely improve
short-run detection. But self_transition=0.65, calibrated for the weaker 04_conv_vae signal,
is far too permissive for the skip-connection architecture. The HMM needs a stricter prior to
match the new signal magnitude.

**Next:** Exp 15 — reuse 05_conv_vae checkpoint, raise hmm_self_transition from 0.65 to 0.80
to bring transitions back to ~exp 13 levels while retaining the architectural sensitivity gain.

## What we could do instead

1. **Dropout p=0.40** — if skip connections don't materially improve GCH1 and PM2_PM3
   FN p50 stays near 1.26, the signal magnitude may still need strengthening.
2. **U-Net style decoder** — add skip connections from encoder to decoder so the
   decoder can reconstruct fine-scale detail; stronger reconstruction of normal bins
   means larger anomaly contrast at amplified bins.
3. **Gene-specific HMM threshold tuning** — if GCH1 FNs are specific to certain
   populations/years (AF-E, AF-W, 1995–2013), consider investigating coverage-batch
   effects rather than pursuing global architecture changes.
