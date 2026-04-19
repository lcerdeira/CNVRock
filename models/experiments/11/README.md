# Experiment 11 — Lower encoder dropout (p=0.20 → p=0.10)

**Status:** Complete 2026-04-14

**Full retraining required (~hours).**

## Hypothesis

Exp 09/10 with the dropout VAE (p=0.20) achieved excellent GCH1 recovery (FNR 0.17→0.07)
but left PM2_PM3 FNR catastrophically high (0.93–0.94). The FN CRR p50=1.43 >> 1.0
confirms the VAE IS producing elevated signal for PM2_PM3 amplifications — the problem
is that the HMM cannot commit to these runs.

HMM self_transition tuning in exp 10 (0.80→0.75) had negligible effect: transitions p50
only rose from 53 to 59, and PM2_PM3 FNR barely moved (0.94→0.93). This eliminates the
HMM threshold as the primary bottleneck. The more likely cause is signal *consistency*:
dropout at p=0.20 introduces stochastic per-bin noise that causes PM2_PM3 amplification
runs — which appear to be very short (1–2 bins) — to look intermittent rather than
sustained. The HMM never sees enough consecutive elevated bins to trigger a transition.

Lowering dropout to p=0.10 reduces this inter-bin noise. The distributed-representation
benefit that recovered GCH1 is preserved (we are not removing dropout entirely), but the
per-bin signal becomes more consistent, giving the HMM a cleaner surface to segment.
HMM self_transition is reset to 0.80 (exp 09 baseline) to isolate the effect of the
architecture change.

## Changes from experiment 10

| Component           | Exp 09/10 | Exp 11 | Rationale                                           |
|---------------------|-----------|--------|-----------------------------------------------------|
| `architecture`      | `02_conv_vae` (p=0.20) | `03_conv_vae` (p=0.10) | Lower dropout for consistent PM2_PM3 signal |
| `hmm_self_transition` | 0.75 (exp 10) | 0.80 | Reset to exp 09 baseline to isolate dropout effect |
| `out_dir`           | `10_lower_self_transition` | `11_lower_dropout` | Reflects actual approach |
| Checkpoint          | reuses exp 09 | retrained from scratch | Architecture changed; new checkpoint needed |
| Everything else     | unchanged | — | Same HMM, CNV caller, and eval config |

## Expected outcome

- **PM2_PM3 FNR**: should fall meaningfully from 0.93. Lower dropout → more consistent
  CRR signal across adjacent bins → HMM at 0.80 can commit to short PM2_PM3 runs. How
  much it falls depends on whether signal consistency was the actual bottleneck; if
  PM2_PM3 FN CRR p50 also drops toward 1.0 (rather than staying at 1.43), that would
  indicate the architecture change is working at the right level.
- **GCH1 FNR**: should hold near 0.07 or improve slightly. The distributed representation
  benefit of dropout is preserved at p=0.10, so the GCH1 gains from exp 09 should not
  regress. A small improvement is possible if GCH1 also benefits from more consistent signal.
- **CRT and MDR1**: should hold near 0.02 and 0.03 — both had large TP CRR values (>1.7)
  well above any noise floor; unaffected by this dropout change.
- **Transitions**: p50 may shift modestly from 59 (less noisy signal → slightly different
  segmentation rhythm); should remain well below exp 08's 123.
- Risk: if GCH1 FNR regresses significantly, p=0.10 is too low and the distributed
  representation property has been undermined. In that case, p=0.15 would be the next
  step.

## What we could do instead

1. **Skip connections in the encoder** — add residual connections between conv blocks to
   preserve spatial fine-structure across the stride-2 downsampling; higher architectural
   complexity, higher potential upside for very short CNV runs.
2. **Modify weighted sampler** — increase `cnv_downsample_ratio` (0.25 → 0.5 or 1.0) to
   expose the model to more CNV-positive training examples, potentially sharpening the
   model's ability to distinguish PM2_PM3 amplification from normal; risk of skewing the
   normal baseline reconstruction.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  | Notes                                                                      |
|---------|------|------|------|----------------------------------------------------------------------------|
| CRT     | 0.69 | 0.89 | 0.53 | Major regression from 0.02 — HMM discarding clear signal (FN CRR p50=1.98) |
| GCH1    | 0.14 | 0.91 | 0.87 | Regression from 0.07 — FN CRR p50=1.23, slight HMM failure                |
| MDR1    | 0.30 | 1.00 | 0.83 | Major regression from 0.03 — HMM discarding clear signal (FN CRR p50=1.97) |
| PM2_PM3 | 0.98 | 0.80 | 0.14 | Further regression from 0.93 — FN CRR p50=1.65, only 8 TPs vs 316 FNs    |

**Where predictions matched:**
- None — all genes regressed relative to exp 09/10 baseline.

**Where predictions diverged:**
- PM2_PM3 FNR did NOT improve; it worsened (0.93 → 0.98).
- CRT and MDR1 suffered severe regressions despite identical HMM settings (self_transition=0.80). In exp 09, these genes had FNR 0.02 and 0.03 respectively with the same HMM. The only change is the architecture.
- All FN CRR p50 >> 1.0 across every gene simultaneously — an HMM-discarding pattern now present everywhere.

**Root cause (corrected):**
The exp 11 failure cannot be cleanly attributed to dropout reducing anomaly signal. Changing dropout changes the entire reconstruction error distribution, which means the HMM's Gaussian emission parameters — fitted on the new distribution — are no longer on the same scale as exp 09/10. A self_transition=0.80 that worked with p=0.20 may simply be miscalibrated for p=0.10's different noise profile. Exp 11 proves that p=0.10 + self_transition=0.80 does not work; it does not prove that p=0.10 is fundamentally harmful. Proper evaluation of p=0.10 would require also re-tuning HMM parameters for that distribution.

→ See experiment 12 (04_conv_vae p=0.30; higher dropout to strengthen anomaly signal; full retraining; self_transition held at 0.80)

## Proposal history

### Original proposal (2026-04-14)

Proposed lowering HMM `self_transition` from 0.75 → 0.70, reusing the exp 09 checkpoint
(no retraining). Rationale: exp 10 (0.80→0.75) had negligible effect on PM2_PM3 FNR
(0.94→0.93), and PM2_PM3 FP=0 indicated headroom to push further. Expected PM2_PM3 FNR
to fall meaningfully at 0.70.

### Feedback received

Moved away from HMM tuning. Suggested exploring training or architecture changes instead:
lower dropout, skip connections, flash attention, or modified weighted sampler.

### Revision (2026-04-14)

Pivoted from HMM tuning to architecture change. The HMM path (0.80→0.75) was already
showing diminishing returns; the feedback confirmed it is time to address the upstream
signal. Selected lower dropout (p=0.20 → p=0.10) as the most targeted single change:
it directly addresses signal consistency without adding architectural complexity. Skip
connections and flash attention remain viable follow-ups if this does not move PM2_PM3 FNR.
`out_dir` renamed from `11_lower_self_transition_2` to `11_lower_dropout` to reflect the
actual experiment. `hmm_self_transition` reset to 0.80 to isolate the dropout effect.
