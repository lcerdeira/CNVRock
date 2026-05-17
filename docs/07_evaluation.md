# 7. Evaluation

## Ground truth

Two sources, joined on **run accession** (`sample_id`):

1. **AMRFinder+** (`assets/amrfinder_gt_expansion.tsv`) — assembly-derived
   gene presence calls for the plasmid panel: `blaKPC`, `blaCTX-M`, `blaNDM`,
   `qnrB1`, `blaOXA-48`, `aac6-Ib-cr`. Also provides a baseline `blaSHV`
   presence call.
2. **Kleborate v3** (`assets/kpsc_expansion_kleborate_gt_runlevel.tsv`) —
   overrides AMRFinder+ for `blaSHV` with **per-sample BLAST copy counts**
   from the assembly, so the chromosomal-CNV metric is "extra copy ≥ 2"
   rather than "presence vs absence".

Both tables are bridged from BioSample-level to run-level identifiers using
`assets/kpsc_expansion_metadata.tsv` (`sample_id` column may be
comma-separated for multi-run BioSamples — we keep the first run).

## Metrics

`models/evaluation/04_kpsc_evaluation.py` computes per-gene:

| Metric | Definition |
|---|---|
| **MCC** | Matthews correlation coefficient. The primary metric. |
| **FNR** | False-negative rate: fraction of true events called as normal. **Primary optimisation target** — minimising missed AMR calls. |
| **PPV** | Precision. Apparent FPs may be real (assembly-derived ground truth is imperfect), so PPV is interpreted with caution. |
| **call_rate** | Fraction of GT-callable samples for which CNVRock produced a non-missing call. |
| **n_eval** | Number of samples in the evaluation. |

For each gene, the eval also reports:

| Diagnostic | What it tells you |
|---|---|
| `CRR/PCN by predicted label` | Whether the underlying signal (`p10/p25/p50/p75/p90` of CRR for chromosomal genes, PCN for plasmid) is correctly stratified by the predicted label. If `pred_normal` and `pred_event` overlap, the threshold is wrong. |
| `delta` (model-added missingness) | Samples where ground truth is callable but the model failed HMM sanity. High delta = HMM convergence problems. |
| Stratification by ST | Per-ST MCC, FNR, PPV. Reveals if any clone is systematically mis-called. |
| Segment-level callability | What % of samples produced ≥ 1 non-CN1 segment. Low callability flags model-collapse modes. |

## Two-sided evaluation

The eval runs **separately** for chromosomal and plasmid genes:

- **Chromosomal:** uses `gene_calls.tsv` and chromosomal-flank CRR.
  Currently scored: `blaSHV` (extra-copy mode).
- **Plasmid:** uses `plasmid_gene_calls.tsv` and absolute PCN.
  Currently scored: `blaKPC`, `blaCTX-M`, `blaNDM`, `qnrB1`, `blaOXA-48`,
  `aac6-Ib-cr`.

## Hold-out subset

If `cnv_eval_holdout_tsv` is set in the config, the evaluation restricts to
the listed run accessions. This is how we report **out-of-distribution
performance** for the manuscript (Phase 1 holdout never seen at any tier of
the scaling study).

## Output

The headline summary `evaluation.txt` follows a fixed format with three blocks:

```
OVERALL                  MCC / FNR / PPV / call_rate / n_eval per gene
MISSINGNESS              per-gene gt_miss / pred_miss / delta
CRR/PCN BY LABEL         signal-conditioned-on-predicted-label percentiles
```

Plus optional **per-ST** and **segment diagnostics** sections.

```{tip}
For the manuscript, the most important rows are the OVERALL MCC values for
each gene at each scaling tier (5K → 10K → 20K → 40K). These populate the
{doc}`08_scaling_study` table.
```
