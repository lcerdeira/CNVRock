#!/usr/bin/env python3
"""
External-truth validation (Phase E reviewer response — "read depth truth").

A reviewer asked what ground truth CNVRock is validated against. There is
no per-bin copy-number ground truth in public short-read data, so we
validate on two complementary axes:

  1. DETECTION accuracy — against AMRFinder+ gene presence/absence calls
     derived independently from genome assemblies. For each plasmid gene
     family we compute the ROC AUC of the CNVRock family-aggregated PCN
     separating AMRFinder-positive from AMRFinder-negative samples. A high
     AUC shows the depth-based PCN axis recovers the assembly-based call.

  2. QUANTIFICATION accuracy — the spike-in recovery RMSE reported in the
     VAE ablation (Figure 2): known copy numbers injected into real data
     and recovered (RMSE 0.71 for the VAE baseline). That is the
     synthetic-ground-truth half and is not recomputed here.

Per-bin copy-number ground truth (long-read / ddPCR) is acknowledged as a
limitation and scheduled for the companion bench-validation study.

Run locally:  python3 analysis/external_truth_validation.py
Output:       data/results/external_validation/auc_summary.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
CALLS = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
GT = REPO / "assets/amrfinder_gt_expansion.tsv"
OUT = REPO / "data/results/external_validation"

# CNVRock PCN columns -> AMRFinder gene-family column. The PCN of a family
# is the max PCN across its allele-variant members.
FAMILY = {
    "blaKPC":        ["pcn_blaKPC-2"],
    "blaNDM":        ["pcn_blaNDM-1", "pcn_blaNDM-5"],
    "blaCTX-M":      ["pcn_blaCTX-M-14", "pcn_blaCTX-M-15",
                      "pcn_blaCTX-M-27", "pcn_blaCTX-M-65"],
    "blaOXA-48-like":["pcn_blaOXA-48", "pcn_blaOXA-181"],
    "qnrB":          ["pcn_qnrB1"],
    "aac6-Ib-cr":    ["pcn_aac6-Ib-cr"],
}


def roc_auc(score, label):
    """ROC AUC via the Mann-Whitney U identity (no sklearn dependency)."""
    pos = score[label == 1]
    neg = score[label == 0]
    if len(pos) < 5 or len(neg) < 5:
        return np.nan, len(pos), len(neg)
    u, _ = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg)), len(pos), len(neg)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    calls = pd.read_csv(CALLS, sep="\t")
    gt = pd.read_csv(GT, sep="\t")

    df = calls.merge(gt, on="sample_id", how="inner",
                     suffixes=("", "_gt"))
    print(f"merged {len(df):,} samples with both CNVRock calls and "
          f"AMRFinder ground truth")

    rows = []
    for fam, pcn_cols in FAMILY.items():
        if fam not in gt.columns:
            continue
        pcn_cols = [c for c in pcn_cols if c in df.columns]
        if not pcn_cols:
            continue
        score = df[pcn_cols].max(axis=1).values
        label = df[fam].astype(float).fillna(0).values
        label = (label >= 1).astype(int)
        ok = ~np.isnan(score)
        auc, n_pos, n_neg = roc_auc(score[ok], label[ok])
        # PCN medians by truth class
        pcn_pos = float(np.nanmedian(score[ok][label[ok] == 1])) \
            if (label[ok] == 1).any() else np.nan
        pcn_neg = float(np.nanmedian(score[ok][label[ok] == 0])) \
            if (label[ok] == 0).any() else np.nan
        rows.append({"gene_family": fam, "AMRFinder_pos": n_pos,
                     "AMRFinder_neg": n_neg,
                     "ROC_AUC": round(auc, 4) if not np.isnan(auc) else np.nan,
                     "PCN_median_pos": round(pcn_pos, 3),
                     "PCN_median_neg": round(pcn_neg, 3)})
        print(f"  {fam:16s} AUC={auc:.4f}  (+{n_pos}/-{n_neg})  "
              f"PCN +{pcn_pos:.2f}/-{pcn_neg:.2f}")

    out = pd.DataFrame(rows)
    out_path = OUT / "auc_summary.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nwrote {out_path}")
    valid = out["ROC_AUC"].dropna()
    if len(valid):
        print(f"mean ROC AUC across {len(valid)} families: {valid.mean():.4f}")


if __name__ == "__main__":
    main()
