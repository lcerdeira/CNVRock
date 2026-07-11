#!/usr/bin/env python3
"""
Threshold-sensitivity / ROC analysis for CNVRock detection calls (10K KpSC).

Reviewer concern: the per-gene absent/amplified thresholds (plasmid PCN
absent 0.20; chromosomal CRR amplified 1.75) could look ad hoc. This script
shows the calls are robust to threshold choice by reporting, per gene family:

  - ROC AUC of the continuous PCN separating AMRFinder+ present vs absent
    (threshold-independent separability), and
  - the MCC-vs-threshold curve, with the operating threshold marked, to show
    it sits on the plateau of near-optimal MCC rather than a sharp peak.

    python3 analysis/threshold_sensitivity.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
CALLS  = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
GT     = REPO / "assets/amrfinder_gt_expansion.tsv"
OUT    = REPO / "data/results/threshold_sensitivity"

# family PCN columns -> GT presence column
FAMILIES = {
    "blaKPC":         (["pcn_blaKPC-2"], "blaKPC"),
    "blaCTX-M":       (["pcn_blaCTX-M-14", "pcn_blaCTX-M-15",
                        "pcn_blaCTX-M-27", "pcn_blaCTX-M-65"], "blaCTX-M"),
    "blaNDM":         (["pcn_blaNDM-1", "pcn_blaNDM-5"], "blaNDM"),
    "blaOXA-48-like": (["pcn_blaOXA-48", "pcn_blaOXA-181"], "blaOXA-48-like"),
    "qnrB":           (["pcn_qnrB1"], "qnrB"),
    "aac6-Ib-cr":     (["pcn_aac6-Ib-cr"], "aac6-Ib-cr"),
}
PCN_ABSENT_THRESH = 0.20   # the operating point


def roc_auc(score, label):
    """AUC via the Mann-Whitney U identity (ties = 0.5)."""
    pos = score[label == 1]; neg = score[label == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    # rank-based AUC
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(len(order)); ranks[order] = np.arange(1, len(order) + 1)
    # average ties
    allv = np.concatenate([pos, neg])
    df = pd.DataFrame({"v": allv, "r": ranks})
    df["r"] = df.groupby("v")["r"].transform("mean")
    r_pos = df["r"].values[:len(pos)]
    auc = (r_pos.sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    return auc


def mcc_at(score, label, thr):
    pred = (score >= thr).astype(int)
    tp = int(((pred == 1) & (label == 1)).sum())
    tn = int(((pred == 0) & (label == 0)).sum())
    fp = int(((pred == 1) & (label == 0)).sum())
    fn = int(((pred == 0) & (label == 1)).sum())
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den > 0 else 0.0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    calls = pd.read_csv(CALLS, sep="\t")
    gt = pd.read_csv(GT, sep="\t")

    thresholds = np.round(np.arange(0.05, 1.55, 0.05), 2)
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    summary = []
    for fam, (cols, gtcol) in FAMILIES.items():
        if gtcol not in gt.columns:
            continue
        d = calls[["sample_id"] + cols].merge(
            gt[["sample_id", gtcol]], on="sample_id",
            how="inner").dropna(subset=[gtcol])
        score = d[cols].fillna(0).sum(axis=1).values.astype(float)
        label = (d[gtcol].astype(float) >= 1).astype(int).values
        auc = roc_auc(score, label)
        mccs = [mcc_at(score, label, t) for t in thresholds]
        best_t = thresholds[int(np.argmax(mccs))]
        mcc_op = mcc_at(score, label, PCN_ABSENT_THRESH)
        summary.append(dict(family=fam, n=len(d), n_pos=int(label.sum()),
                            roc_auc=round(auc, 3),
                            mcc_at_operating=round(mcc_op, 3),
                            best_threshold=best_t,
                            best_mcc=round(max(mccs), 3)))
        ax.plot(thresholds, mccs, marker="o", ms=3, lw=1.4,
                label=f"{fam} (AUC {auc:.2f})")

    ax.axvline(PCN_ABSENT_THRESH, ls="--", color="#b00", lw=1.2,
               label=f"operating PCN threshold {PCN_ABSENT_THRESH}")
    ax.set_xlabel("PCN presence threshold", fontsize=10)
    ax.set_ylabel("Matthews correlation coefficient (MCC)", fontsize=10)
    ax.set_title("Detection robustness: MCC vs PCN threshold, per gene family\n"
                 "(10K tier vs AMRFinder+ presence)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT / "threshold_sensitivity.png", dpi=300)

    res = pd.DataFrame(summary).sort_values("roc_auc", ascending=False)
    res.to_csv(OUT / "threshold_sensitivity.tsv", sep="\t", index=False)
    print(res.to_string(index=False))
    print(f"\nMean ROC AUC across families: {res['roc_auc'].mean():.3f}")
    print(f"Saved {OUT}/threshold_sensitivity.png and .tsv")


if __name__ == "__main__":
    main()
