#!/usr/bin/env python3
"""
Temporal trend in AMR gene amplification (limitation 2, first-pass).

The per-gene λ/μ ratios (§4) are steady-state estimates. Here we ask the
downstream surveillance question directly: is the *prevalence of
amplification among carriers* changing over calendar time? For each gene we
fit a logistic regression of the amplified indicator on the (centred)
collection year among carriers, and report the per-year odds ratio with a
2 000-resample bootstrap 95% CI. Collection year is recovered from ENA
sample metadata (7 448/10 000 isolates have a usable year, 2013–2024).

This is a first-pass, calendar-time trend and is NOT adjusted for the
sequence-type / geography confounding that also affects the per-ST analysis
(§3.6): which clones and countries were sampled shifts across years, so a
trend may reflect ascertainment as well as biology. A fully confounder-
adjusted hierarchical model remains future work.

    python3 analysis/temporal_amplification.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO  = Path(__file__).resolve().parents[1]
YC    = REPO / "data/results/expansion_10k_year_country.tsv"
PLASM = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
CHROM = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT   = REPO / "data/results/temporal_amplification"
N_BOOT = 2000

GENES = {
    "blaKPC":      (["pcn_blaKPC-2"], 0.20, 1.50),
    "blaCTX-M":    (["pcn_blaCTX-M-15", "pcn_blaCTX-M-14",
                     "pcn_blaCTX-M-27", "pcn_blaCTX-M-65"], 0.20, 1.50),
    "blaSHV(chr)": (["crr_blaSHV"], None, 1.75),   # chromosomal CRR
}


def logistic_fit(x, y, iters=200, lr=0.3):
    """Tiny gradient-descent logistic (intercept + slope); standardise x."""
    xs = (x - x.mean()) / (x.std() + 1e-9)
    b0, b1 = 0.0, 0.0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(b0 + b1 * xs)))
        g0 = np.mean(p - y); g1 = np.mean((p - y) * xs)
        b0 -= lr * g0; b1 -= lr * g1
    return b1 / (x.std() + 1e-9)          # slope on the ORIGINAL year scale


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    yc = pd.read_csv(YC, sep="\t").dropna(subset=["year"])
    yc["year"] = yc["year"].astype(int)
    yc = yc[(yc.year >= 2013) & (yc.year <= 2024)]
    plasm = pd.read_csv(PLASM, sep="\t")
    chrom = pd.read_csv(CHROM, sep="\t")[["sample_id", "crr_blaSHV"]]
    d = yc.merge(plasm, on="sample_id").merge(chrom, on="sample_id", how="left")

    rng = np.random.default_rng(42)
    rows = []
    fig, ax = plt.subplots(figsize=(8, 5))
    for gene, (cols, absent, amp_t) in GENES.items():
        s = d[cols].fillna(0).sum(axis=1) if len(cols) > 1 else d[cols[0]].fillna(0)
        carrier = (s > absent) if absent is not None else d["crr_blaSHV"].notna()
        sub = pd.DataFrame({"year": d["year"], "amp": (s >= amp_t).astype(int)})[carrier].dropna()
        if len(sub) < 50:
            continue
        yr = sub["year"].values.astype(float); am = sub["amp"].values.astype(float)
        slope = logistic_fit(yr, am)
        boot = np.array([logistic_fit(*(lambda idx: (yr[idx], am[idx]))(
            rng.integers(0, len(yr), len(yr)))) for _ in range(N_BOOT)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        or_yr = np.exp(slope)
        rows.append(dict(gene=gene, n_carriers=len(sub), n_amp=int(am.sum()),
                         or_per_year=round(or_yr, 3),
                         ci_lo=round(np.exp(lo), 3), ci_hi=round(np.exp(hi), 3),
                         ci_excludes_1=bool(lo > 0 or hi < 0)))
        # prevalence-by-year for the plot
        prev = sub.groupby("year")["amp"].agg(["mean", "size"])
        prev = prev[prev["size"] >= 8]
        ax.plot(prev.index, prev["mean"] * 100, marker="o", ms=4, lw=1.5,
                label=f"{gene} (OR/yr {or_yr:.2f})")

    ax.set_xlabel("collection year", fontsize=10)
    ax.set_ylabel("amplification prevalence among carriers (%)", fontsize=10)
    ax.set_title("Temporal trend in AMR gene amplification (first-pass, "
                 "unadjusted)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT / "temporal_amplification.png", dpi=300)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "temporal_amplification.tsv", sep="\t", index=False)
    print(res.to_string(index=False))
    print(f"\nSaved {OUT}/temporal_amplification.png and .tsv")


if __name__ == "__main__":
    main()
