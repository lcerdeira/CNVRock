#!/usr/bin/env python3
"""
Porin copy-number-loss analysis (Phase E reviewer response).

A reviewer noted that resistance mechanisms beyond gene amplification —
notably outer-membrane porin loss (ompK35/ompK36) — also raise carbapenem
MICs. Point mutations in these porins are out of CNVRock's scope (it is a
CNV tool, not an SNP caller), but porin gene LOSS / IS-disruption is a
copy-number event CNVRock can see.

This script tests, at the 10K tier, whether the VAE copy-ratio at the
ompK35 and ompK36 loci is LOWER in carbapenem-resistant samples than in
susceptible ones — i.e. whether the genome-wide CNV signal independently
recovers porin-loss as a carbapenem-resistance mechanism.

OmpK35  NC_016845.1:1,904,308-1,905,386  (KPHS_18370/18380) -> bins 1904-1905
OmpK36  NC_016845.1:3,727,882-3,728,985  (KPHS_37010)       -> bins 3727-3728

Run locally:  python3 analysis/porin_cnv_analysis.py
Output:       data/results/porin_cnv/porin_cnv_summary.tsv  + figure
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
STORE = REPO / "data/inputs/KpSC-expansion-10k-mq20-1000bp-npy"
RESULTS = REPO / "data/results/33_kpsc_expansion_10k"
META = REPO / "assets/kpsc_expansion_metadata_runlevel.tsv"
CABBAGE = REPO / "assets/cabbage_kpsc_phenotypes.tsv"
OUT = REPO / "data/results/porin_cnv"

PORINS = {
    "ompK35": list(range(1904, 1906)),   # KPHS_18370/18380
    "ompK36": list(range(3727, 3729)),   # KPHS_37010
}
CARBAPENEMS = ["meropenem", "imipenem", "ertapenem", "doripenem"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # ── observed depth (chromosome) + VAE reconstruction -> copy-ratio ────
    sample_ids = np.array([str(s) for s in
                           np.load(RESULTS / "sample_ids.npy", allow_pickle=True)])
    counts = np.load(STORE / "counts.npy").astype(np.float32)
    store_sids = np.array([str(s) for s in
                           np.load(STORE / "sample_ids.npy", allow_pickle=True)])
    contigs = np.array([str(c) for c in np.load(STORE / "contigs.npy",
                                                allow_pickle=True)])
    m = contigs == "NC_016845.1"
    if m.any():
        counts = counts[:, m]
    row = {s: i for i, s in enumerate(store_sids)}
    idx = [row.get(s, -1) for s in sample_ids]
    x = np.full((len(sample_ids), counts.shape[1]), np.nan, np.float32)
    keep = [i for i, r in enumerate(idx) if r >= 0]
    x[keep] = counts[[idx[i] for i in keep]]
    xhat = np.load(RESULTS / "reconstructions.npy").astype(np.float32)
    nb = min(x.shape[1], xhat.shape[1])
    with np.errstate(invalid="ignore", divide="ignore"):
        r = x[:, :nb] / np.where(xhat[:, :nb] == 0, np.nan, xhat[:, :nb])

    # per-sample porin copy-ratio = median over the porin's bins
    porin_cr = {}
    for name, bins in PORINS.items():
        bins = [b for b in bins if b < nb]
        porin_cr[name] = np.nanmedian(r[:, bins], axis=1)

    # ── carbapenem phenotype from CABBAGE ────────────────────────────────
    meta = pd.read_csv(META, sep="\t")[
        ["sample_id", "sample_accession"]].rename(
        columns={"sample_accession": "biosample"})
    cab = pd.read_csv(CABBAGE, sep="\t",
                      usecols=["BioSample_ID", "antibiotic_name", "ast_standard",
                               "updated_phenotype", "resistance_phenotype"])
    cab["pheno"] = (cab["updated_phenotype"].fillna(cab["resistance_phenotype"])
                    .astype(str).str.lower().str[:1].str.upper())
    cab = cab[cab["pheno"].isin(["R", "S"])]
    ph = cab.rename(columns={"BioSample_ID": "biosample"}).merge(
        meta, on="biosample", how="inner")
    sid_idx = pd.Series(np.arange(len(sample_ids)), index=sample_ids)
    ph["row"] = ph["sample_id"].map(sid_idx)
    ph = ph.dropna(subset=["row"])
    ph["row"] = ph["row"].astype(int)

    # ── test: is porin copy-ratio LOWER in R than S? ─────────────────────
    rows = []
    for porin, cr in porin_cr.items():
        for drug in CARBAPENEMS:
            g = ph[ph["antibiotic_name"] == drug]
            ri = g.loc[g["pheno"] == "R", "row"].values
            si = g.loc[g["pheno"] == "S", "row"].values
            cr_r = cr[ri]; cr_r = cr_r[~np.isnan(cr_r)]
            cr_s = cr[si]; cr_s = cr_s[~np.isnan(cr_s)]
            if len(cr_r) < 8 or len(cr_s) < 8:
                rows.append({"porin": porin, "antibiotic": drug,
                             "n_R": len(cr_r), "n_S": len(cr_s),
                             "cr_median_R": np.nan, "cr_median_S": np.nan,
                             "mwu_p_R_lt_S": np.nan,
                             "frac_R_CNloss": np.nan, "frac_S_CNloss": np.nan})
                continue
            # one-sided: copy-ratio in R is LESS than in S (porin loss)
            _, p = stats.mannwhitneyu(cr_r, cr_s, alternative="less")
            rows.append({
                "porin": porin, "antibiotic": drug,
                "n_R": len(cr_r), "n_S": len(cr_s),
                "cr_median_R": round(float(np.median(cr_r)), 3),
                "cr_median_S": round(float(np.median(cr_s)), 3),
                "mwu_p_R_lt_S": p,
                # CN-loss = copy-ratio below 0.5 (heterozygous-style loss /
                # IS-disruption leaves partial depth)
                "frac_R_CNloss": round(float(np.mean(cr_r < 0.5)), 3),
                "frac_S_CNloss": round(float(np.mean(cr_s < 0.5)), 3),
            })
    summ = pd.DataFrame(rows)
    summ_path = OUT / "porin_cnv_summary.tsv"
    summ.to_csv(summ_path, sep="\t", index=False)
    print(summ.to_string(index=False))
    print(f"\nwrote {summ_path}")

    # ── figure: porin copy-ratio R vs S per carbapenem ───────────────────
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.0))
    rng = np.random.default_rng(7)
    for i, porin in enumerate(PORINS):
        cr = porin_cr[porin]
        for j, drug in enumerate(CARBAPENEMS):
            ax = axes[i, j]
            g = ph[ph["antibiotic_name"] == drug]
            ri = g.loc[g["pheno"] == "R", "row"].values
            si = g.loc[g["pheno"] == "S", "row"].values
            cr_r = cr[ri]; cr_r = cr_r[~np.isnan(cr_r)]
            cr_s = cr[si]; cr_s = cr_s[~np.isnan(cr_s)]
            if len(cr_r) < 8 or len(cr_s) < 8:
                ax.set_title(f"{porin} × {drug}\n(n too small)", fontsize=8)
                ax.set_xticks([])
                continue
            bp = ax.boxplot([np.clip(cr_s, 0, 2), np.clip(cr_r, 0, 2)],
                            positions=[0, 1], widths=0.55, patch_artist=True,
                            showfliers=False,
                            medianprops=dict(color="black", lw=1.4))
            for patch, c in zip(bp["boxes"], ["#7e8a99", "#c1272d"]):
                patch.set_facecolor(c); patch.set_alpha(0.55)
            ax.scatter(rng.uniform(-0.16, 0.16, len(cr_s)),
                       np.clip(cr_s, 0, 2), s=7, color="#3a3a3a",
                       alpha=0.4, edgecolors="none")
            ax.scatter(1 + rng.uniform(-0.16, 0.16, len(cr_r)),
                       np.clip(cr_r, 0, 2), s=7, color="#7a0000",
                       alpha=0.5, edgecolors="none")
            _, p = stats.mannwhitneyu(cr_r, cr_s, alternative="less")
            ax.axhline(0.5, color="#b45f06", ls="--", lw=0.8)
            ptxt = f"p = {p:.1e}" if p < 0.001 else f"p = {p:.3f}"
            ax.set_title(f"{porin} × {drug}\n{ptxt}", fontsize=8,
                         fontweight="bold")
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"S\n(n={len(cr_s)})", f"R\n(n={len(cr_r)})"],
                               fontsize=7.5)
            if j == 0:
                ax.set_ylabel(f"{porin}\ncopy-ratio (capped 2×)", fontsize=8)
            ax.grid(axis="y", lw=0.3, alpha=0.4)
    fig.suptitle("Porin (ompK35 / ompK36) copy-ratio by carbapenem phenotype. "
                 "Dashed line = copy-ratio 0.5 (CN-loss threshold); "
                 "one-sided Mann-Whitney U tests copy-ratio(R) < copy-ratio(S).",
                 fontsize=9.5, y=1.02)
    fig.tight_layout()
    fig_path = OUT / "porin_cnv_figure.png"
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
