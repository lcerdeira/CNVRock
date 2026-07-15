#!/usr/bin/env python3
"""
Multi-reference chromosomal CNV calling for K. variicola blaLEN (limitation 5a).

Chromosomal blaSHV is called only for K. pneumoniae / K. quasipneumoniae
because K. variicola carries the LEN-family homolog (blaLEN) at the syntenic
locus; its reads cross-map onto the HS11286 blaSHV coordinates and produce
spurious calls (§3.3). We demonstrate a species-appropriate reference lifts
this restriction: the retained short reads of the 10 K-tier K. variicola
isolates were re-aligned to a K. variicola reference (NC_011283.1;
hpc/kvariicola_blalen.sh) and the blaLEN copy-ratio computed from depth.
The blaSHV CDS maps to two paralogous LEN loci, so the family copy-ratio is
the summed depth of both loci over the chromosomal mean.

This script compares the recovered blaLEN CRR against the spurious blaSHV CRR
the same isolates receive on the HS11286 reference.

    python3 analysis/kvariicola_multiref.py
"""
from __future__ import annotations
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
LR_DIR = REPO / "data/results/kvariicola_multiref"
GCALLS = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
AMP = 1.75


def main() -> None:
    parts = []
    for f in sorted(glob.glob(str(LR_DIR / "*.tsv"))):
        if Path(f).name in ("blalen_crr.tsv", "_all.tsv"):
            continue
        parts.append(pd.read_csv(f, sep="\t", header=None,
                                 names=["sample_id", "len1", "len2",
                                        "chrom_mean", "crr_percall"]))
    if not parts:
        raise SystemExit("No K. variicola results — run hpc/kvariicola_blalen.sh first.")
    d = pd.concat(parts, ignore_index=True)
    d["blaLEN_crr"] = (d["len1"] + d["len2"]) / d["chrom_mean"]      # family (2 paralogs)
    d = d[np.isfinite(d["blaLEN_crr"])].copy()
    print(f"K. variicola isolates with a blaLEN call: {len(d)}")
    print(f"  family CRR median {d['blaLEN_crr'].median():.2f}, "
          f"single-copy (0.7–1.3): {((d.blaLEN_crr>=0.7)&(d.blaLEN_crr<=1.3)).sum()}, "
          f"amplified (>{AMP}): {(d.blaLEN_crr>AMP).sum()}")

    gc = pd.read_csv(GCALLS, sep="\t")[["sample_id", "crr_blaSHV"]]
    m = d.merge(gc, on="sample_id", how="inner").dropna()
    rho, _ = stats.spearmanr(m["blaLEN_crr"], m["crr_blaSHV"])
    print(f"\nSame {len(m)} isolates on HS11286: spurious blaSHV CRR median "
          f"{m['crr_blaSHV'].median():.2f}; ρ(blaLEN vs blaSHV) = {rho:.2f}")

    d.to_csv(LR_DIR / "blalen_crr.tsv", sep="\t", index=False)

    # ── figure: (A) blaLEN CRR histogram; (B) blaLEN vs spurious blaSHV ──
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))
    axA.hist(d["blaLEN_crr"].clip(upper=4), bins=40, color="#2ca02c", alpha=0.8)
    axA.axvline(1.0, ls=":", color="#555", label="single copy")
    axA.axvline(AMP, ls="--", color="#b00", label=f"amplified ≥ {AMP}")
    axA.set_xlabel("blaLEN family copy-ratio (K. variicola reference)", fontsize=10)
    axA.set_ylabel("K. variicola isolates", fontsize=10)
    axA.set_title(f"(A) blaLEN is callable with a species reference (n = {len(d)})",
                  fontsize=10, fontweight="bold")
    axA.legend(fontsize=8, frameon=False)

    axB.scatter(m["crr_blaSHV"], m["blaLEN_crr"], s=22, alpha=0.6,
                color="#1f4e79", edgecolors="white", linewidths=0.4)
    lim = max(m["crr_blaSHV"].max(), m["blaLEN_crr"].max()) * 1.05
    axB.plot([0, lim], [0, lim], ls="--", lw=1, color="#888", label="y = x")
    axB.axhline(AMP, ls=":", lw=0.8, color=(0.7, 0.2, 0.2))
    axB.axvline(AMP, ls=":", lw=0.8, color=(0.7, 0.2, 0.2))
    axB.set_xlabel("spurious blaSHV CRR (HS11286 reference)", fontsize=10)
    axB.set_ylabel("blaLEN CRR (K. variicola reference)", fontsize=10)
    axB.set_title(f"(B) HS11286 blaSHV is a distorted proxy for blaLEN (ρ = {rho:.2f})",
                  fontsize=10, fontweight="bold")
    axB.legend(fontsize=8, frameon=False)
    for ax in (axA, axB):
        ax.grid(lw=0.3, alpha=0.4)
    fig.suptitle("Multi-reference chromosomal calling resolves the K. variicola "
                 "species restriction", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(LR_DIR / "kvariicola_multiref.png", dpi=300)
    print(f"\nSaved {LR_DIR}/kvariicola_multiref.png and blalen_crr.tsv")


if __name__ == "__main__":
    main()
