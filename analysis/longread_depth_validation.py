#!/usr/bin/env python3
"""
Long-read DEPTH validation of CNVRock chromosomal blaSHV CRR.

This is the orthogonal, real-data validation of CNVRock's copy-number VALUES
(not just presence). Short-read assembly collapses tandem arrays and cannot
serve as ground truth; long-read DEPTH does not — every tandem copy maps to
the same reference position, so long-read coverage at the blaSHV locus scales
with copy number exactly as short-read depth does, but from an independent
platform and pipeline (minimap2 map-ont vs BWA-MEM + Conv-VAE baseline).

Inputs (produced by hpc/longread_ont_depth.sh, one TSV per ONT run):
    data/results/longread_depth/<run>.tsv  ->  sample, run, shv_depth,
                                               chrom_mean, longread_CRR

This script joins those to CNVRock's per-sample short-read crr_blaSHV
(gene_calls.tsv) via sample_accession -> sample_id, and reports:
  - Spearman and Pearson correlation of long-read CRR vs CNVRock CRR
  - concordance of the amplified call (CRR >= 1.75) between the two platforms
  - a scatter plot

    python3 analysis/longread_depth_validation.py
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

REPO       = Path(__file__).resolve().parents[1]
LR_DIR     = REPO / "data/results/longread_depth"
META       = REPO / "assets/kpsc_expansion_metadata_runlevel.tsv"
GENE_CALLS = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT        = REPO / "data/results/longread_depth"
AMP_THRESH = 1.75


def main() -> None:
    parts = []
    for f in sorted(glob.glob(str(LR_DIR / "*.tsv"))):
        if Path(f).name == "longread_depth_joined.tsv":
            continue
        try:
            row = pd.read_csv(f, sep="\t", header=None,
                              names=["sample_accession", "run", "shv_depth",
                                     "chrom_mean", "lr_crr"])
            parts.append(row)
        except Exception:
            continue
    if not parts:
        raise SystemExit("No long-read depth results yet in "
                         f"{LR_DIR} — run hpc/longread_ont_depth.sh first.")
    lr = pd.concat(parts, ignore_index=True)
    lr = lr[pd.to_numeric(lr["lr_crr"], errors="coerce").notna()].copy()
    lr["lr_crr"] = lr["lr_crr"].astype(float)
    print(f"Long-read runs with a CRR: {len(lr)}")

    # map sample_accession -> CNVRock sample_id (run-level)
    meta = pd.read_csv(META, sep="\t", dtype=str)[["sample_id", "sample_accession"]]
    gc = pd.read_csv(GENE_CALLS, sep="\t")[["sample_id", "crr_blaSHV"]]
    j = (lr.merge(meta, on="sample_accession", how="left")
           .merge(gc, on="sample_id", how="inner")
           .dropna(subset=["crr_blaSHV", "lr_crr"]))
    j = j.rename(columns={"crr_blaSHV": "sr_crr"})
    print(f"Joined to CNVRock short-read CRR: {len(j)} isolates")
    if len(j) < 5:
        raise SystemExit("Too few joined isolates for a correlation.")

    rho, p_rho = stats.spearmanr(j["lr_crr"], j["sr_crr"])
    r, p_r     = stats.pearsonr(j["lr_crr"], j["sr_crr"])
    print(f"\nSpearman ρ = {rho:.3f} (p = {p_rho:.2e})")
    print(f"Pearson  r = {r:.3f} (p = {p_r:.2e})")

    # amplified-call concordance
    lr_amp = j["lr_crr"] >= AMP_THRESH
    sr_amp = j["sr_crr"] >= AMP_THRESH
    tp = int((lr_amp & sr_amp).sum()); tn = int((~lr_amp & ~sr_amp).sum())
    fp = int((~lr_amp & sr_amp).sum()); fn = int((lr_amp & ~sr_amp).sum())
    print(f"\nAmplified-call (CRR≥{AMP_THRESH}) concordance vs long-read:")
    print(f"  both amp {tp} | both normal {tn} | SR-only {fp} | LR-only {fn}")
    n_lr_amp = int(lr_amp.sum())
    if n_lr_amp:
        print(f"  of {n_lr_amp} long-read-amplified isolates, "
              f"CNVRock recovers {tp} ({100*tp/n_lr_amp:.0f}%)")

    j.to_csv(OUT / "longread_depth_joined.tsv", sep="\t", index=False)

    # ── scatter ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(j["lr_crr"], j["sr_crr"], s=26, alpha=0.6, color="#1f4e79",
               edgecolors="white", linewidths=0.4)
    lim = max(j["lr_crr"].max(), j["sr_crr"].max()) * 1.05
    ax.plot([0, lim], [0, lim], ls="--", lw=1, color="#888", label="y = x")
    ax.axhline(AMP_THRESH, ls=":", lw=0.8, color=(0.8, 0.2, 0.2))
    ax.axvline(AMP_THRESH, ls=":", lw=0.8, color=(0.8, 0.2, 0.2))
    ax.set_xlabel("Long-read (ONT) blaSHV copy-ratio (minimap2 depth)", fontsize=10)
    ax.set_ylabel("CNVRock short-read blaSHV CRR (Conv-VAE baseline)", fontsize=10)
    ax.set_title(f"Long-read depth validation of blaSHV CRR\n"
                 f"n = {len(j)}   Spearman ρ = {rho:.2f}   Pearson r = {r:.2f}",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.legend(fontsize=9, frameon=False)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT / "longread_depth_validation.png", dpi=300)
    print(f"\nSaved {OUT}/longread_depth_validation.png and joined TSV")


if __name__ == "__main__":
    main()
