#!/usr/bin/env python3
"""
Per-call uncertainty for the chromosomal blaSHV copy-ratio (A7).

This bootstraps THE ESTIMATOR THE PIPELINE ACTUALLY CALLS. An earlier version
of this script bootstrapped a different quantity — raw gene depth over the
median of a *local* +/-100 kb window — on the mistaken premise that the
chromosomal call does not use the VAE reconstruction. It does. The caller
(models/cnv/06_gene_cnv_caller.py) computes a *double* ratio:

    copy_ratio[i,b] = counts[i,b] / safe_recon[i,b]      # VAE-normalised bin
    CRR[i]          = mean(copy_ratio over gene bins)
                    / mean(copy_ratio over flank bins)

where safe_recon falls back to the sample mean wherever the reconstruction is
below hmm_low_cov_threshold, and the flank is every chromosomal bin OUTSIDE
gene +/- cnv_flank_padding (i.e. the rest of the chromosome, not a local
window). Bootstrapping the wrong estimator put 87 of the 171 published
amplification calls outside the analysis entirely.

Uncertainty model. The gene sits on ~2 kb, so its copy-ratio is dominated by
Poisson sampling of the reads in those few bins; we resample them exactly.
The flank spans thousands of bins, so by the CLT its mean is Gaussian with a
variance we propagate analytically rather than resampling (which would cost
n_samples x n_boot x n_flank draws for a negligible correction).

A call whose 95 % CI excludes 1.0 is a confident amplification.

    python3 analysis/percall_uncertainty.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
STORE  = REPO / "data/inputs/KpSC-expansion-10k-mq20-1000bp-npy"
RESULT = REPO / "data/results/33_kpsc_expansion_10k"
OUT    = REPO / "data/results/percall_uncertainty"

CHROM              = "NC_016845.1"
SHV_START, SHV_END = 2549403, 2550263      # blaSHV CDS (hardcoded in the caller)
FLANK_PADDING      = 100_000               # cfg cnv_flank_padding
LOW_COV            = 10                    # cfg hmm_low_cov_threshold
AMP                = 1.75                  # cfg cnv_crr_amp_threshold
N_BOOT             = 1000
SEED               = 42


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # Chromosomal blaSHV is only meaningful for K. pneumoniae /
    # K. quasipneumoniae: K. variicola carries blaLEN at the syntenic locus and
    # its reads cross-map onto the HS11286 blaSHV coordinates (see
    # analysis/kvariicola_multiref.py). The manuscript restricts the call to
    # those two species, so this analysis must use the same scope.
    meta = pd.read_csv(REPO / "assets/kpsc_expansion_metadata_runlevel.tsv",
                       sep="\t")[["sample_id", "Species"]]
    eligible = set(meta.loc[meta["Species"].astype(str)
                   .str.contains("pneumoniae|quasipneumoniae", na=False),
                   "sample_id"])

    counts  = np.load(STORE / "counts.npy").astype(np.float64)
    contigs = pd.DataFrame(np.load(STORE / "contigs.npy", allow_pickle=True))
    contigs.columns = ["chrom", "start", "end"][:contigs.shape[1]]
    recons  = np.load(RESULT / "reconstructions.npy").astype(np.float64)
    ids     = np.load(RESULT / "sample_ids.npy", allow_pickle=True)
    assert counts.shape == recons.shape, (counts.shape, recons.shape)

    # ── reproduce the caller's normalisation exactly ────────────────────────
    sample_mean = counts.mean(axis=1, keepdims=True)
    safe_recon  = np.where(recons >= LOW_COV, recons, sample_mean)
    copy_ratio  = counts / (safe_recon + 1e-6)

    chroms = contigs["chrom"].values.astype(str)
    starts = contigs["start"].values.astype(float)
    chrom_mask = chroms == CHROM
    s          = starts[chrom_mask]
    gene_mask  = (s >= SHV_START) & (s <= SHV_END)
    flank_mask = (s < SHV_START - FLANK_PADDING) | (s > SHV_END + FLANK_PADDING)
    print(f"gene bins: {int(gene_mask.sum())} | flank bins: {int(flank_mask.sum())} "
          f"| chromosomal bins: {int(chrom_mask.sum())}")

    cr_chrom = copy_ratio[:, chrom_mask]
    cnt_c    = counts[:, chrom_mask]
    rec_c    = safe_recon[:, chrom_mask]

    mean_gene  = np.nanmean(cr_chrom[:, gene_mask],  axis=1)
    mean_flank = np.nanmean(cr_chrom[:, flank_mask], axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        crr = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)

    # ── uncertainty ─────────────────────────────────────────────────────────
    # gene: exact Poisson resampling of the few gene-bin counts
    g_cnt = cnt_c[:, gene_mask]                       # (n, G)
    g_rec = rec_c[:, gene_mask] + 1e-6
    G     = g_cnt.shape[1]
    # flank: Gaussian by CLT.  var(mean) = (1/K^2) * sum(lambda_j / r_j^2)
    f_cnt = cnt_c[:, flank_mask]
    f_rec = rec_c[:, flank_mask] + 1e-6
    K     = f_cnt.shape[1]
    flank_sd = np.sqrt(np.nansum(f_cnt / f_rec**2, axis=1)) / K

    rng = np.random.default_rng(SEED)
    lo = np.full(len(ids), np.nan)
    hi = np.full(len(ids), np.nan)
    ok = np.isfinite(crr) & (mean_flank > 0)
    for i in np.flatnonzero(ok):
        gb = rng.poisson(g_cnt[i], size=(N_BOOT, G)) / g_rec[i]     # (B, G)
        gb = gb.mean(axis=1)
        fb = rng.normal(mean_flank[i], flank_sd[i], size=N_BOOT)
        fb = np.where(fb > 0, fb, np.nan)
        b  = gb / fb
        lo[i], hi[i] = np.nanpercentile(b, [2.5, 97.5])

    d = pd.DataFrame({"sample_id": ids, "shv_crr": crr, "ci_lo": lo, "ci_hi": hi})
    d["ci_width"] = d["ci_hi"] - d["ci_lo"]
    d = d[np.isfinite(d["shv_crr"]) & np.isfinite(d["ci_width"])]
    n_all = len(d)
    d = d[d["sample_id"].isin(eligible)].reset_index(drop=True)
    print(f"species scope: {len(d)} of {n_all} samples are "
          f"K. pneumoniae / K. quasipneumoniae")
    d.to_csv(OUT / "blashv_percall_ci.tsv", sep="\t", index=False)

    amp  = d[d["shv_crr"] >= AMP]
    conf = amp[amp["ci_lo"] > 1.0]
    print(f"\nn samples with a CI: {len(d)}")
    print(f"median CI width (all): {d['ci_width'].median():.3f}")
    print(f"amplified calls (CRR >= {AMP}): {len(amp)}")
    print(f"  of which CI excludes 1.0: {len(conf)} "
          f"({100*len(conf)/max(len(amp),1):.0f} %)")
    print(f"  median CI width, amplified: {amp['ci_width'].median():.3f}")

    # ── agreement with the pipeline's own call (should now be exact) ────────
    gc = pd.read_csv(RESULT / "gene_calls.tsv", sep="\t")[["sample_id", "crr_blaSHV"]]
    m  = d.merge(gc, on="sample_id", how="inner").dropna(subset=["crr_blaSHV"])
    delta = (m["shv_crr"] - m["crr_blaSHV"]).abs()
    n_amp_pipe = int((m["crr_blaSHV"] >= AMP).sum())
    print(f"\nagreement with pipeline crr_blaSHV (n = {len(m)}):")
    print(f"  max |difference| = {delta.max():.2e}   median = {delta.median():.2e}")
    print(f"  amplified: pipeline {n_amp_pipe}  vs  this script {len(amp)}")

    # ── figure ──────────────────────────────────────────────────────────────
    top = d.sort_values("shv_crr", ascending=False).head(60).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    x = np.arange(len(top))
    ax.errorbar(x, top["shv_crr"],
                yerr=[top["shv_crr"] - top["ci_lo"], top["ci_hi"] - top["shv_crr"]],
                fmt="o", ms=3, lw=0.8, color="#1f4e79", ecolor="#9bb8d3", capsize=1.5)
    ax.axhline(AMP, ls="--", color="#b00", lw=1, label=f"amplified >= {AMP}")
    ax.axhline(1.0, ls=":", color="#888", lw=1, label="single copy")
    ax.set_xlabel("isolate (top 60 by blaSHV CRR)", fontsize=10)
    ax.set_ylabel("blaSHV copy-ratio  (95 % bootstrap CI)", fontsize=10)
    ax.set_title("Per-call blaSHV CRR with depth-bootstrap confidence intervals",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT / "percall_uncertainty.png", dpi=300)
    print(f"\nSaved {OUT}/blashv_percall_ci.tsv and percall_uncertainty.png")


if __name__ == "__main__":
    main()
