#!/usr/bin/env python3
"""
Per-call uncertainty for the chromosomal blaSHV copy-ratio (A7).

Monte-Carlo dropout perturbs the VAE reconstruction (the expected-depth
baseline), but the chromosomal gene call is a DEPTH ratio (gene vs flank)
that does not use the reconstruction — so dropout does not yield a clean
per-call confidence interval. The dominant source of per-call uncertainty
is instead the Poisson sampling noise of the read depth over the ~2 kb
blaSHV locus. We quantify it directly with a parametric bootstrap on the
raw 1 kb read counts (mq ≥ 20 chromosomal store):

  CRR = mean(depth over blaSHV bins) / median(depth over local flank bins)

Per sample, 1 000 bootstrap replicates draw the gene read count from a
Poisson and resample the flank bins, giving a 95 % CI on the CRR. This is a
genuine per-call confidence interval usable for clinical triage: a call
whose CI excludes 1.0 is a confident amplification.

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
GCALLS = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT    = REPO / "data/results/percall_uncertainty"
CHROM  = "NC_016845.1"
SHV_START, SHV_END = 2549403, 2550263          # blaSHV CDS
FLANK_HALF = 100_000                            # ±100 kb local flank
N_BOOT = 1000
AMP = 1.75


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts  = np.load(STORE / "counts.npy")                       # (n, L) uint32
    ids     = np.load(STORE / "sample_ids.npy", allow_pickle=True)
    contigs = np.load(STORE / "contigs.npy", allow_pickle=True)

    def in_range(lo, hi):
        return np.array([i for i, c in enumerate(contigs)
                         if str(c[0]) == CHROM and int(c[1]) >= lo and int(c[2]) <= hi])
    gene_bins  = in_range(SHV_START - 999, SHV_END + 999)          # bins over the CDS
    flank_bins = np.array([b for b in in_range(SHV_START - FLANK_HALF, SHV_END + FLANK_HALF)
                           if b not in set(gene_bins.tolist())])
    print(f"blaSHV gene bins: {list(gene_bins)} | flank bins: {len(flank_bins)}")

    rng = np.random.default_rng(42)
    g_reads = counts[:, gene_bins].sum(axis=1).astype(float)       # total gene reads
    n_gbin  = len(gene_bins)
    flank_c = counts[:, flank_bins].astype(float)
    flank_med = np.median(flank_c, axis=1)
    crr = (g_reads / n_gbin) / np.where(flank_med > 0, flank_med, np.nan)

    lo = np.full(len(ids), np.nan); hi = np.full(len(ids), np.nan)
    for i in range(len(ids)):
        fm = flank_med[i]
        if not np.isfinite(crr[i]) or fm <= 0:
            continue
        gb = rng.poisson(g_reads[i], N_BOOT) / n_gbin              # Poisson gene depth
        fb = np.median(rng.choice(flank_c[i], size=(N_BOOT, flank_c.shape[1]),
                                  replace=True), axis=1)           # resampled flank
        fb = np.where(fb > 0, fb, np.nan)
        b = gb / fb
        lo[i], hi[i] = np.nanpercentile(b, [2.5, 97.5])

    d = pd.DataFrame({"sample_id": ids, "shv_crr": crr,
                      "ci_lo": lo, "ci_hi": hi})
    d["ci_width"] = d["ci_hi"] - d["ci_lo"]
    d = d[np.isfinite(d["shv_crr"]) & np.isfinite(d["ci_width"])]
    d.to_csv(OUT / "blashv_percall_ci.tsv", sep="\t", index=False)

    amp = d[d["shv_crr"] >= AMP]
    # a call is a "confident amplification" when its 95% CI excludes 1.0
    conf = amp[amp["ci_lo"] > 1.0]
    print(f"\nn samples: {len(d)}")
    print(f"median CRR MC-CI width: {d['ci_width'].median():.3f}")
    print(f"amplified calls (CRR ≥ {AMP}): {len(amp)}")
    print(f"  of which CI excludes 1.0 (confident amplification): "
          f"{len(conf)} ({100*len(conf)/max(len(amp),1):.0f}%)")
    print(f"  median CI width, amplified: {amp['ci_width'].median():.3f}")

    # cross-check bootstrap CRR against the pipeline's crr_blaSHV
    try:
        gc = pd.read_csv(GCALLS, sep="\t")[["sample_id", "crr_blaSHV"]]
        m = d.merge(gc, on="sample_id", how="inner").dropna()
        from scipy import stats
        rho, _ = stats.spearmanr(m["shv_crr"], m["crr_blaSHV"])
        print(f"\nbootstrap CRR vs pipeline crr_blaSHV: Spearman ρ = {rho:.2f} "
              f"(n = {len(m)})  [sanity check]")
    except Exception as e:
        print("cross-check skipped:", e)

    # ── figure: CRR with 95% CI, ranked (amplified region) ──────────────
    top = d.sort_values("shv_crr", ascending=False).head(60).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    x = np.arange(len(top))
    ax.errorbar(x, top["shv_crr"],
                yerr=[top["shv_crr"] - top["ci_lo"], top["ci_hi"] - top["shv_crr"]],
                fmt="o", ms=3, lw=0.8, color="#1f4e79", ecolor="#9bb8d3",
                capsize=1.5)
    ax.axhline(AMP, ls="--", color="#b00", lw=1, label=f"amplified ≥ {AMP}")
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
