#!/usr/bin/env python3
"""
VAE-baseline ablation (Phase E reviewer response) — v2.

A reviewer asked whether the convolutional-VAE expected-depth baseline
actually beats the simple normalisation used by targeted estimators such
as CCNE (Jiang et al. 2022). We compare three ways of turning observed
1 kb depth x into a copy-ratio r:

  A. SIMPLE        r = x / median(x)              per-sample scalar median
  B. HOUSEKEEPING  r = x / mean(x at MLST genes)  CCNE-style single-locus
  C. VAE           r = x / x_hat                  CNVRock learned baseline

Three rigorous, non-circular metrics:

  1. NOISE FLOOR — MAD of r across the genome (lower = the baseline
     removed more technical depth structure; r sits tighter around 1).

  2. SPIKE-IN RECOVERY RMSE — inject a known copy number f∈{1.5,2,3,4}
     into random 30-bin regions of random samples (the validation design
     CCNE itself used), recompute r, and measure how accurately each
     baseline recovers f. Directly comparable to CCNE's reported RMSE.

  3. FALSE-AMPLIFICATION RATE — in NON-spiked control bins, the fraction
     called amplified (r > 1.75). A noisy baseline manufactures spurious
     amplifications; lower is better (specificity).

Run on HPC:  sbatch hpc/run_vae_ablation.sh
Output:      data/results/vae_ablation/ablation_summary.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# K. pneumoniae MLST housekeeping set (Institut Pasteur scheme), HS11286
# NC_016845.1 CDS-midpoint positions (bp).
#
# CORRECTED 2026-07-22. The previous values were labelled as the MLST set but
# were wrong by megabases (gapA was off by ~1 Mb, infB by ~3.5 Mb, rpoB by
# ~4.6 Mb), so baseline B was in fact 7 arbitrary chromosomal bins rather than
# the housekeeping set it claimed to be. Positions below were read from the
# HS11286 feature table by product annotation, since that assembly uses
# KPHS_* locus tags rather than gene names.
#
# tonB is omitted: the HS11286 annotation contains only "TonB-dependent
# receptor" entries, which are different genes, and no unambiguous tonB. Six
# loci are used rather than guessing a seventh.
MLST_GENES = {
    "gapA": 2133023, "infB": 4733635, "mdh": 2815555,
    "pgi":  288228,  "phoE": 1076695, "rpoB": 229368,
}
AMP_THRESHOLD = 1.75


def per_sample_median_norm(x):
    med = np.nanmedian(x, axis=1, keepdims=True)
    return x / np.where(med == 0, np.nan, med)


def housekeeping_ratio(x, bin_size=1000):
    hk_bins = [int(p // bin_size) for p in MLST_GENES.values()]
    hk_bins = [b for b in hk_bins if 0 <= b < x.shape[1]]
    hk = np.nanmean(x[:, hk_bins], axis=1, keepdims=True)
    return x / np.where(hk == 0, np.nan, hk)


def noise_floor(r):
    """Median per-sample MAD of r around its own median — technical noise."""
    mad = np.nanmedian(np.abs(r - np.nanmedian(r, axis=1, keepdims=True)),
                       axis=1)
    return float(np.nanmedian(mad))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store-dir", type=Path,
                    default=Path("data/inputs/KpSC-expansion-10k-mq20-1000bp-npy"))
    ap.add_argument("--results-dir", type=Path,
                    default=Path("data/results/33_kpsc_expansion_10k"))
    ap.add_argument("--n-spikes", type=int, default=400,
                    help="number of (sample,region) spike-in tests per CN level")
    ap.add_argument("--region-bins", type=int, default=30,
                    help="width of each injected amplification, in 1 kb bins")
    ap.add_argument("--out-dir", type=Path,
                    default=Path("data/results/vae_ablation"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ── load observed depth (chromosome bins) + VAE reconstruction ────────
    sample_ids = np.array([str(s) for s in
                           np.load(args.results_dir / "sample_ids.npy",
                                   allow_pickle=True)])
    counts = np.load(args.store_dir / "counts.npy").astype(np.float32)
    store_sids = np.array([str(s) for s in
                           np.load(args.store_dir / "sample_ids.npy",
                                   allow_pickle=True)])
    contigs_p = args.store_dir / "contigs.npy"
    if contigs_p.exists():
        contigs = np.array([str(c) for c in np.load(contigs_p, allow_pickle=True)])
        m = contigs == "NC_016845.1"
        if m.any():
            counts = counts[:, m]
    row = {s: i for i, s in enumerate(store_sids)}
    idx = [row.get(s, -1) for s in sample_ids]
    x = np.full((len(sample_ids), counts.shape[1]), np.nan, np.float32)
    keep = [i for i, r in enumerate(idx) if r >= 0]
    x[keep] = counts[[idx[i] for i in keep]]
    xhat = np.load(args.results_dir / "reconstructions.npy").astype(np.float32)
    nb = min(x.shape[1], xhat.shape[1])
    x, xhat = x[:, :nb], xhat[:, :nb]
    print(f"observed depth {x.shape}, matched {len(keep)}/{len(sample_ids)}")

    def ratios(obs):
        with np.errstate(invalid="ignore", divide="ignore"):
            return {"A_simple_median": per_sample_median_norm(obs),
                    "B_housekeeping":  housekeeping_ratio(obs),
                    "C_vae":           obs / np.where(xhat == 0, np.nan, xhat)}

    # ── metric 1: noise floor on un-spiked data ───────────────────────────
    base_r = ratios(x)
    nf = {k: noise_floor(v) for k, v in base_r.items()}

    # ── metric 3: false-amplification rate on un-spiked data ──────────────
    # un-spiked genome should be CN≈1 almost everywhere; any bin>1.75 is a
    # baseline-induced spurious amplification.
    fa = {k: float(np.nanmean(v > AMP_THRESHOLD)) for k, v in base_r.items()}

    # ── metric 2: spike-in recovery RMSE ──────────────────────────────────
    cn_levels = [1.5, 2.0, 3.0, 4.0]
    valid_rows = np.array(keep)
    w = args.region_bins
    # avoid the MLST housekeeping bins and edges when placing spikes
    hk_bins = set(int(p // 1000) for p in MLST_GENES.values())
    spike_recovery = {k: [] for k in base_r}
    spike_truth = []
    for f in cn_levels:
        for _ in range(args.n_spikes):
            s = int(rng.choice(valid_rows))
            start = int(rng.integers(w, nb - w))
            if any(b in hk_bins for b in range(start, start + w)):
                continue
            obs = x.copy()
            obs[s, start:start + w] *= f
            rr = ratios(obs)
            spike_truth.append(f)
            for k in base_r:
                rec = np.nanmedian(rr[k][s, start:start + w])
                spike_recovery[k].append(rec)
    spike_truth = np.array(spike_truth)
    rmse = {}
    for k in base_r:
        rec = np.array(spike_recovery[k])
        ok = ~np.isnan(rec)
        rmse[k] = float(np.sqrt(np.mean((rec[ok] - spike_truth[ok]) ** 2)))

    # ── assemble ──────────────────────────────────────────────────────────
    rows = []
    for k in base_r:
        rows.append({"baseline": k,
                     "noise_floor_MAD": round(nf[k], 4),
                     "spikein_recovery_RMSE": round(rmse[k], 4),
                     "false_amp_rate": round(fa[k], 5)})
    out = pd.DataFrame(rows)
    out_path = args.out_dir / "ablation_summary.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nwrote {out_path}")
    print(out.to_string(index=False))
    # winner per metric (lower is better for all three)
    print("\nbest per metric (lower=better):")
    for col in ["noise_floor_MAD", "spikein_recovery_RMSE", "false_amp_rate"]:
        best = out.loc[out[col].idxmin(), "baseline"]
        print(f"  {col:24s} -> {best}")


if __name__ == "__main__":
    main()
