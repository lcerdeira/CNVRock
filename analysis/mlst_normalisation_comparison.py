#!/usr/bin/env python3
"""
Reviewer-proposed MLST-locus normalisation vs the learned VAE baseline.

A reviewer proposed replacing the Conv-VAE expected-depth baseline with a
simpler reference: take the 7 MLST loci, compute their median read depth, and
normalise the locus of interest against that value. This script implements
that protocol exactly and scores it against the pipeline's own baseline on the
same data, so the comparison is a measurement rather than an argument.

Test organism is S. aureus (exp 47) — one of the two organisms the reviewer
named, and the one where we have a ground truth (AMRFinder+ mecA presence,
1,698 MRSA / 1,264 MSSA out of 2,962 isolates).

Two arms, identical inputs, identical gene coordinates:

  VAE  (pipeline)   CRR = mean_g(counts/recon) / mean_f(counts/recon)
                    where g = gene bins, f = chromosomal flank bins.

  MLST (reviewer)   CRR = mean_g(counts) / median(counts at the 7 MLST bins)

Scored on:
  1. mecA presence detection (MCC / FNR / PPV vs AMRFinder+), and the AUC,
     which is threshold-free and so does not favour either arm's scale.
  2. Noise floor: per-sample MAD of the normalised profile. This is the
     quantity that determines whether subtle (1.5-2x) CNV is resolvable.

    python3 analysis/mlst_normalisation_comparison.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
STORE  = REPO / "data/inputs/saureus-USA300-mq20-1000bp-npy"
RESULT = REPO / "data/results/47_saureus"
GT     = REPO / "assets/amrfinder_gt_saureus.tsv"

CHROM      = "NC_007793.1"
MECA_START, MECA_END = 39127, 41136        # SCCmec-borne mecA
FLANK_PAD  = 100_000                        # cfg cnv_flank_padding
LOW_COV    = 10                             # cfg hmm_low_cov_threshold

# S. aureus 7-gene MLST scheme, coordinates read from the USA300 FPR3757
# (NC_007793.1) annotation. glpF is annotated as "MIP/aquaporin family
# protein" rather than by gene name.
MLST_LOCI = {
    "arcC":   (68890,   69819),
    "yqiL":   (268835,  270019),
    "pta":    (644279,  645265),
    "tpi":    (845808,  846569),
    "gmk":    (1206797, 1207420),
    "glpF":   (1312304, 1313122),
    "aroE":   (1704228, 1705034),
}


def mad(profile: np.ndarray) -> float:
    """Median per-sample MAD around each sample's own median."""
    med = np.nanmedian(profile, axis=1, keepdims=True)
    return float(np.nanmedian(np.nanmedian(np.abs(profile - med), axis=1)))


def metrics(truth: np.ndarray, score: np.ndarray, thr: float) -> dict:
    from sklearn.metrics import matthews_corrcoef, roc_auc_score
    ok = np.isfinite(score)
    y, s = truth[ok], score[ok]
    pred = (s >= thr).astype(int)
    tp = int(((y == 1) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    return {
        "auc": round(float(roc_auc_score(y, s)), 4),
        "mcc": round(float(matthews_corrcoef(y, pred)), 3),
        "fnr": round(fn / max(tp + fn, 1), 3),
        "ppv": round(tp / max(tp + fp, 1), 3),
        "n":   int(ok.sum()),
    }


def best_mcc(truth: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    """Best achievable MCC over a threshold sweep — gives each arm its own
    optimum so the comparison does not hinge on a threshold chosen for one."""
    from sklearn.metrics import matthews_corrcoef
    ok = np.isfinite(score)
    y, s = truth[ok], score[ok]
    grid = np.quantile(s[np.isfinite(s)], np.linspace(0.01, 0.99, 199))
    best = max(((matthews_corrcoef(y, (s >= t).astype(int)), t) for t in grid),
               key=lambda p: p[0])
    return round(float(best[0]), 3), round(float(best[1]), 3)


def main() -> None:
    counts  = np.load(STORE / "counts.npy").astype(np.float64)
    contigs = pd.DataFrame(np.load(STORE / "contigs.npy", allow_pickle=True))
    contigs.columns = ["chrom", "start", "end"][:contigs.shape[1]]
    recons  = np.load(RESULT / "reconstructions.npy").astype(np.float64)
    ids     = np.load(RESULT / "sample_ids.npy", allow_pickle=True)

    chroms = contigs["chrom"].values.astype(str)
    starts = contigs["start"].values.astype(float)
    cmask  = chroms == CHROM
    s      = starts[cmask]
    cnt_c  = counts[:, cmask]

    gene_mask  = (s >= MECA_START) & (s <= MECA_END)
    flank_mask = (s < MECA_START - FLANK_PAD) | (s > MECA_END + FLANK_PAD)

    # ── arm 1: the pipeline's VAE baseline ──────────────────────────────────
    sample_mean = counts.mean(axis=1, keepdims=True)
    safe_recon  = np.where(recons >= LOW_COV, recons, sample_mean)
    cr          = (counts / (safe_recon + 1e-6))[:, cmask]
    vae_crr = np.nanmean(cr[:, gene_mask], axis=1) / np.nanmean(cr[:, flank_mask], axis=1)
    vae_profile = cr / np.nanmedian(cr, axis=1, keepdims=True)

    # ── arm 2: the reviewer's 7-locus MLST median ───────────────────────────
    mlst_bins = []
    for name, (lo, hi) in MLST_LOCI.items():
        idx = np.flatnonzero((s >= lo - 1000) & (s <= hi + 1000))
        if len(idx) == 0:
            print(f"  WARNING: no bin for MLST locus {name}")
            continue
        mlst_bins.append(idx)
    per_locus = np.stack([cnt_c[:, idx].mean(axis=1) for idx in mlst_bins], axis=1)
    mlst_ref  = np.median(per_locus, axis=1)          # median across the 7 loci
    print(f"MLST reference: {len(mlst_bins)} loci, "
          f"median depth {np.median(mlst_ref):.0f}x")

    denom = np.where(mlst_ref > 0, mlst_ref, np.nan)
    mlst_crr = np.nanmean(cnt_c[:, gene_mask], axis=1) / denom
    mlst_profile = cnt_c / denom[:, None]

    # ── ground truth ────────────────────────────────────────────────────────
    gt = pd.read_csv(GT, sep="\t")[["sample_id", "mecA"]]
    order = pd.DataFrame({"sample_id": [str(i) for i in ids]}).merge(
        gt, on="sample_id", how="left")
    truth = order["mecA"].values.astype(float)
    keep = np.isfinite(truth)
    truth = truth[keep]

    print(f"\ncohort: {int(keep.sum()):,} isolates "
          f"({int(truth.sum()):,} MRSA / {int((truth == 0).sum()):,} MSSA)\n")

    rows = []
    for label, crr, prof, thr in [
            ("VAE (pipeline)",   vae_crr[keep],  vae_profile[keep],  0.50),
            ("MLST 7-locus",     mlst_crr[keep], mlst_profile[keep], 0.50)]:
        m = metrics(truth, crr, thr)
        bm, bt = best_mcc(truth, crr)
        m.update(arm=label, noise_floor=round(mad(prof), 4),
                 best_mcc=bm, best_thr=bt)
        rows.append(m)

    res = pd.DataFrame(rows)[["arm", "auc", "mcc", "best_mcc", "best_thr",
                              "fnr", "ppv", "noise_floor", "n"]]
    print(res.to_string(index=False))

    out = RESULT / "mlst_normalisation_comparison.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
