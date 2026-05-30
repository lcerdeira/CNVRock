#!/usr/bin/env python3
"""
Intrinsic dimensionality of the VAE latent space (exp33, 10K KpSC).

Methods:
  1. Two-NN estimator (Facco et al., Nature Communications 2017) — the
     maximum-likelihood estimator of intrinsic dimension from the ratio
     r2/r1 of distances to the 2nd vs 1st nearest neighbour.
  2. PCA cumulative explained variance — standard linear lower bound.
  3. Neighbourhood structure visualisation (per-gene CN colouring).

Question: are 10 latent dimensions justified?
If the intrinsic dimension is 3–4, the remaining dimensions encode noise
that the HMM segmenter subsequently tries to separate as signal.

Run locally (latents.npy available):
    python3 analysis/latent_intrinsic_dim.py
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr

REPO     = Path(__file__).resolve().parents[1]
LATENTS  = REPO / "data/results/33_kpsc_expansion_10k/latents.npy"
CALLS    = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
PLASMID  = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
OUT      = REPO / "data/results/latent_intrinsic_dim"


# ── Two-NN intrinsic dimension estimator ─────────────────────────────────────

def twonn_dim(X: np.ndarray, n_subsample: int = 2000) -> tuple[float, float]:
    """
    Facco et al. 2017 two-NN maximum-likelihood estimator.
    Returns (d_hat, 95% CI half-width).
    """
    rng = np.random.default_rng(42)
    idx = rng.choice(len(X), size=min(n_subsample, len(X)), replace=False)
    X_sub = X[idx]
    nbrs = NearestNeighbors(n_neighbors=2, metric="euclidean").fit(X_sub)
    dists, _ = nbrs.kneighbors(X_sub)
    r1 = dists[:, 0]          # distance to 1st NN
    r2 = dists[:, 1]          # distance to 2nd NN
    # exclude points where r1 == 0 (duplicates)
    mask = r1 > 0
    mu = r2[mask] / r1[mask]  # must be > 1 by definition
    mu = mu[mu > 1]
    n = len(mu)
    # MLE for Pareto: d = n / sum(log(mu))
    d_hat = n / np.sum(np.log(mu))
    # Bootstrap 95% CI
    bs_d = []
    for _ in range(200):
        idx_b = rng.integers(0, n, size=n)
        bs_d.append(n / np.sum(np.log(mu[idx_b])))
    ci = 1.96 * np.std(bs_d)
    return d_hat, ci


# ── PCA variance analysis ─────────────────────────────────────────────────────

def pca_analysis(X: np.ndarray) -> np.ndarray:
    pca = PCA(n_components=X.shape[1])
    pca.fit(X)
    return pca.explained_variance_ratio_


# ── Latent–gene correlation ───────────────────────────────────────────────────

def latent_gene_correlation(Z: np.ndarray, calls_path: Path,
                             plasmid_path: Path) -> dict:
    """
    Pearson r between each latent dimension and per-gene CRR/PCN values.
    Reveals which dimensions encode resistance signal vs technical variation.
    """
    import pandas as pd
    chrom = pd.read_csv(calls_path, sep="\t")
    plasm = pd.read_csv(plasmid_path, sep="\t")
    merged = chrom.merge(plasm, on="sample_id", how="inner")

    # Use only samples present in latents (may differ from full eval set)
    ids = np.load(str(calls_path).replace("gene_calls.tsv", "sample_ids.npy"),
                  allow_pickle=True)
    id_to_idx = {sid: i for i, sid in enumerate(ids)}
    merged["_idx"] = merged["sample_id"].map(id_to_idx)
    merged = merged.dropna(subset=["_idx"])
    merged["_idx"] = merged["_idx"].astype(int)

    Z_sub = Z[merged["_idx"].values]

    genes = {
        "crr_blaSHV":   chrom.columns,
        "pcn_blaKPC-2": plasm.columns,
        "pcn_blaNDM-1": plasm.columns,
        "pcn_blaCTX-M-15": plasm.columns,
    }
    results = {}
    for gene_col in ["crr_blaSHV", "pcn_blaKPC-2", "pcn_blaNDM-1",
                      "pcn_blaCTX-M-15"]:
        all_cols = list(merged.columns)
        if gene_col not in all_cols:
            continue
        vals = merged[gene_col].fillna(0).astype(float).values
        corrs = [pearsonr(Z_sub[:, d], vals)[0] for d in range(Z_sub.shape[1])]
        results[gene_col] = corrs
    return results


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    Z = np.load(str(LATENTS))
    print(f"Latent space: {Z.shape}  (n_samples × n_dims)")

    # 1. Two-NN intrinsic dimension
    print("\n── Two-NN intrinsic dimension ──────────────────────────────")
    d_hat, ci = twonn_dim(Z)
    print(f"  d̂ = {d_hat:.2f} ± {ci:.2f}  (95% CI)")
    if d_hat < Z.shape[1] * 0.6:
        print(f"  → The 10-dim latent space is OVER-DIMENSIONED.")
        print(f"    True intrinsic dimension ≈ {d_hat:.0f}. "
              f"~{Z.shape[1] - round(d_hat):.0f} dimensions encode noise/redundancy.")
    else:
        print(f"  → 10 dimensions are approximately justified (d̂ ≈ {d_hat:.1f}).")

    # 2. PCA
    print("\n── PCA cumulative explained variance ───────────────────────")
    evr = pca_analysis(Z)
    cumvar = np.cumsum(evr)
    for k in [1, 2, 3, 4, 5, 7, 10]:
        print(f"  PC1–{k}: {100*cumvar[k-1]:.1f}% variance explained")
    n_90 = int(np.searchsorted(cumvar, 0.90)) + 1
    print(f"  → {n_90} PCs explain ≥90% of latent variance.")

    # 3. Latent–gene correlation
    print("\n── Latent dimension → gene signal correlation ───────────────")
    if CALLS.exists() and PLASMID.exists():
        corr_map = latent_gene_correlation(Z, CALLS, PLASMID)
        for gene, corrs in corr_map.items():
            best_dim = int(np.argmax(np.abs(corrs)))
            best_r   = corrs[best_dim]
            print(f"  {gene:20s}: best dim={best_dim+1}  r={best_r:+.3f}")

    # 4. Figures
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # PCA scree
    axes[0].bar(range(1, 11), 100 * evr, color="#2176ae", alpha=0.85)
    axes[0].plot(range(1, 11), 100 * cumvar, "o-", color="#e63946", lw=1.5)
    axes[0].axhline(90, color="grey", ls="--", lw=0.8)
    axes[0].set_xlabel("Principal component", fontsize=9)
    axes[0].set_ylabel("% variance", fontsize=9)
    axes[0].set_title(
        f"Latent space PCA — intrinsic dim d̂ = {d_hat:.1f} ± {ci:.1f}",
        fontsize=9)

    # Correlation heatmap
    if corr_map:
        import matplotlib.cm as cm
        gene_labels = list(corr_map.keys())
        mat = np.array([corr_map[g] for g in gene_labels])
        im = axes[1].imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1,
                            aspect="auto")
        axes[1].set_xticks(range(10))
        axes[1].set_xticklabels([f"z{i+1}" for i in range(10)], fontsize=7)
        axes[1].set_yticks(range(len(gene_labels)))
        axes[1].set_yticklabels(gene_labels, fontsize=8)
        axes[1].set_title("Pearson r: latent dim × gene CRR/PCN", fontsize=9)
        plt.colorbar(im, ax=axes[1], shrink=0.8)

    plt.tight_layout()
    out_fig = OUT / "latent_intrinsic_dim.png"
    plt.savefig(out_fig, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved: {out_fig}")

    # Summary
    print(f"\n{'='*55}")
    print(f"SUMMARY")
    print(f"  Declared latent dims:  10")
    print(f"  Two-NN intrinsic dim:  {d_hat:.2f} ± {ci:.2f}")
    print(f"  PCs for 90% variance:  {n_90}")
    print(f"  Implication: {'over-parameterised — consider d=' + str(round(d_hat)) if d_hat < 7 else 'well-calibrated'}")


if __name__ == "__main__":
    main()
