"""
Static screenshot generator for KpSC diagnostics (no Streamlit required).

Generates key diagnostic plots from experiment outputs and saves them as PNGs.

Usage (from repo root):
    python3 diagnostics/generate_screenshots.py --exp 29
    python3 diagnostics/generate_screenshots.py --exp 29 --sample ERR123456
"""

import argparse
import importlib
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../models"))

_REPO_ROOT   = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
_EXPS_DIR    = os.path.join(_REPO_ROOT, "models", "experiments")
_OUT_DIR     = os.path.join(_REPO_ROOT, "assets", "screenshots")

CN_COLORS = ["#313695", "#74add1", "#888888", "#fdae61", "#d73027", "#a50026"]


def _resolve(cfg, key, config_dir):
    v = cfg.get(key)
    if v and not os.path.isabs(v):
        return os.path.normpath(os.path.join(config_dir, v))
    return v


def load_exp(exp_id):
    config_dir  = os.path.join(_EXPS_DIR, str(exp_id))
    config_path = os.path.join(config_dir, "config.yaml")
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    out_dir    = _resolve(cfg, "out_dir",    config_dir)
    store_path = _resolve(cfg, "store_path", config_dir)
    return cfg, config_dir, out_dir, store_path


def load_data(out_dir, store_path):
    latents         = np.load(os.path.join(out_dir,    "latents.npy"))
    reconstructions = np.load(os.path.join(out_dir,    "reconstructions.npy"))
    sample_ids      = np.load(os.path.join(out_dir,    "sample_ids.npy"),    allow_pickle=True)
    counts          = np.load(os.path.join(store_path, "counts.npy"))
    contigs_raw     = np.load(os.path.join(store_path, "contigs.npy"),       allow_pickle=True)
    store_ids       = np.load(os.path.join(store_path, "sample_ids.npy"),    allow_pickle=True)
    segs_path       = os.path.join(out_dir, "segments.parquet")
    segments        = pd.read_parquet(segs_path) if os.path.exists(segs_path) else None
    gene_calls      = pd.read_csv(os.path.join(out_dir, "gene_calls.tsv"), sep="\t", index_col=0)
    plasmid_path    = os.path.join(out_dir, "plasmid_gene_calls.tsv")
    plasmid_calls   = pd.read_csv(plasmid_path, sep="\t", index_col=0) if os.path.exists(plasmid_path) else None
    contigs_df      = pd.DataFrame(contigs_raw)
    counts_df       = pd.DataFrame(counts, index=store_ids)
    return {
        "latents":       pd.DataFrame(latents,         index=sample_ids),
        "reconstructions": pd.DataFrame(reconstructions, index=sample_ids),
        "counts":        counts_df,
        "contigs":       contigs_df,
        "segments":      segments,
        "gene_calls":    gene_calls,
        "plasmid_calls": plasmid_calls,
        "sample_ids":    sample_ids,
    }


# ---------------------------------------------------------------------------
# Plot 1: PCA of latent space
# ---------------------------------------------------------------------------

def plot_pca_latents(data, out_path, highlight=None):
    pca    = PCA(n_components=2)
    coords = pca.fit_transform(data["latents"].values)
    var    = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(coords[:, 0], coords[:, 1], s=6, alpha=0.4, c="#888888")

    if highlight and highlight in data["latents"].index:
        idx = list(data["latents"].index).index(highlight)
        ax.scatter(coords[idx, 0], coords[idx, 1], s=80, c="red", zorder=5, label=highlight)
        ax.legend(fontsize=8)

    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title("Latent space PCA — KpSC exp 29")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Read-depth track for one sample (chromosome 1)
# ---------------------------------------------------------------------------

def plot_depth_track(data, sample_id, chrom, out_path):
    contigs = data["contigs"]
    mask    = contigs["chrom"] == chrom
    if not mask.any():
        print(f"  Chrom {chrom} not found — skipping depth track")
        return

    sample_counts = data["counts"].loc[sample_id].values
    sample_recon  = data["reconstructions"].loc[sample_id].values

    chrom_idx = contigs.index[mask].tolist()
    starts    = contigs.loc[mask, "start"].values.astype(float)
    raw_c     = sample_counts[chrom_idx]
    recon_c   = sample_recon[chrom_idx]
    ratio     = raw_c.astype(float) / (recon_c.astype(float) + 1e-6)

    # Insert NaN at large gaps
    gaps = np.where(np.diff(starts) > 1000)[0] + 1
    xs   = np.insert(starts, gaps, np.nan)
    ys   = np.insert(ratio,  gaps, np.nan)
    yr   = np.insert(raw_c.astype(float),   gaps, np.nan)
    yrec = np.insert(recon_c.astype(float), gaps, np.nan)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
    ax1.plot(xs, ys, lw=0.7, c="#777777")
    ax1.set_ylim(-0.1, 5.1)
    ax1.set_ylabel("copy ratio")
    ax1.set_title(f"{sample_id} — {chrom}")
    ax1.axhline(1.0, lw=0.5, ls="--", c="black", alpha=0.5)

    if data["segments"] is not None:
        segs = data["segments"]
        segs = segs[(segs["sample_id"] == sample_id) & (segs["chrom"] == chrom)]
        for _, seg in segs.iterrows():
            cn = int(min(seg["cn"], len(CN_COLORS) - 1))
            ax1.plot([seg["x0"], seg["x1"]], [seg["cn"], seg["cn"]],
                     c=CN_COLORS[cn], lw=2.5)

    ax2.plot(xs, yr,   lw=1.0, c="green", label="input", alpha=0.8)
    ax2.plot(xs, yrec, lw=1.0, c="red",   label="recon", alpha=0.8)
    ax2.set_ylabel("read depth")
    ax2.set_xlabel("genomic position")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 3: PCN distributions per plasmid gene
# ---------------------------------------------------------------------------

def plot_pcn_distributions(data, out_path):
    calls = data["plasmid_calls"]
    if calls is None:
        print("  No plasmid calls — skipping PCN plot")
        return

    pcn_cols = [c for c in calls.columns if c.startswith("pcn_")]
    if not pcn_cols:
        return

    n = len(pcn_cols)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.array(axes).flatten()

    for i, col in enumerate(pcn_cols):
        gene = col.replace("pcn_", "")
        ax   = axes[i]
        vals = calls[col].dropna()
        n_present = (calls.get(gene, pd.Series(dtype=int)) >= 1).sum()
        ax.hist(vals[vals < 5], bins=50, color="#4a90d9", alpha=0.75, edgecolor="none")
        ax.set_title(f"{gene}\n({n_present} present)", fontsize=9)
        ax.set_xlabel("PCN", fontsize=8)
        ax.set_ylabel("n samples", fontsize=8)
        ax.tick_params(labelsize=7)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Plasmid copy number distributions — KpSC exp 29", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 4: Latent vector for one sample
# ---------------------------------------------------------------------------

def plot_latent_bar(data, sample_id, out_path):
    vals   = data["latents"].loc[sample_id].values.astype(float)
    labels = [f"z{i+1}" for i in range(len(vals))]
    abs_max = np.abs(vals).max() or 1.0
    norm    = plt.Normalize(vmin=0, vmax=abs_max)
    colors  = plt.cm.Reds(norm(np.abs(vals)))

    fig, ax = plt.subplots(figsize=(8, 2))
    bars = ax.bar(labels, vals, color=colors, edgecolor="none")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylim(-abs_max * 1.3, abs_max * 1.3)
    ax.set_ylabel("latent value", fontsize=9)
    ax.set_title(f"Latent vector — {sample_id}", fontsize=9)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp",    default="29",   help="Experiment folder name")
    ap.add_argument("--sample", default=None,   help="Sample ID to highlight/detail")
    ap.add_argument("--chrom",  default="NC_016845.1", help="Chromosome for depth track")
    args = ap.parse_args()

    cfg, config_dir, out_dir, store_path = load_exp(args.exp)
    print(f"Loading exp {args.exp}: {out_dir}")
    data = load_data(out_dir, store_path)
    print(f"  {len(data['sample_ids'])} samples, {data['latents'].shape[1]} latent dims")

    os.makedirs(_OUT_DIR, exist_ok=True)
    prefix = os.path.join(_OUT_DIR, f"exp{args.exp}")

    sample_id = args.sample
    if sample_id is None:
        # Pick a KPC-positive sample for the detail plots
        if data["plasmid_calls"] is not None and "blaKPC-2" in data["plasmid_calls"].columns:
            pos = data["plasmid_calls"][data["plasmid_calls"]["blaKPC-2"] >= 1].index.tolist()
            sample_id = pos[0] if pos else data["sample_ids"][0]
        else:
            sample_id = data["sample_ids"][0]
    print(f"  Detail sample: {sample_id}")

    plot_pca_latents(data,   f"{prefix}_pca_latents.png", highlight=sample_id)
    plot_pcn_distributions(data, f"{prefix}_pcn_distributions.png")
    plot_latent_bar(data, sample_id, f"{prefix}_latent_bar_{sample_id}.png")
    plot_depth_track(data, sample_id, args.chrom, f"{prefix}_depth_{sample_id}.png")

    print(f"\nAll screenshots saved to {_OUT_DIR}/")


if __name__ == "__main__":
    main()
