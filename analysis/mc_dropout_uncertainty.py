#!/usr/bin/env python3
"""
MC Dropout uncertainty quantification for CNVRock copy-ratio estimates.

During standard inference, model.eval() disables dropout → deterministic output.
MC Dropout keeps dropout ACTIVE (model.train()) and samples N forward passes
per sample, giving a distribution over reconstructions and therefore over CRRs.

Per-bin CRR uncertainty = std(x / x̂_i for i=1..N) across MC samples.
Per-gene CRR uncertainty = std(mean_CRR_i) across MC samples.

Applications:
  1. Flag high-uncertainty amplification calls for clinical review
  2. Estimate whether the 4 C. auris 'amplification-only WT' isolates
     have high or low reconstruction uncertainty (confidence in the signal)
  3. Report per-gene calibration: is blaKPC more confidently called than blaSHV?

Run on HPC (needs model checkpoint + GPU):
    /home/lshlt19/miniconda3/envs/cnvrock/bin/python3 \
        analysis/mc_dropout_uncertainty.py --exp 33 --n-mc 50
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
import torch
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "models"))


def get_model_and_store(exp_id: int, hpc_root: Path) -> tuple:
    """Load VAE checkpoint and data store for a given experiment."""
    import yaml, importlib
    cfg_path = hpc_root / "models/experiments" / str(exp_id) / "config.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())

    # Resolve paths
    store_path = (cfg_path.parent / cfg["store_path"]).resolve()
    out_dir    = (cfg_path.parent / cfg["out_dir"]).resolve()
    ckpt_path  = out_dir / "checkpoint.pth"

    # Load architecture — module name starts with digit, use importlib
    sys.path.insert(0, str(hpc_root / "models"))
    arch_mod = importlib.import_module("architectures.06_conv_vae")
    ConvVAE  = arch_mod.ConvVAE

    n_bins_raw = np.load(str(store_path / "counts.npy"),
                         mmap_mode="r").shape[1]
    model = ConvVAE(latent_dim=cfg["latent_dim"], n_bins_raw=n_bins_raw)
    state = torch.load(str(ckpt_path), map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    return model, store_path, cfg


def mc_reconstruct(model: torch.nn.Module,
                   x: torch.Tensor,
                   n_mc: int = 50) -> torch.Tensor:
    """
    N MC-Dropout forward passes.
    Returns tensor of shape (n_mc, batch_size, n_bins).
    """
    model.train()   # KEEP DROPOUT ACTIVE — essential for MC Dropout
    with torch.no_grad():
        # forward returns dict {"recon": (B,L), "z": (mu, logvar)}
        recs = torch.stack([model(x)["recon"] for _ in range(n_mc)], dim=0)
    return recs  # (n_mc, B, L)


def run(exp_id: int, hpc_root: Path, n_mc: int = 50,
        n_samples: int = 500) -> None:
    out_dir = hpc_root / "data/results" / f"mc_dropout_exp{exp_id}"
    out_dir.mkdir(parents=True, exist_ok=True)

    model, store_path, cfg = get_model_and_store(exp_id, hpc_root)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    print(f"Model loaded  |  device={device}  |  n_mc={n_mc}")

    # Load counts + sample IDs
    counts  = np.load(str(store_path / "counts.npy"), mmap_mode="r")
    medians = np.load(str(store_path / "medians.npy"))
    ids     = np.load(str(store_path / "sample_ids.npy"), allow_pickle=True)

    # Normalise (same as training)
    X = counts.astype(np.float32) / (medians[:, None] + 1e-6)

    # Sub-sample for speed
    rng = np.random.default_rng(42)
    idx = rng.choice(len(X), size=min(n_samples, len(X)), replace=False)
    X_sub = torch.tensor(X[idx], device=device)  # (n_sub, L)
    ids_sub = ids[idx]

    print(f"Running MC Dropout on {len(idx)} samples × {n_mc} passes…")

    # Process in batches of 64
    bs = 64
    all_mean = []
    all_std  = []
    for i in range(0, len(X_sub), bs):
        # ConvVAE encoder does .unsqueeze(1) internally — pass (B, L)
        batch = X_sub[i:i+bs]               # (B, L)
        recs  = mc_reconstruct(model, batch, n_mc)  # (n_mc, B, L)
        # CRR per MC sample: x / x̂  (broadcast)
        x_obs = batch                               # (B, L)
        crr_mc = x_obs.unsqueeze(0) / (recs + 1e-6)  # (n_mc, B, L)
        all_mean.append(crr_mc.mean(0).cpu().numpy())  # (B, L)
        all_std.append(crr_mc.std(0).cpu().numpy())    # (B, L) uncertainty
        if i % 200 == 0:
            print(f"  {i}/{len(X_sub)}")

    mean_crr = np.concatenate(all_mean, axis=0)   # (n_sub, L)
    std_crr  = np.concatenate(all_std, axis=0)    # (n_sub, L)

    # Save
    np.save(str(out_dir / "mc_crr_mean.npy"), mean_crr)
    np.save(str(out_dir / "mc_crr_std.npy"), std_crr)
    np.save(str(out_dir / "mc_sample_ids.npy"), ids_sub)
    print(f"Saved MC arrays to {out_dir}/")

    # ── Per-gene CRR uncertainty ──────────────────────────────────────
    # Load gene coords to extract per-gene bin indices
    import yaml
    cfg_path = hpc_root / "models/experiments" / str(exp_id) / "config.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())

    # Load contigs to get bin-to-coordinate mapping
    contigs = np.load(str(store_path / "contigs.npy"), allow_pickle=True)

    gene_uncertainty = []
    gene_coords_path = hpc_root / "assets/kpsc_genome_gene_coords.tsv"
    if not gene_coords_path.exists():
        gene_coords_path = hpc_root / "assets/kpsc_genes_of_interest.tsv"

    # Use blaSHV coordinates as a known example
    # blaSHV locus: NC_016845.1:2549403-2550263
    blashv_start, blashv_end = 2549403, 2550263
    blashv_bins = [i for i, c in enumerate(contigs)
                   if (hasattr(c, '__iter__') and
                       str(c[0]) == "NC_016845.1" and
                       int(c[1]) >= blashv_start and
                       int(c[2]) <= blashv_end + 1000)]

    if blashv_bins:
        gene_mean_crr = mean_crr[:, blashv_bins].mean(axis=1)
        gene_std_crr  = std_crr[:, blashv_bins].mean(axis=1)
        print(f"\n── blaSHV CRR uncertainty ─────────────────────────────────")
        print(f"  Mean CRR: {gene_mean_crr.mean():.3f} ± {gene_std_crr.mean():.4f}")
        print(f"  Samples with CRR > 1.75 (amplified): "
              f"{(gene_mean_crr > 1.75).sum()} / {len(gene_mean_crr)}")

        # Flag high-uncertainty calls
        amp_idx = np.where(gene_mean_crr > 1.75)[0]
        if len(amp_idx):
            print(f"  Amplified sample uncertainties (top 5):")
            top = amp_idx[np.argsort(gene_std_crr[amp_idx])[-5:][::-1]]
            for t in top:
                print(f"    {ids_sub[t]}: CRR={gene_mean_crr[t]:.2f} "
                      f"± {gene_std_crr[t]:.4f}")

    # Overall uncertainty statistics
    genome_uncertainty = std_crr.mean(axis=1)
    print(f"\n── Genome-wide CRR uncertainty ─────────────────────────────")
    print(f"  Mean per-bin std across all samples: {std_crr.mean():.4f}")
    print(f"  99th percentile: {np.percentile(genome_uncertainty, 99):.4f}")
    print(f"  Samples with high uncertainty (>p95): "
          f"{(genome_uncertainty > np.percentile(genome_uncertainty, 95)).sum()}")

    # Save per-sample summary
    summary = pd.DataFrame({
        "sample_id": ids_sub,
        "mean_genome_uncertainty": genome_uncertainty,
        "max_bin_uncertainty":     std_crr.max(axis=1),
    })
    summary.sort_values("mean_genome_uncertainty", ascending=False,
                        inplace=True)
    summary.to_csv(out_dir / "mc_uncertainty_summary.tsv", sep="\t",
                   index=False)
    print(f"\nSaved {out_dir}/mc_uncertainty_summary.tsv")
    print("\nTop 10 highest-uncertainty samples:")
    print(summary.head(10).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", type=int, default=33,
                        help="Experiment ID (default: 33 = KpSC 10K)")
    parser.add_argument("--n-mc", type=int, default=50,
                        help="Number of MC Dropout passes (default: 50)")
    parser.add_argument("--n-samples", type=int, default=500,
                        help="Number of samples to evaluate (default: 500)")
    parser.add_argument("--hpc-root", type=str,
                        default="/home/lshlt19/CNVRock",
                        help="Path to CNVRock repo on HPC")
    args = parser.parse_args()

    hpc_root = Path(args.hpc_root)
    run(args.exp, hpc_root, n_mc=args.n_mc, n_samples=args.n_samples)


if __name__ == "__main__":
    main()
