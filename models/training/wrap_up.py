"""
Wrap-up script: runs after a model has finished training.

Steps
-----
1. Inference — encode every sample with mu (deterministic), decode, denormalise.
   Writes: latents.npy, reconstructions.npy, sample_ids.npy
2. HMM segmentation — see hmm/  (versioned)
   Writes: segments.parquet  [sample_id, chrom, x0, x1, cn, confidence]
3. CNV calling — see cnv/  (versioned)
   Writes: gene_calls.tsv
4. Evaluation — see evaluation/  (versioned, optional)
   Writes: evaluation.txt

Usage (standalone — re-run wrap-up on an existing checkpoint):
    python -m training.wrap_up path/to/config.yaml [path/to/checkpoint.pth]

Imported by train.py:
    from training.wrap_up import run_inference
"""

import argparse
import os
import sys

import numpy as np
import torch
from torch.utils.data import DataLoader


# ---------------------------------------------------------------------------
# 1. Inference
# ---------------------------------------------------------------------------

def run_inference(model, dataset, device, out_dir, batch_size=128):
    """Encode every sample with mu (deterministic) and save outputs.

    Outputs written to out_dir:
        latents.npy          — (n_samples, latent_dim)  mu vectors
        reconstructions.npy  — (n_samples, n_bins)      raw count space (denormalised)
        sample_ids.npy       — (n_samples,)             sample ID strings
    """
    os.makedirs(out_dir, exist_ok=True)
    model.eval()

    dl = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    all_mu    = []
    all_recon = []

    with torch.no_grad():
        for batch in dl:
            x     = batch.to(device)
            mu, _ = model.enc(x)
            recon = model.dec(mu)                       # deterministic: use mu
            recon_denorm = torch.pow(2, recon) - 1      # inverse of log2(count+1)
            all_mu.append(mu.cpu().numpy())
            all_recon.append(recon_denorm.cpu().numpy())

    latents = np.concatenate(all_mu,    axis=0)         # (n_samples, latent_dim)
    recons  = np.concatenate(all_recon, axis=0)         # (n_samples, n_bins)

    np.save(os.path.join(out_dir, "latents.npy"),         latents)
    np.save(os.path.join(out_dir, "reconstructions.npy"), recons)
    np.save(os.path.join(out_dir, "sample_ids.npy"),      np.array(dataset.sample_ids))

    print(f"Saved latents         {latents.shape} → {out_dir}/latents.npy",       flush=True)
    print(f"Saved reconstructions {recons.shape}  → {out_dir}/reconstructions.npy", flush=True)
    n = len(dataset.sample_ids)
    print(f"Saved sample_ids      ({n},) → {out_dir}/sample_ids.npy",             flush=True)


# ---------------------------------------------------------------------------
# Entry point (standalone re-run of wrap-up on an existing checkpoint)
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Re-run inference, HMM fitting, CNV calling, and evaluation on an existing checkpoint."
    )
    parser.add_argument("config",     help="Path to experiment config.yaml")
    parser.add_argument("checkpoint", nargs="?",
                        help="Path to checkpoint.pth (default: out_dir/checkpoint.pth)")
    args = parser.parse_args()

    # Add models/ to path so package imports work when invoked as a module
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import importlib  # noqa: PLC0415
    import yaml       # noqa: PLC0415

    from training.dataset import ReadCountDataset  # noqa: PLC0415

    config_dir = os.path.dirname(os.path.abspath(args.config))
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    # Load versioned components named in config
    ConvVAE             = importlib.import_module(f"architectures.{cfg['architecture']}").ConvVAE
    run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples
    run_cnv_calls       = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
    run_evaluation      = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

    def resolve(path):
        return path if os.path.isabs(path) else os.path.join(config_dir, path)

    store_path      = resolve(cfg["store_path"])
    out_dir         = resolve(cfg["out_dir"])
    checkpoint_path = args.checkpoint or os.path.join(out_dir, "checkpoint.pth")

    # ── device ──────────────────────────────────────────────────────────────
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    print(f"Device: {device}", flush=True)

    # ── load model ──────────────────────────────────────────────────────────
    model = ConvVAE(latent_dim=cfg["latent_dim"]).to(device)
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    print(f"Loaded checkpoint: {checkpoint_path}", flush=True)

    # ── inference ───────────────────────────────────────────────────────────
    ds = ReadCountDataset(store_path, normalise=cfg["normalise"])
    run_inference(model, ds, device, out_dir, batch_size=cfg["batch_size"])

    # ── HMM segmentation ────────────────────────────────────────────────────
    print("Fitting HMM segments...", flush=True)
    run_hmm_all_samples(store_path, out_dir, cfg)

    # ── CNV calling ─────────────────────────────────────────────────────────
    print("Calling gene CNVs...", flush=True)
    run_cnv_calls(store_path, out_dir, cfg)

    # ── Evaluation (optional — requires pf9_gt_path in config) ──────────────
    if cfg.get("pf9_gt_path"):
        cfg_resolved = dict(cfg)
        cfg_resolved["pf9_gt_path"] = resolve(cfg["pf9_gt_path"])
        if cfg.get("pf9_meta_path"):
            cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
        run_evaluation(out_dir, cfg_resolved)
    else:
        print("Skipping evaluation (pf9_gt_path not set in config).", flush=True)


if __name__ == "__main__":
    main()
