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
import inspect
import os
import shutil
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
    sample_offset = 0

    # Per-sample medians for correct denormalization (None if normalise=False)
    medians = dataset._medians  # shape (n_samples,) or None

    with torch.no_grad():
        for batch in dl:
            bsz   = batch.shape[0]
            x     = batch.to(device)
            mu, _ = model.enc(x)
            recon = model.dec(mu)                       # in log2(count/median+1) space
            # Denormalize: invert log2(count/median+1) → count/median → count
            recon_np = recon.cpu().numpy()
            if medians is not None:
                batch_medians = medians[sample_offset:sample_offset + bsz, np.newaxis]
                recon_raw = (np.power(2, recon_np) - 1) * batch_medians
            else:
                recon_raw = np.power(2, recon_np) - 1
            all_mu.append(mu.cpu().numpy())
            all_recon.append(recon_raw)
            sample_offset += bsz

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
    ds = ReadCountDataset(store_path, normalise=cfg["normalise"])
    n_bins_raw = ds.counts.shape[1]
    if "n_bins_raw" in inspect.signature(ConvVAE.__init__).parameters:
        model = ConvVAE(latent_dim=cfg["latent_dim"], n_bins_raw=n_bins_raw).to(device)
    else:
        model = ConvVAE(latent_dim=cfg["latent_dim"]).to(device)
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    print(f"Loaded checkpoint: {checkpoint_path}", flush=True)

    # ── inference ───────────────────────────────────────────────────────────
    run_inference(model, ds, device, out_dir, batch_size=cfg["batch_size"])

    # ── HMM segmentation ────────────────────────────────────────────────────
    print("Fitting HMM segments...", flush=True)
    run_hmm_all_samples(store_path, out_dir, cfg)

    # ── CNV calling (chromosomal) ────────────────────────────────────────────
    print("Calling gene CNVs...", flush=True)
    run_cnv_calls(store_path, out_dir, cfg)

    # ── CNV calling (plasmid genes) — optional ───────────────────────────────
    if cfg.get("plasmid_store_path") and cfg.get("plasmid_gene_coords_path"):
        print("Calling plasmid gene CNVs...", flush=True)
        run_plasmid_cnv_calls = importlib.import_module(
            f"cnv.{cfg.get('plasmid_cnv', '07_plasmid_cnv_caller')}"
        ).run_plasmid_cnv_calls
        plasmid_cfg = dict(cfg)
        plasmid_cfg["store_path"]               = store_path   # already resolved above
        plasmid_cfg["plasmid_store_path"]       = resolve(cfg["plasmid_store_path"])
        plasmid_cfg["plasmid_gene_coords_path"] = resolve(cfg["plasmid_gene_coords_path"])
        run_plasmid_cnv_calls(out_dir, plasmid_cfg)
    else:
        print("Skipping plasmid CNV calling (plasmid_store_path not set).", flush=True)

    # ── Evaluation (optional) ────────────────────────────────────────────────
    if cfg.get("pf9_gt_path"):
        cfg_resolved = dict(cfg)
        cfg_resolved["pf9_gt_path"] = resolve(cfg["pf9_gt_path"])
        if cfg.get("pf9_meta_path"):
            cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
        run_evaluation(out_dir, cfg_resolved)
    elif cfg.get("kpsc_gt_path"):
        cfg_resolved = dict(cfg)
        cfg_resolved["kpsc_gt_path"] = resolve(cfg["kpsc_gt_path"])
        if cfg.get("kpsc_meta_path"):
            cfg_resolved["kpsc_meta_path"] = resolve(cfg["kpsc_meta_path"])
        run_evaluation(out_dir, cfg_resolved)
    else:
        print("Skipping evaluation (no gt_path set in config).", flush=True)

    # ── Storage report ───────────────────────────────────────────────────────
    out_bytes = sum(
        os.path.getsize(os.path.join(dp, f))
        for dp, _, files in os.walk(out_dir)
        for f in files
    )
    disk       = shutil.disk_usage(out_dir)
    out_gb     = out_bytes      / 1024 ** 3
    free_gb    = disk.free      / 1024 ** 3
    total_gb   = disk.total     / 1024 ** 3
    used_pct   = 100 * disk.used / disk.total
    print(
        f"\nStorage — output dir: {out_gb:.2f} GB  |  "
        f"disk free: {free_gb:.1f} / {total_gb:.1f} GB  ({used_pct:.1f}% used)",
        flush=True,
    )


if __name__ == "__main__":
    main()
