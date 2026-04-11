"""
Central training entry point.

Usage:
    python train.py path/to/config.yaml

All paths in the config are resolved relative to the config file's directory,
so experiments stay self-contained regardless of where you invoke this script.
"""

import argparse
import os
import sys

import yaml
import torch
from torch.utils.data import DataLoader

# Allow running this file from any working directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architectures import ConvVAE
from training.dataset import ReadCountDataset
from training.trainer import train_vae
from training.wrap_up import run_hmm_all_samples, run_inference


def get_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main():
    parser = argparse.ArgumentParser(description="Train a model from a config file.")
    parser.add_argument("config", help="Path to experiment config.yaml")
    args = parser.parse_args()

    config_dir = os.path.dirname(os.path.abspath(args.config))
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    def resolve(path):
        return path if os.path.isabs(path) else os.path.join(config_dir, path)

    store_path = resolve(cfg["store_path"])
    out_dir    = resolve(cfg["out_dir"])

    lsf_jobid = os.getenv("LSB_JOBID")
    if lsf_jobid:
        print(f"LSF job: {lsf_jobid}\n{'='*50}", flush=True)

    device = get_device()
    print(f"Device: {device}", flush=True)

    ds = ReadCountDataset(store_path, normalise=cfg.get("normalise", True))

    num_workers = 8 if device.type == "cuda" else 0
    dl = DataLoader(ds, batch_size=cfg["batch_size"], shuffle=True, num_workers=num_workers)

    model     = ConvVAE(latent_dim=cfg["latent_dim"]).to(device)
    optimiser = torch.optim.Adam(
        model.parameters(),
        lr=cfg.get("lr", 1e-3),
        weight_decay=cfg.get("weight_decay", 1e-5),
    )

    print(f"samples={len(ds)} | latent_dim={cfg['latent_dim']}", flush=True)

    os.makedirs(out_dir, exist_ok=True)
    checkpoint_path = os.path.join(out_dir, "checkpoint.pth")
    log_path        = os.path.join(out_dir, "training_log.json")

    train_vae(
        model, dl, optimiser,
        epochs          = cfg["epochs"],
        max_beta        = cfg.get("max_beta", 1.0),
        warmup_epochs   = cfg.get("warmup_epochs", 400),
        patience        = cfg.get("patience", 50),
        device          = device,
        model_save_path = checkpoint_path,
        log_path        = log_path,
    )

    if os.path.exists(checkpoint_path):
        model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
        print("Loaded best checkpoint for inference.", flush=True)

    run_inference(model, ds, device, out_dir, batch_size=cfg["batch_size"])

    print("Fitting HMM segments...", flush=True)
    run_hmm_all_samples(store_path, out_dir)


if __name__ == "__main__":
    main()
