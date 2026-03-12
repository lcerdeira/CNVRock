#!/usr/bin/env python
"""
train.py — run a training experiment.

Usage
-----
    python scripts/train.py configs/experiments/exp_001_baseline.yaml
"""

import argparse
import json
import shutil
from pathlib import Path

import yaml


def load_config(config_path: Path) -> dict:
    """Load a config file, merging with its _base_ if specified."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    base_key = "_base_"
    if base_key in cfg:
        base_path = config_path.parent / cfg.pop(base_key)
        with open(base_path) as f:
            base_cfg = yaml.safe_load(f)
        # Shallow merge: experiment keys override base keys
        base_cfg.update(cfg)
        cfg = base_cfg

    return cfg


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a deep_gw_cnv model.")
    parser.add_argument("config", type=Path, help="Path to experiment YAML config.")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Derive output directory from the config filename
    exp_name = args.config.stem           # e.g. "exp_001_baseline"
    output_dir = Path(cfg.get("output_dir", "experiments")) / exp_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Freeze a copy of the resolved config for reproducibility
    frozen_cfg_path = output_dir / "config.yaml"
    with open(frozen_cfg_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)
    print(f"Config snapshot written to {frozen_cfg_path}")

    # --- Lazily import heavy dependencies so the script fails fast on bad configs
    from deep_gw_cnv.utils import seed_everything

    seed_everything(cfg.get("seed", 42))

    from deep_gw_cnv.datasets import build_loaders
    from deep_gw_cnv.models import REGISTRY

    train_loader, val_loader = build_loaders(cfg)

    ModelClass = REGISTRY[cfg["model"]]
    model = ModelClass(cfg.get("model_cfg", {}))

    # TODO: implement the training loop (optimiser, scheduler, loss, logging)
    print(f"Model: {cfg['model']}  |  dataset: {cfg['dataset']}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
