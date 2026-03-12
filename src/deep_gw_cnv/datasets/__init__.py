"""
Dataset registry — maps config string names to Dataset classes.

Usage
-----
from deep_gw_cnv.datasets import REGISTRY, build_loaders

train_loader, val_loader = build_loaders(cfg)
"""

from torch.utils.data import DataLoader

REGISTRY: dict = {}


def build_loaders(cfg: dict) -> tuple[DataLoader, DataLoader]:
    """Instantiate train and validation DataLoaders from *cfg*."""
    DatasetClass = REGISTRY[cfg["dataset"]]
    train_ds = DatasetClass(cfg, split="train")
    val_ds = DatasetClass(cfg, split="val")
    loader_cfg = cfg.get("dataloader", {})
    train_loader = DataLoader(train_ds, shuffle=True, **loader_cfg)
    val_loader = DataLoader(val_ds, shuffle=False, **loader_cfg)
    return train_loader, val_loader
