"""
Model registry — maps config string names to nn.Module classes.

Usage
-----
from deep_gw_cnv.models import REGISTRY

ModelClass = REGISTRY[cfg["model"]]
model = ModelClass(cfg["model_cfg"])
"""

REGISTRY: dict = {}
