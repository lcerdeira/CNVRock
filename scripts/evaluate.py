#!/usr/bin/env python
"""
evaluate.py — evaluate a trained model checkpoint.

Usage
-----
    python scripts/evaluate.py configs/experiments/exp_001_baseline.yaml [--split test]
"""

import argparse
import json
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate a deep_gw_cnv checkpoint.")
    parser.add_argument("config", type=Path, help="Path to experiment YAML config.")
    parser.add_argument(
        "--split",
        default="val",
        choices=["val", "test"],
        help="Data split to evaluate on (default: val).",
    )
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    exp_name = args.config.stem
    output_dir = Path(cfg.get("output_dir", "experiments")) / exp_name
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # TODO: load checkpoint, run inference, compute metrics
    metrics: dict = {}
    metrics_path = results_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"Metrics written to {metrics_path}")


if __name__ == "__main__":
    main()
