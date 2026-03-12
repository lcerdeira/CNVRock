#!/usr/bin/env python
"""
predict.py — run inference on new data with a trained checkpoint.

Usage
-----
    python scripts/predict.py configs/experiments/exp_001_baseline.yaml \\
        --input path/to/inputs/ --output path/to/predictions/
"""

import argparse
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description="Run inference with a deep_gw_cnv model.")
    parser.add_argument("config", type=Path, help="Path to experiment YAML config.")
    parser.add_argument("--input", type=Path, required=True, help="Input data path.")
    parser.add_argument("--output", type=Path, required=True, help="Where to write predictions.")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    args.output.mkdir(parents=True, exist_ok=True)

    # TODO: load checkpoint, build dataset from args.input, write predictions
    print(f"Predictions will be written to {args.output}")


if __name__ == "__main__":
    main()
