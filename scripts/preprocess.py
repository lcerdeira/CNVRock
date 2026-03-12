#!/usr/bin/env python
"""
preprocess.py — convert raw data into the processed format consumed by Dataset classes.

Usage
-----
    python scripts/preprocess.py --raw data/raw/ --out data/processed/
"""

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess raw data.")
    parser.add_argument("--raw", type=Path, default=Path("data/raw/"), help="Raw data directory.")
    parser.add_argument("--out", type=Path, default=Path("data/processed/"), help="Output directory.")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    # TODO: implement preprocessing pipeline
    print(f"Raw data: {args.raw}")
    print(f"Output:   {args.out}")


if __name__ == "__main__":
    main()
