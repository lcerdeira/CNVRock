#!/usr/bin/env python3
"""
Merge per-sample Kleborate TSV outputs into a single table.

Usage
-----
python3 data/setup/merge_kleborate_expansion.py \
    --kleborate-dir data/raw/kleborate_expansion/ \
    --assembly-tsv  assets/kpsc_expansion_assembly_urls.tsv \
    --out           assets/kpsc_expansion_kleborate_gt.tsv

Reports how many samples are done vs missing, then writes the merged TSV.
The output can be used by get_kleborate_gt.py for ground-truth extraction.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--kleborate-dir", default="data/raw/kleborate_expansion",
                   help="Directory of per-sample Kleborate TSVs [%(default)s]")
    p.add_argument("--assembly-tsv",  default="assets/kpsc_expansion_assembly_urls.tsv",
                   help="Assembly URL TSV (col 1 = sample_accession) [%(default)s]")
    p.add_argument("--out", default="assets/kpsc_expansion_kleborate_gt.tsv",
                   help="Output merged TSV [%(default)s]")
    p.add_argument("--min-done-frac", type=float, default=0.95,
                   help="Warn if fewer than this fraction of samples are done [%(default)s]")
    return p.parse_args()


def main():
    args = parse_args()
    kleb_dir = Path(args.kleborate_dir)
    asm_tsv  = Path(args.assembly_tsv)

    # Expected samples
    asm_df = pd.read_csv(asm_tsv, sep="\t")
    expected = set(asm_df.iloc[:, 0].astype(str))
    print(f"Expected samples: {len(expected):,}")

    # Available TSVs
    tsv_files = sorted(kleb_dir.glob("*.tsv"))
    done = {f.stem for f in tsv_files}
    print(f"Done:             {len(done):,}  ({100*len(done)/max(len(expected),1):.1f}%)")

    missing = expected - done
    if missing:
        frac_done = len(done) / len(expected)
        if frac_done < args.min_done_frac:
            print(f"WARNING: only {frac_done:.1%} complete — "
                  f"{len(missing):,} samples missing. Re-run array job to fill gaps.",
                  file=sys.stderr)
        else:
            print(f"  ({len(missing):,} samples missing — within acceptable threshold)")

    if not tsv_files:
        print("ERROR: no TSV files found in", kleb_dir, file=sys.stderr)
        sys.exit(1)

    # Read and concatenate
    print(f"\nReading {len(tsv_files):,} TSV files …", flush=True)
    dfs = []
    for i, f in enumerate(tsv_files):
        try:
            df = pd.read_csv(f, sep="\t")
            dfs.append(df)
        except Exception as e:
            print(f"  WARNING: skipping {f.name}: {e}", file=sys.stderr)
        if (i + 1) % 5000 == 0:
            print(f"  {i+1:,} / {len(tsv_files):,} …", flush=True)

    merged = pd.concat(dfs, ignore_index=True)
    print(f"Merged: {len(merged):,} rows, {len(merged.columns)} columns")
    print(f"Columns: {list(merged.columns[:10])} …")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput → {out_path}")

    # Quick AMR gene summary
    amr_cols = [c for c in merged.columns
                if any(g in c for g in ["blaKPC", "blaNDM", "blaCTX", "blaOXA",
                                        "qnrB", "aac", "blaTEM", "blaSHV"])]
    if amr_cols:
        print("\nAMR gene presence (% of samples, Kleborate):")
        for col in sorted(amr_cols)[:15]:
            n_pos = (merged[col].notna() & (merged[col] != "-")).sum()
            print(f"  {col:30s}: {n_pos:6,} / {len(merged):,}  "
                  f"({100*n_pos/max(len(merged),1):.1f}%)")


if __name__ == "__main__":
    main()
