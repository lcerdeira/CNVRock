#!/usr/bin/env python3
"""
Build a CNVRock NPY store from KpSC GATK CollectReadCounts TSVs.

Replaces the legacy P. falciparum-specific readcounts_to_npy.py with:
  - CLI args (no hardcoded paths)
  - Manifest-driven sample selection (for the 5k/10k/20k/80k tiers)
  - Configurable contig filter — splits the combined chromosome + plasmid
    count files produced by hpc/aspera_subset_pipeline.sh into either a
    chromosome-only store or a plasmid-only store.

Inputs
------
    --counts-dir DIR        Folder with {ACC}.counts.tsv files
    --manifest FILE         TSV with `accession` column (subset to include)
    --out-dir DIR           Output directory for NPY store
    --keep-contigs CSV      Comma-separated contig names to keep
                             (e.g. NC_016845.1 for chrom store;
                              MK552109.1,MZ606384.2,... for plasmid store)
    --workers N             Parallel readers (default 16)

Outputs (compatible with CNVRock loader)
----------------------------------------
    contigs.npy     structured array [(chrom, start, end)] per bin
    counts.npy      uint32 (n_samples, n_bins)
    sample_ids.npy  object array of accessions, same row order as counts

Example — build the 5k chromosome store:
    python3 data/setup/readcounts_to_npy_kpsc.py \\
        --counts-dir data/raw/readcounts_subset \\
        --manifest   assets/kpsc_expansion_subset_5k.tsv \\
        --out-dir    data/inputs/KpSC-expansion-5k-1000bp-npy \\
        --keep-contigs NC_016845.1

Example — build the 5k plasmid store:
    python3 data/setup/readcounts_to_npy_kpsc.py \\
        --counts-dir data/raw/readcounts_subset \\
        --manifest   assets/kpsc_expansion_subset_5k.tsv \\
        --out-dir    data/inputs/KpSC-expansion-5k-plasmid-1000bp-npy \\
        --keep-contigs MK552109.1,MZ606384.2,NZ_CP031850.1,NZ_KX236178.1,NC_016980.1,CP006662.2,JN626286.1,CP034201.2,CP113224.1,CP030319.1,CP138680.1
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

REQUIRED_COLUMNS = {"CONTIG", "START", "END", "COUNT"}


def read_counts_tsv(path: Path, keep_contigs: set[str]) -> pd.DataFrame:
    """Parse one GATK CollectReadCounts TSV. Returns rows for keep_contigs only."""
    if not path.exists():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")

    with open(path) as f:
        skip = sum(1 for line in f if line.startswith("@"))

    df = pd.read_csv(path, sep="\t", skiprows=skip)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    df = df[df["CONTIG"].isin(keep_contigs)].copy()
    if df.empty:
        raise ValueError(
            f"{path}: no rows match keep_contigs (saw {df['CONTIG'].unique().tolist()})"
        )

    df["START"] = df["START"].astype(np.uint32)
    df["END"]   = df["END"].astype(np.uint32)
    df["COUNT"] = df["COUNT"].astype(np.uint32)
    if df.isnull().any().any():
        raise ValueError(f"{path}: NaNs after parsing — file may be truncated")
    return df.rename(columns={"CONTIG": "CHROM"})[["CHROM", "START", "END", "COUNT"]]


def validate_all(counts_dir: Path, sample_ids: list[str],
                 keep_contigs: set[str], workers: int) -> int:
    """Read every file, check shape consistency. Returns the common bin count."""
    print(f"Validating {len(sample_ids):,} files …")
    errors: list[str] = []
    bin_counts: dict[str, int] = {}

    def one(sid: str):
        path = counts_dir / f"{sid}.counts.tsv"
        try:
            df = read_counts_tsv(path, keep_contigs)
            return sid, len(df), None
        except Exception as exc:                    # noqa: BLE001
            return sid, None, str(exc)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for sid, nb, err in tqdm(pool.map(one, sample_ids), total=len(sample_ids)):
            if err:
                errors.append(f"  {sid}: {err}")
            else:
                bin_counts[sid] = nb

    if errors:
        raise RuntimeError(
            f"{len(errors)} files failed validation:\n"
            + "\n".join(errors[:20])
            + (f"\n  … +{len(errors) - 20} more" if len(errors) > 20 else "")
        )

    shapes = set(bin_counts.values())
    if len(shapes) > 1:
        freq = Counter(bin_counts.values())
        modal = freq.most_common(1)[0][0]
        outliers = [f"  {sid}: {n}" for sid, n in bin_counts.items() if n != modal]
        raise RuntimeError(
            f"Inconsistent bin counts: {dict(freq)}.\n"
            f"Outliers (vs modal={modal}):\n" + "\n".join(outliers[:20])
        )

    n_bins = shapes.pop()
    print(f"Pre-flight OK: {len(sample_ids):,} samples × {n_bins:,} bins")
    return n_bins


def build_store(counts_dir: Path, sample_ids: list[str], keep_contigs: set[str],
                out_dir: Path, workers: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    n_bins = validate_all(counts_dir, sample_ids, keep_contigs, workers)
    n_samples = len(sample_ids)

    # Anchor contigs structured array from the first sample (already validated)
    first = read_counts_tsv(counts_dir / f"{sample_ids[0]}.counts.tsv", keep_contigs)
    contigs = np.empty(
        n_bins,
        dtype=[("chrom", object), ("start", np.uint32), ("end", np.uint32)],
    )
    contigs["chrom"] = first["CHROM"].values
    contigs["start"] = first["START"].values
    contigs["end"]   = first["END"].values
    np.save(out_dir / "contigs.npy", contigs)
    np.save(out_dir / "sample_ids.npy",
            np.array(sample_ids, dtype=object))
    print(f"Wrote contigs.npy ({n_bins} bins) + sample_ids.npy ({n_samples})")

    counts = np.zeros((n_samples, n_bins), dtype=np.uint32)

    def one(args):
        idx, sid = args
        df = read_counts_tsv(counts_dir / f"{sid}.counts.tsv", keep_contigs)
        if len(df) != n_bins:
            raise RuntimeError(
                f"bin mismatch for {sid}: expected {n_bins}, got {len(df)}"
            )
        return idx, df["COUNT"].values

    jobs = list(enumerate(sample_ids))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for idx, vals in tqdm(pool.map(one, jobs), total=n_samples):
            counts[idx, :] = vals

    np.save(out_dir / "counts.npy", counts)
    print(f"Wrote counts.npy {counts.shape} → {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-dir", required=True, type=Path,
                    help="Folder containing {ACC}.counts.tsv files")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="TSV with 'accession' column (subset to include)")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="Output directory for the NPY store")
    ap.add_argument("--keep-contigs", required=True,
                    help="Comma-separated contig names to include")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    keep_contigs = {c.strip() for c in args.keep_contigs.split(",") if c.strip()}
    print(f"Keep contigs: {sorted(keep_contigs)}")

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    if "accession" not in manifest.columns:
        sys.exit(f"ERROR: --manifest missing 'accession' column: {args.manifest}")

    # Filter to samples that actually have a count file on disk
    sample_ids = manifest["accession"].tolist()
    have = sorted(
        sid for sid in sample_ids
        if (args.counts_dir / f"{sid}.counts.tsv").exists()
    )
    missing = sorted(set(sample_ids) - set(have))
    if missing:
        print(f"WARNING: {len(missing):,} samples from manifest have no count "
              f"file in {args.counts_dir} — skipping. First few: {missing[:5]}",
              file=sys.stderr)
    if not have:
        sys.exit("ERROR: no samples have count files; aborting.")
    print(f"Building store from {len(have):,} samples "
          f"(manifest had {len(sample_ids):,}).")

    build_store(args.counts_dir, have, keep_contigs, args.out_dir, args.workers)


if __name__ == "__main__":
    main()
