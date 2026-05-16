#!/usr/bin/env python3
"""
Build a plasmid-gene NPY store from KpSC 1 kb count TSVs.

Different from `readcounts_to_npy_kpsc.py` (which keeps every 1 kb bin):
this script AGGREGATES 1 kb bins by gene region defined in
`plasmid_gene_coords.tsv`, producing one count per gene per sample. This
matches the existing CNVRock plasmid store schema (12 bins / contigs entries,
each representing one AMR gene of interest).

For each sample × gene:
    count = sum of `COUNT` across all 1 kb bins where
            bin.CHROM == gene.contig
            AND bin.END > gene.start
            AND bin.START < gene.end
(overlap test, half-open intervals).

Inputs
------
    --counts-dir DIR   {ACC}.counts.tsv files (1 kb bins for chrom + plasmids)
    --manifest FILE    TSV with `accession` column (sample subset)
    --gene-coords FILE plasmid_gene_coords.tsv (gene, contig, start, end, …)
    --out-dir DIR      Output store directory
    --workers N        Parallel readers (default 16)

Outputs (CNVRock plasmid store schema)
--------------------------------------
    contigs.npy     structured array [(chrom, start, end)] per gene
    counts.npy      uint32 (n_samples, n_genes)
    sample_ids.npy  object array of accessions
    genes.npy       object array of gene names (NEW vs old store — for
                     traceability; CNVRock loader ignores extras)
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

REQUIRED_BIN_COLS = {"CONTIG", "START", "END", "COUNT"}
REQUIRED_GENE_COLS = {"gene", "contig", "start", "end"}


def read_bins(path: Path) -> pd.DataFrame:
    """Parse a GATK CollectReadCounts TSV. Returns CHROM/START/END/COUNT."""
    if not path.exists():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")

    with open(path) as f:
        skip = sum(1 for line in f if line.startswith("@"))
    df = pd.read_csv(path, sep="\t", skiprows=skip)

    missing = REQUIRED_BIN_COLS - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    df = df.rename(columns={"CONTIG": "CHROM"})
    df["START"] = df["START"].astype(np.uint32)
    df["END"]   = df["END"].astype(np.uint32)
    df["COUNT"] = df["COUNT"].astype(np.uint32)
    return df[["CHROM", "START", "END", "COUNT"]]


def aggregate_one(path: Path, genes: pd.DataFrame) -> np.ndarray:
    """Return per-gene summed COUNT for one sample (shape: (n_genes,))."""
    bins = read_bins(path)
    out = np.zeros(len(genes), dtype=np.uint32)
    # Index bins by chrom once for speed
    by_chrom = {c: g.reset_index(drop=True) for c, g in bins.groupby("CHROM")}
    for i, row in enumerate(genes.itertuples(index=False)):
        chrom = row.contig
        if chrom not in by_chrom:
            continue                                   # leave as 0
        sub = by_chrom[chrom]
        # Overlap: bin.END > gene.start AND bin.START < gene.end
        mask = (sub["END"] > row.start) & (sub["START"] < row.end)
        out[i] = sub.loc[mask, "COUNT"].sum()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--gene-coords", required=True, type=Path,
                    help="plasmid_gene_coords.tsv")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    # ── Gene coords ───────────────────────────────────────────────────────────
    genes = pd.read_csv(args.gene_coords, sep="\t")
    missing = REQUIRED_GENE_COLS - set(genes.columns)
    if missing:
        sys.exit(f"ERROR: gene coords missing columns {missing}")
    genes["start"] = genes["start"].astype(np.uint32)
    genes["end"]   = genes["end"].astype(np.uint32)
    n_genes = len(genes)
    print(f"Loaded {n_genes} gene targets:")
    for r in genes.itertuples(index=False):
        print(f"  {r.gene:<14s} {r.contig:<16s} {r.start:>9d}–{r.end:<9d}")

    # ── Manifest → sample list ────────────────────────────────────────────────
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    if "accession" not in manifest.columns:
        sys.exit(f"ERROR: manifest missing 'accession' column: {args.manifest}")

    sample_ids = manifest["accession"].tolist()
    have = sorted(
        sid for sid in sample_ids
        if (args.counts_dir / f"{sid}.counts.tsv").exists()
    )
    missing = sorted(set(sample_ids) - set(have))
    if missing:
        print(f"WARNING: {len(missing):,} samples have no count file — skipping. "
              f"First few: {missing[:5]}", file=sys.stderr)
    if not have:
        sys.exit("ERROR: no samples have count files; aborting.")
    n_samples = len(have)
    print(f"\nBuilding store from {n_samples:,} samples (manifest had "
          f"{len(sample_ids):,}).")

    # ── Aggregate ─────────────────────────────────────────────────────────────
    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts = np.zeros((n_samples, n_genes), dtype=np.uint32)

    def one(args_):
        idx, sid = args_
        try:
            return idx, aggregate_one(args.counts_dir / f"{sid}.counts.tsv", genes)
        except Exception as exc:                      # noqa: BLE001
            raise RuntimeError(f"Failed sample '{sid}': {exc}") from exc

    jobs = list(enumerate(have))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for idx, vec in tqdm(pool.map(one, jobs), total=n_samples,
                             desc="Aggregating"):
            counts[idx, :] = vec

    # ── Write store ───────────────────────────────────────────────────────────
    contigs = np.empty(
        n_genes,
        dtype=[("chrom", object), ("start", np.uint32), ("end", np.uint32)],
    )
    contigs["chrom"] = genes["contig"].values
    contigs["start"] = genes["start"].values
    contigs["end"]   = genes["end"].values
    np.save(args.out_dir / "contigs.npy", contigs)
    np.save(args.out_dir / "sample_ids.npy",
            np.array(have, dtype=object))
    np.save(args.out_dir / "counts.npy", counts)
    np.save(args.out_dir / "genes.npy",
            np.array(genes["gene"].tolist(), dtype=object))

    print(f"\nDone. counts shape: {counts.shape}  →  {args.out_dir}")
    print("Per-gene mean depth:")
    for i, gene in enumerate(genes["gene"]):
        nz = (counts[:, i] > 0).sum()
        print(f"  {gene:<14s} mean={counts[:, i].mean():8.1f}  "
              f"non-zero={nz:,}/{n_samples:,}")


if __name__ == "__main__":
    main()
