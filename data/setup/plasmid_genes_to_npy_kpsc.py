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
    --counts-dir DIR    {ACC}.counts.tsv files (1 kb bins for chrom + plasmids)
    --manifest FILE     TSV with `accession` column (sample subset)
    --gene-coords FILE  plasmid_gene_coords.tsv (gene, contig, start, end, …)
    --families FILE     Optional plasmid_gene_families.tsv (family, members)
                         When given, output is per-FAMILY (sum bins across all
                         allele variants in each family). Without it, output is
                         per-individual-gene.
    --out-dir DIR       Output store directory
    --workers N         Parallel readers (default 16)

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
    ap.add_argument("--families", default=None, type=Path,
                    help="Optional plasmid_gene_families.tsv. If given, "
                         "output store is per-FAMILY (sum across allele "
                         "variants). Without it, output is per individual gene.")
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

    # ── Gene families (optional) ──────────────────────────────────────────────
    family_membership: dict[str, list[int]] = {}    # family_name -> [gene row idxs]
    family_order: list[str] = []
    if args.families is not None:
        fam_df = pd.read_csv(args.families, sep="\t")
        if "family" not in fam_df.columns or "members" not in fam_df.columns:
            sys.exit("ERROR: --families file needs `family` and `members` columns")
        gene_to_row = {g: i for i, g in enumerate(genes["gene"].tolist())}
        for r in fam_df.itertuples(index=False):
            members = [m.strip() for m in str(r.members).split(",") if m.strip()]
            idxs = [gene_to_row[m] for m in members if m in gene_to_row]
            if not idxs:
                print(f"  WARNING: family {r.family} has no members in "
                      f"gene-coords; skipping", file=sys.stderr)
                continue
            family_membership[r.family] = idxs
            family_order.append(r.family)
        print(f"\nLoaded {len(family_order)} families:")
        for fam in family_order:
            members = [genes["gene"].iloc[i] for i in family_membership[fam]]
            print(f"  {fam:<18s} ← {', '.join(members)}")

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

    # ── Aggregate by family if requested ──────────────────────────────────────
    if family_membership:
        n_fam = len(family_order)
        fam_counts = np.zeros((n_samples, n_fam), dtype=np.uint32)
        # contigs.npy needs a representative (chrom, start, end) per family — we
        # use the first member's coords. CNVRock loader doesn't index into this
        # for gene-level PCN; it just needs a per-bin record.
        fam_contigs = np.empty(
            n_fam,
            dtype=[("chrom", object), ("start", np.uint32), ("end", np.uint32)],
        )
        for j, fam in enumerate(family_order):
            idxs = family_membership[fam]
            fam_counts[:, j] = counts[:, idxs].sum(axis=1)
            first = idxs[0]
            fam_contigs[j] = (genes["contig"].iloc[first],
                              genes["start"].iloc[first],
                              genes["end"].iloc[first])

        np.save(args.out_dir / "contigs.npy", fam_contigs)
        np.save(args.out_dir / "sample_ids.npy",
                np.array(have, dtype=object))
        np.save(args.out_dir / "counts.npy", fam_counts)
        np.save(args.out_dir / "genes.npy",
                np.array(family_order, dtype=object))

        print(f"\nDone. counts shape: {fam_counts.shape}  →  {args.out_dir}")
        print("Per-family mean depth:")
        for j, fam in enumerate(family_order):
            nz = int((fam_counts[:, j] > 0).sum())
            print(f"  {fam:<18s} mean={fam_counts[:, j].mean():9.1f}  "
                  f"non-zero={nz:,}/{n_samples:,}")
        return

    # ── Write per-gene store ──────────────────────────────────────────────────
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
