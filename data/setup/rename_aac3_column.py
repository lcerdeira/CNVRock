#!/usr/bin/env python3
"""
One-time migration: rename the aac3-IIa gene column to aac3-II in the
plasmid NPY data store and in plasmid_gene_coords.tsv.

Run ONCE on HPC after pulling the updated source files:
  python3 data/setup/rename_aac3_column.py

Idempotent: if aac3-II already exists (and aac3-IIa does not), exits cleanly.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]

# The two NPY stores that may have the old column
STORES = [
    REPO / "data/inputs/KpSC-plasmid-1000bp-npy",
    REPO / "data/inputs/KpSC-expansion-plasmid-1000bp-npy",
]

COORDS = REPO / "assets/plasmid_refs/plasmid_gene_coords.tsv"
FASTA  = REPO / "assets/HS11286_extended.fasta"

OLD = "aac3-IIa"
NEW = "aac3-II"


def rename_npy_store(store: Path) -> None:
    meta_path = store / "meta.tsv"
    if not meta_path.exists():
        # Try .npz-style: look for a genes_index.tsv / columns.tsv
        for alt in ["genes_index.tsv", "columns.tsv", "gene_index.tsv"]:
            alt_path = store / alt
            if alt_path.exists():
                meta_path = alt_path
                break
        else:
            print(f"  {store.name}: no meta file found — skipping")
            return

    meta = pd.read_csv(meta_path, sep="\t")
    col_candidates = [c for c in meta.columns if "gene" in c.lower() or c == "name"]
    gene_col = col_candidates[0] if col_candidates else meta.columns[0]

    if NEW in meta[gene_col].values:
        print(f"  {store.name}: '{NEW}' already present — nothing to do")
        return
    if OLD not in meta[gene_col].values:
        print(f"  {store.name}: '{OLD}' not found — nothing to rename")
        return

    # Back up meta
    bak = meta_path.with_suffix(".tsv.bak_aac3rename")
    shutil.copy2(meta_path, bak)

    meta[gene_col] = meta[gene_col].str.replace(
        f"^{OLD}$", NEW, regex=True)
    meta.to_csv(meta_path, sep="\t", index=False)
    print(f"  {store.name}: renamed '{OLD}' -> '{NEW}' in {meta_path.name}")

    # Rename any per-gene npy files (pattern: aac3-IIa.npy → aac3-II.npy)
    old_npy = store / f"{OLD}.npy"
    new_npy = store / f"{NEW}.npy"
    if old_npy.exists() and not new_npy.exists():
        old_npy.rename(new_npy)
        print(f"    renamed {old_npy.name} -> {new_npy.name}")


def rename_coords(coords: Path) -> None:
    if not coords.exists():
        print(f"  {coords}: not found — skipping")
        return
    df = pd.read_csv(coords, sep="\t", dtype=str)
    if NEW in df["gene"].values:
        print(f"  {coords.name}: '{NEW}' already present — nothing to do")
        return
    if OLD not in df["gene"].values:
        print(f"  {coords.name}: '{OLD}' not found — nothing to rename")
        return
    bak = coords.with_suffix(".tsv.bak_aac3rename")
    shutil.copy2(coords, bak)
    df["gene"] = df["gene"].str.replace(f"^{OLD}$", NEW, regex=True)
    df.to_csv(coords, sep="\t", index=False)
    print(f"  {coords.name}: renamed '{OLD}' -> '{NEW}'")


def rename_fasta_header(fasta: Path) -> None:
    if not fasta.exists():
        print(f"  {fasta.name}: not found — skipping")
        return
    text = fasta.read_text()
    if f">{NEW}" in text:
        print(f"  {fasta.name}: >{NEW} already present — nothing to do")
        return
    if f">{OLD}" not in text:
        print(f"  {fasta.name}: >{OLD} not found — nothing to rename")
        return
    bak = fasta.with_suffix(".fasta.bak_aac3rename")
    shutil.copy2(fasta, bak)
    fasta.write_text(text.replace(f">{OLD}", f">{NEW}", 1))
    print(f"  {fasta.name}: renamed >{OLD} -> >{NEW} header")


def main() -> None:
    print(f"Renaming '{OLD}' → '{NEW}' across data pipeline…\n")

    for store in STORES:
        if store.exists():
            rename_npy_store(store)
        else:
            print(f"  {store.name}: store dir not found — skipping")

    rename_coords(COORDS)
    rename_fasta_header(FASTA)

    print("\nDone. Re-run evaluation (exp 37 / exp 33) after this patch.")


if __name__ == "__main__":
    main()
