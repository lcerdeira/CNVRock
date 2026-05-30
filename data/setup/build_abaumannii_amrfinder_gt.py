#!/usr/bin/env python3
"""
Build AMRFinder ground-truth for the A. baumannii cohort.

Source: AllTheBacteria pre-computed AMRFinderPlus parquet
        (~1.1 GB, ~/Library/Application Support/atb/data/amrfinderplus.parquet)

Cohort: assets/abaumannii_subset.tsv  (1,454 samples, biosample + accession)

Output: assets/abaumannii_amrfinder_gt.tsv

Evaluable genes (those that are absent in some strains):
  blaOXA-23       — acquired OXA-type carbapenemase; chromosomal in AB5075-UW
                    but absent in many other A. baumannii (CN=0 expected).
                    ~32% prevalence in this cohort.
  blaOXA-24-like  — acquired carbapenemase (OXA-24/40/72 family); Subclass=CARBAPENEM.
  blaOXA-58-like  — acquired carbapenemase (OXA-58 family);  Subclass=CARBAPENEM.

NOT evaluable as presence GT (intrinsic, present in all A. baumannii):
  blaOXA-51-like  — intrinsic OXA carbapenemase (OXA-66, OXA-69, OXA-71 …)
  adeA/B/C/R/S    — AdeABC efflux system (chromosomal, intrinsic)
  adeI/J/K        — AdeIJK efflux system (chromosomal, intrinsic)

Run locally (ATB parquet required):
  python3 data/setup/build_abaumannii_amrfinder_gt.py
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ATB = Path(os.path.expanduser(
    "~/Library/Application Support/atb/data/amrfinderplus.parquet"))
MANIFEST = REPO / "assets/abaumannii_subset.tsv"
OUT = REPO / "assets/abaumannii_amrfinder_gt.tsv"

# OXA-51-like intrinsic variants present in almost all A. baumannii — used as
# a denominator-check but NOT output as GT columns.
OXA51_LIKE = {"blaOXA-51", "blaOXA-66", "blaOXA-69", "blaOXA-64", "blaOXA-65",
              "blaOXA-71", "blaOXA-82", "blaOXA-90", "blaOXA-91", "blaOXA-92",
              "blaOXA-94", "blaOXA-95", "blaOXA-98", "blaOXA-100", "blaOXA-104",
              "blaOXA-106", "blaOXA-109", "blaOXA-113", "blaOXA-120", "blaOXA-121",
              "blaOXA-126", "blaOXA-132", "blaOXA-144", "blaOXA-172", "blaOXA-208",
              "blaOXA-223", "blaOXA-237", "blaOXA-252", "blaOXA-254", "blaOXA-259",
              "blaOXA-263", "blaOXA-312", "blaOXA-313", "blaOXA-314", "blaOXA-316",
              "blaOXA-317", "blaOXA-340", "blaOXA-371", "blaOXA-374", "blaOXA-378",
              "blaOXA-402", "blaOXA-407", "blaOXA-414", "blaOXA-441", "blaOXA-500",
              "blaOXA-510", "blaOXA-523", "blaOXA-525", "blaOXA-528", "blaOXA-529",
              "blaOXA-530", "blaOXA-531", "blaOXA-532"}

# Genes to include in the GT file.
# Key   = output column name (must match gene name in gene_calls.tsv)
# Value = callable that returns True for a given (element_symbol, subclass) pair
GENE_DEFS = {
    # Exact match — blaOXA-23 is a single well-defined variant
    "blaOXA-23": lambda sym, sub: sym == "blaOXA-23",
    # OXA-24/40/72 family — acquired carbapenems, NOT OXA-51-like
    "blaOXA-24-like": lambda sym, sub: (
        sub == "CARBAPENEM"
        and sym.startswith("blaOXA")
        and sym not in OXA51_LIKE
        and sym != "blaOXA-23"
        and sym != "blaOXA-58"
        and any(sym.startswith(p) for p in [
            "blaOXA-24", "blaOXA-40", "blaOXA-72",
        ])
    ),
    # OXA-58 family
    "blaOXA-58-like": lambda sym, sub: sym.startswith("blaOXA-58"),
}


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    bs_to_acc = dict(zip(manifest["biosample"], manifest["accession"]))
    bs_set = set(bs_to_acc.keys())
    print(f"A. baumannii cohort: {len(bs_set):,} biosamples")

    print("Loading ATB AMRFinderPlus parquet …")
    amr = pd.read_parquet(
        ATB, columns=["Name", "Element symbol", "Subclass"])
    abr = amr[amr["Name"].isin(bs_set)].copy()
    print(f"ATB rows for cohort: {len(abr):,} across {abr['Name'].nunique():,} biosamples")

    results: dict[str, pd.Series] = {"biosample": sorted(bs_set)}
    for gene, pred in GENE_DEFS.items():
        mask = abr.apply(
            lambda row: pred(row["Element symbol"], row["Subclass"]), axis=1)
        hits = abr[mask].groupby("Name").size().gt(0).astype(int)
        results[gene] = (
            pd.Series(results["biosample"])
            .map(hits)
            .fillna(0)
            .astype(int)
            .values
        )

    df = pd.DataFrame(results)
    df.insert(0, "sample_id", df["biosample"].map(bs_to_acc))

    # Samples with no AMRFinder rows at all → already 0 for all genes; confirm.
    n_no_rows = len(bs_set) - abr["Name"].nunique()
    if n_no_rows:
        print(f"  {n_no_rows} biosamples had no AMRFinder rows (all genes → 0)")

    df = df.sort_values("sample_id").reset_index(drop=True)
    df.to_csv(OUT, sep="\t", index=False)

    print(f"\nWrote {OUT}  ({len(df):,} rows)")
    print("\nGene prevalence:")
    for gene in GENE_DEFS:
        n = int(df[gene].sum())
        pct = 100 * n / len(df)
        print(f"  {gene:<20s}: {n:5,} / {len(df):,}  ({pct:.1f}%)")


if __name__ == "__main__":
    main()
