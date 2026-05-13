#!/usr/bin/env python3
"""
Select a stratified, representative subset of the KpSC expansion cohort
for FASTQ download / CNVRock training.

Output:
    assets/kpsc_expansion_subset_5k.tsv     — selected (accession, layout, r1_url, r2_url)
    assets/kpsc_expansion_subset_5k_meta.tsv — selected accession with key Kleborate cols
                                                (for downstream stratified eval)

Inputs:
    assets/kpsc_expansion_kleborate_gt.tsv  — 88,128 rows, Kleborate v3 ground truth
                                                (one row per BioSample, identified by `strain`)
    assets/kpsc_expansion_metadata.tsv      — bridges BioSample → run accessions
                                                (sample_id may be comma-separated for multi-run BioSamples)
    assets/ena_url_manifest.tsv             — 80,496 rows, FASTQ URLs (one row per run)

Stratification:
    1. Restrict to KpSC core species: K. pneumoniae, K. quasipneumoniae (both subsp.),
       K. variicola subsp. variicola, K. africana.
       (Drops ~2,500 non-KpSC samples — E. coli, K. oxytoca, etc.)
    2. Bridge BioSample → run accession via metadata.sample_id (first run only when
       a BioSample has multiple runs — keeps one FASTQ pair per sample).
    3. Inner-join with ENA manifest (must have FASTQ URLs).
    4. Stratify by species × Bla_Carb_acquired presence (carbapenemase carrier yes/no)
       to keep carbapenemase carriers (CNVRock's primary signal) well-represented.
    5. Within each stratum, cap any single ST at MAX_PER_ST samples
       (avoid overweighting common clones like ST258).
    6. Sample within strata proportional to stratum size, biased 1.5× toward
       carbapenemase carriers.

Reproducibility: deterministic via fixed random seed (--seed, default 42).
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np

KPSC_CORE_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae subsp. quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. similipneumoniae",
    "Klebsiella variicola subsp. variicola",
    "Klebsiella variicola subsp. tropica",
    "Klebsiella africana",
}

MAX_PER_ST = 150       # cap on samples from any single ST
CARB_OVERSAMPLE = 1.5  # weight applied to carbapenemase-carrying samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kleborate", default="assets/kpsc_expansion_kleborate_gt.tsv")
    ap.add_argument("--metadata", default="assets/kpsc_expansion_metadata.tsv")
    ap.add_argument("--ena-manifest", default="assets/ena_url_manifest.tsv")
    ap.add_argument("--target-n", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-tsv", default="assets/kpsc_expansion_subset_5k.tsv")
    ap.add_argument("--out-meta", default="assets/kpsc_expansion_subset_5k_meta.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"Loading Kleborate ground truth: {args.kleborate}")
    kleb = pd.read_csv(args.kleborate, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(kleb):,} rows")

    print(f"Loading metadata bridge: {args.metadata}")
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    # If multiple runs per BioSample (comma-separated), keep the first.
    meta["run_acc"] = meta["sample_id"].str.split(",").str[0]
    print(f"  {len(meta):,} rows  (bridges BioSample → run_acc)")

    print(f"Loading ENA URL manifest: {args.ena_manifest}")
    ena = pd.read_csv(args.ena_manifest, sep="\t", dtype=str)
    print(f"  {len(ena):,} rows")

    # 1. Restrict to KpSC core species
    in_kpsc = kleb["species"].isin(KPSC_CORE_SPECIES)
    print(f"\nKpSC core species filter:  {in_kpsc.sum():,} / {len(kleb):,} retained")
    kleb = kleb[in_kpsc].copy()

    # 2. Attach run accession via metadata bridge
    kleb = kleb.merge(meta[["sample_accession", "run_acc"]],
                      left_on="strain", right_on="sample_accession", how="inner")
    print(f"After BioSample→run join:  {len(kleb):,}")

    # 3. Inner join with ENA URL manifest
    merged = kleb.merge(ena, left_on="run_acc", right_on="accession", how="inner")
    print(f"After ENA manifest join:   {len(merged):,} samples have FASTQ URLs")

    # 3. Compute strata: species × has_carb
    merged["has_carb"] = (
        merged["Bla_Carb_acquired"].fillna("-").ne("-")
        & merged["Bla_Carb_acquired"].fillna("-").ne("")
    )
    print(f"\nCarbapenemase carriers: {merged['has_carb'].sum():,} / {len(merged):,}")
    print("\nStratum sizes:")
    print(merged.groupby(["species", "has_carb"]).size())

    # 4. Cap each ST. Replace blank STs with sentinel so groupby works.
    merged["ST_filled"] = merged["ST"].fillna("UNKNOWN").replace("", "UNKNOWN")

    # 5. Stratified sample with ST cap and carb oversampling
    selected_idx = []
    for (species, has_carb), group in merged.groupby(["species", "has_carb"]):
        # Cap each ST within stratum
        capped_idx = []
        for st, st_group in group.groupby("ST_filled"):
            if len(st_group) > MAX_PER_ST:
                pick = rng.choice(st_group.index, size=MAX_PER_ST, replace=False)
                capped_idx.extend(pick)
            else:
                capped_idx.extend(st_group.index)
        selected_idx.extend(
            (idx, CARB_OVERSAMPLE if has_carb else 1.0) for idx in capped_idx
        )

    # Build pool with weights
    pool_idx = np.array([i for i, _ in selected_idx])
    pool_w = np.array([w for _, w in selected_idx], dtype=float)
    pool_w = pool_w / pool_w.sum()  # rebalance to sum to 1

    target = min(args.target_n, len(pool_idx))
    if target == 0:
        print("ERROR: empty sampling pool — check joins", file=sys.stderr)
        sys.exit(1)
    print(f"\nPost-cap pool: {len(pool_idx):,}  →  sampling {target:,}")
    chosen = rng.choice(pool_idx, size=target, replace=False, p=pool_w)
    subset = merged.loc[chosen].copy()

    # 6. Verify representativity
    print("\nSubset breakdown (species × has_carb):")
    print(subset.groupby(["species", "has_carb"]).size())
    print(f"\nSubset unique STs: {subset['ST_filled'].nunique():,}")
    print(f"Carb carriers in subset: {subset['has_carb'].sum():,} ({100*subset['has_carb'].mean():.1f}%)")

    # 7. Write outputs
    out_dl = subset[["accession", "layout", "r1_url", "r2_url"]].sort_values("accession")
    out_dl.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"\nWrote download manifest:   {args.out_tsv}  ({len(out_dl):,} rows)")

    meta_cols = [
        "accession", "species", "ST", "Bla_acquired",
        "Bla_ESBL_acquired", "Bla_Carb_acquired", "AGly_acquired",
        "Flq_acquired", "Tet_acquired", "Sul_acquired",
        "virulence_score", "resistance_score",
    ]
    meta_cols = [c for c in meta_cols if c in subset.columns]
    subset[meta_cols].sort_values("accession").to_csv(
        args.out_meta, sep="\t", index=False
    )
    print(f"Wrote subset metadata:     {args.out_meta}")


if __name__ == "__main__":
    main()
