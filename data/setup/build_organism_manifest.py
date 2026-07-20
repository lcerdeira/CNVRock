#!/usr/bin/env python3
"""
Build a CNVRock sample manifest for an arbitrary organism.

Generalises build_abaumannii_manifest.py so E. coli and S. aureus (and any
later organism) reuse one code path. For each organism it:

  1. Queries the ENA portal API for read_run records, filtering server-side
     to PAIRED Illumina WGS. E. coli alone has millions of runs on ENA, so
     filtering client-side would mean downloading a payload measured in
     hundreds of megabytes.
  2. Intersects on BioSample with the cached CABBAGE AST table.
  3. Keeps one run per BioSample (the deepest) with resolvable paired FASTQ.
  4. Stratifies by resistance status on a priority antibiotic and draws a
     seeded subsample of --n samples, so the tier is not dominated by
     whichever phenotype happens to be over-sequenced.
  5. Writes assets/<name>_subset.tsv in the format the Aspera download
     pipeline consumes.

Stratification is by phenotype rather than by ST because MLST is not known
until the assemblies are in hand, whereas AST is available at manifest time.

Usage (needs internet):
    python3 data/setup/build_organism_manifest.py --organism ecoli   --n 3000
    python3 data/setup/build_organism_manifest.py --organism saureus --n 3000
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

import numpy as np
import pandas as pd

REPO    = Path(__file__).resolve().parents[2]
CABBAGE = REPO / "assets/cabbage_cache/phenotype_2025-12.parquet"

ORGANISMS = {
    "ecoli":   dict(taxon=562,  cabbage="Escherichia coli",
                    out="ecoli_subset.tsv",   priority="ceftriaxone"),
    "saureus": dict(taxon=1280, cabbage="Staphylococcus aureus",
                    out="saureus_subset.tsv", priority="oxacillin"),
}

FIELDS = ("run_accession,sample_accession,secondary_sample_accession,"
          "library_layout,library_strategy,fastq_ftp,read_count,"
          "instrument_platform")


def query_ena(taxon: int) -> pd.DataFrame:
    q = (f'tax_eq({taxon}) AND library_layout="PAIRED"'
         ' AND instrument_platform="ILLUMINA" AND library_strategy="WGS"')
    url = ("https://www.ebi.ac.uk/ena/portal/api/search"
           f"?result=read_run&query={quote(q)}&fields={FIELDS}"
           "&format=tsv&limit=0")
    print(f"  querying ENA (taxon {taxon}, PAIRED Illumina WGS) …", flush=True)
    with urlopen(url, timeout=1800) as resp:
        return pd.read_csv(io.BytesIO(resp.read()), sep="\t", dtype=str)


def split_fastq(ftp: str) -> tuple[str, str]:
    urls = ["https://" + p if not p.startswith("http") else p
            for p in str(ftp).split(";") if p]
    r1 = next((u for u in urls if "_1.fastq" in u), "")
    r2 = next((u for u in urls if "_2.fastq" in u), "")
    return r1, r2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, choices=sorted(ORGANISMS))
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    cfg = ORGANISMS[args.organism]
    out = REPO / "assets" / cfg["out"]

    ena = query_ena(cfg["taxon"])
    print(f"  ENA returned {len(ena):,} read_run records")

    cab = pd.read_parquet(CABBAGE)
    cab = cab[cab["organism"].astype(str).str.contains(cfg["cabbage"], na=False)]
    biosamples = set(cab["BioSample_ID"].dropna())
    print(f"  CABBAGE {cfg['cabbage']}: {len(biosamples):,} BioSamples with AST")

    runs = ena[ena["sample_accession"].isin(biosamples)].copy()
    print(f"  ENA runs matching a CABBAGE sample: {len(runs):,}")

    runs = runs[runs["fastq_ftp"].notna() & (runs["fastq_ftp"] != "")].copy()
    runs["read_count"] = pd.to_numeric(runs["read_count"],
                                       errors="coerce").fillna(0)
    runs = (runs.sort_values("read_count", ascending=False)
                .drop_duplicates(subset=["sample_accession"]))
    r1r2 = runs["fastq_ftp"].map(split_fastq)
    runs["r1_url"] = [x[0] for x in r1r2]
    runs["r2_url"] = [x[1] for x in r1r2]
    runs = runs[(runs["r1_url"] != "") & (runs["r2_url"] != "")]
    print(f"  unique BioSamples with resolvable paired FASTQ: {len(runs):,}")

    # ── stratify by phenotype on the priority antibiotic ────────────────────
    pri = cab[cab["antibiotic_name"].astype(str).str.lower()
              == cfg["priority"]][["BioSample_ID", "resistance_phenotype"]]
    pri = pri.drop_duplicates(subset=["BioSample_ID"])
    runs = runs.merge(pri, left_on="sample_accession",
                      right_on="BioSample_ID", how="left")
    runs["stratum"] = (runs["resistance_phenotype"].astype(str)
                       .str.lower().str[:1].replace({"n": "u"}).fillna("u"))
    print("  strata (%s): %s" % (cfg["priority"],
          runs["stratum"].value_counts().to_dict()))

    if len(runs) > args.n:
        rng = np.random.default_rng(args.seed)
        frac = args.n / len(runs)
        picks = []
        for _, grp in runs.groupby("stratum", sort=True):
            k = max(1, int(round(len(grp) * frac)))
            idx = rng.choice(len(grp), size=min(k, len(grp)), replace=False)
            picks.append(grp.iloc[np.sort(idx)])
        runs = pd.concat(picks).head(args.n)
        print(f"  stratified subsample -> {len(runs):,}")

    manifest = pd.DataFrame({
        "accession":  runs["run_accession"].values,
        "layout":     "PAIRED",
        "r1_url":     runs["r1_url"].values,
        "r2_url":     runs["r2_url"].values,
        "biosample":  runs["sample_accession"].values,
        "read_count": runs["read_count"].astype(int).values,
        "stratum":    runs["stratum"].values,
    }).sort_values("accession")
    manifest.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}  ({len(manifest):,} samples)")


if __name__ == "__main__":
    main()
