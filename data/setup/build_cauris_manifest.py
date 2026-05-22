#!/usr/bin/env python3
"""
Build the Candida auris sample manifest for CNVRock.

Unlike the bacterial cohorts there is no CABBAGE AST table for a fungus, so
the manifest is simply every public C. auris paired-end Illumina run; the
phenotype layer (ERG11/TAC1b/FKS1 resistance mutations + clade) is derived
downstream from the same WGS, not from an external AST file (see
paper/cauris_experiment_skeleton.md).

  1. Query the ENA portal for all C. auris (NCBI taxon 498019) read_run
     records.
  2. Keep PAIRED Illumina runs with FASTQ available; one run per BioSample
     (deepest by read_count).
  3. Write assets/cauris_subset.tsv (accession, layout, r1_url, r2_url,
     biosample, read_count) — the format the Aspera pipeline consumes.

Usage (needs internet — run on the HPC login node):
    python3 data/setup/build_cauris_manifest.py
"""
from __future__ import annotations

import io
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "assets/cauris_subset.tsv"

ENA_TAXON = 498019  # Candida auris (Candidozyma auris)
ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/search"
    f"?result=read_run&query=tax_eq({ENA_TAXON})"
    "&fields=run_accession,sample_accession,library_layout,library_strategy,"
    "fastq_ftp,read_count,instrument_platform&format=tsv&limit=0"
)


def split_fastq(ftp):
    parts = [p for p in str(ftp).split(";") if p]
    urls = ["https://" + p if not p.startswith("http") else p for p in parts]
    r1 = next((u for u in urls if "_1.fastq" in u), "")
    r2 = next((u for u in urls if "_2.fastq" in u), "")
    return r1, r2


def main():
    print(f"Querying ENA for all C. auris (taxon {ENA_TAXON}) read runs …")
    with urlopen(ENA_URL, timeout=600) as resp:
        ena = pd.read_csv(io.BytesIO(resp.read()), sep="\t", dtype=str)
    print(f"  ENA returned {len(ena):,} read_run records")

    runs = ena[
        (ena["library_layout"] == "PAIRED")
        & (ena["fastq_ftp"].notna() & (ena["fastq_ftp"] != ""))
        & (ena["instrument_platform"].str.upper() == "ILLUMINA")
        & (ena["library_strategy"].str.upper() == "WGS")
    ].copy()
    runs["read_count"] = pd.to_numeric(runs["read_count"], errors="coerce").fillna(0)
    print(f"  PAIRED Illumina WGS runs with FASTQ: {len(runs):,}")

    runs = (runs.sort_values("read_count", ascending=False)
                .drop_duplicates(subset=["sample_accession"]))
    print(f"  unique BioSamples: {len(runs):,}")

    r1r2 = runs["fastq_ftp"].map(split_fastq)
    runs["r1_url"] = [x[0] for x in r1r2]
    runs["r2_url"] = [x[1] for x in r1r2]
    runs = runs[(runs["r1_url"] != "") & (runs["r2_url"] != "")]
    print(f"  with resolvable paired FASTQ URLs: {len(runs):,}")

    manifest = pd.DataFrame({
        "accession":  runs["run_accession"].values,
        "layout":     "PAIRED",
        "r1_url":     runs["r1_url"].values,
        "r2_url":     runs["r2_url"].values,
        "biosample":  runs["sample_accession"].values,
        "read_count": runs["read_count"].astype(int).values,
    }).sort_values("accession")
    manifest.to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT}  ({len(manifest):,} samples)")


if __name__ == "__main__":
    main()
