#!/usr/bin/env python3
"""
Build the Acinetobacter baumannii sample manifest for CNVRock.

The CABBAGE A. baumannii AST table gives BioSample IDs (2,369 unique) but
not sequencing-run accessions. This script:

  1. Queries the ENA portal API for every A. baumannii (NCBI taxon 470)
     read_run record — run_accession, sample_accession, library_layout,
     fastq_ftp, read_count.
  2. Intersects on BioSample with the CABBAGE A. baumannii table.
  3. Keeps PAIRED runs with FASTQ available; one run per BioSample
     (highest read_count).
  4. Writes assets/abaumannii_subset.tsv with an `accession` column —
     the same format the Aspera download pipeline consumes.

Usage (needs internet — run on the HPC login node):
    python3 data/setup/build_abaumannii_manifest.py
"""
from __future__ import annotations

import io
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
CABBAGE = REPO / "assets/cabbage_abaumannii_phenotypes.tsv"
OUT = REPO / "assets/abaumannii_subset.tsv"

ENA_TAXON = 470  # Acinetobacter baumannii
ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/search"
    f"?result=read_run&query=tax_eq({ENA_TAXON})"
    "&fields=run_accession,sample_accession,secondary_sample_accession,"
    "library_layout,library_strategy,fastq_ftp,read_count,instrument_platform"
    "&format=tsv&limit=0"
)


def main():
    print(f"Querying ENA for all A. baumannii (taxon {ENA_TAXON}) read runs …")
    with urlopen(ENA_URL, timeout=600) as resp:
        ena = pd.read_csv(io.BytesIO(resp.read()), sep="\t", dtype=str)
    print(f"  ENA returned {len(ena):,} read_run records")

    cab = pd.read_csv(CABBAGE, sep="\t", dtype=str)
    biosamples = set(cab["BioSample_ID"].dropna())
    secondary = set(cab["SRA_accession"].dropna())
    print(f"  CABBAGE A. baumannii: {len(biosamples):,} BioSamples, "
          f"{len(secondary):,} secondary accessions")

    # match on either the BioSample (SAMEA/SAMN) or the secondary (ERS/SRS)
    hit = (ena["sample_accession"].isin(biosamples) |
           ena["secondary_sample_accession"].isin(secondary))
    runs = ena[hit].copy()
    print(f"  ENA runs matching a CABBAGE sample: {len(runs):,}")

    # Illumina paired short reads with FASTQ available
    runs = runs[
        (runs["library_layout"] == "PAIRED")
        & (runs["fastq_ftp"].notna() & (runs["fastq_ftp"] != ""))
        & (runs["instrument_platform"].str.upper() == "ILLUMINA")
    ].copy()
    runs["read_count"] = pd.to_numeric(runs["read_count"], errors="coerce").fillna(0)
    print(f"  PAIRED Illumina runs with FASTQ: {len(runs):,}")

    # one run per BioSample — the deepest
    runs = (runs.sort_values("read_count", ascending=False)
                .drop_duplicates(subset=["sample_accession"]))
    print(f"  unique BioSamples with a usable run: {len(runs):,}")

    # ENA fastq_ftp is ';'-separated ftp paths; the download pipeline expects
    # https:// R1 / R2 URL columns (accession, layout, r1_url, r2_url).
    def split_fastq(ftp):
        parts = [p for p in str(ftp).split(";") if p]
        urls = ["https://" + p if not p.startswith("http") else p
                for p in parts]
        # keep the _1 / _2 paired files; ignore a bare unpaired file
        r1 = next((u for u in urls if "_1.fastq" in u), "")
        r2 = next((u for u in urls if "_2.fastq" in u), "")
        return r1, r2

    r1r2 = runs["fastq_ftp"].map(split_fastq)
    runs["r1_url"] = [x[0] for x in r1r2]
    runs["r2_url"] = [x[1] for x in r1r2]
    runs = runs[(runs["r1_url"] != "") & (runs["r2_url"] != "")]
    print(f"  runs with resolvable paired FASTQ URLs: {len(runs):,}")

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
