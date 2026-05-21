#!/usr/bin/env python3
"""
Build the Acinetobacter baumannii AB5075-UW reference for CNVRock.

AB5075-UW (BioProject PRJNA254856) is the community MDR reference strain:
its chromosome carries the AdeABC and AdeIJK RND efflux operons and the
intrinsic blaOXA-51-like carbapenemase — the chromosomal CNV targets.

This script:
  1. Resolves the AB5075-UW chromosome + native plasmid nuccore accessions
     via NCBI Entrez.
  2. Downloads them and writes assets/abaumannii_ref/AB5075.fasta.
  3. Reports contig names + lengths for the interval-list / gene-coords step.

Acquired carbapenemase contigs (blaOXA-23/-24-40/-58, blaNDM-1, armA) are
appended in a separate step (build_abaumannii_extended.py), mirroring the
HS11286_extended workflow for KpSC.

Usage (cnvrock env has Biopython):
    python3 data/setup/build_abaumannii_reference.py
"""
from __future__ import annotations

import time
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = "louise.cerdeira@gmail.com"

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "assets/abaumannii_ref"
OUT_FASTA = OUT_DIR / "AB5075.fasta"

# AB5075-UW complete genome — RefSeq accessions (BioProject PRJNA254856).
# Chromosome + 3 native plasmids. (CP008710-712 are NOT AB5075 — they
# belong to other organisms — and are deliberately excluded.)
AB5075_ACCESSIONS = [
    "NZ_CP008706.1",   # chromosome (~3.84 Mb) — AdeABC/AdeIJK, blaOXA-51-like
    "NZ_CP008707.1",   # plasmid p1AB5075 (~83.6 kb, AbaR resistance island)
    "NZ_CP008708.1",   # plasmid p2AB5075 (~8.7 kb)
    "NZ_CP008709.1",   # plasmid p3AB5075 (~1.9 kb)
]


def fetch(acc: str) -> SeqIO.SeqRecord:
    for attempt in range(1, 4):
        try:
            with Entrez.efetch(db="nuccore", id=acc, rettype="fasta",
                               retmode="text") as h:
                rec = SeqIO.read(h, "fasta")
            return rec
        except Exception as exc:                       # noqa: BLE001
            print(f"  attempt {attempt} failed for {acc}: {exc}")
            time.sleep(5)
    raise RuntimeError(f"could not fetch {acc}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    for acc in AB5075_ACCESSIONS:
        print(f"Fetching {acc} ...")
        rec = fetch(acc)
        records.append(rec)
        print(f"  {rec.id}  {len(rec.seq):,} bp  — {rec.description[:70]}")
        time.sleep(1)
    SeqIO.write(records, OUT_FASTA, "fasta")
    total = sum(len(r.seq) for r in records)
    print(f"\nwrote {OUT_FASTA}")
    print(f"  {len(records)} contigs, {total:,} bp total")
    print("Next: append acquired-resistance gene contigs, then "
          "bwa index + gatk PreprocessIntervals (1 kb bins).")


if __name__ == "__main__":
    main()
