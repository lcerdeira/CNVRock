#!/usr/bin/env python3
"""
Re-annotate plasmid_gene_coords.tsv by BLAST-ing canonical AMR gene reference
sequences against HS11286_extended.fasta.

Why this exists
---------------
Diagnosed 2026-05-17: existing plasmid_gene_coords.tsv had drifted coordinates
for blaKPC-2, blaCTX-M-15, blaTEM-1, aac6-Ib-cr (and likely others). Reads
mapped correctly to the plasmid contigs at MQ=20 but to positions far from
the declared gene regions, producing zero counts for those genes. qnrB1
coordinates were correct, confirming the issue is per-gene, not pipeline-wide.

Workflow
--------
1. Fetch each gene's canonical nucleotide sequence from NCBI Entrez (one
   well-curated allele accession per gene). One-time download to
   assets/plasmid_refs/gene_references.fasta.
2. Build a BLAST database from HS11286_extended.fasta.
3. blastn each gene against the reference; pick the top-scoring hit per gene.
4. Write the corrected plasmid_gene_coords.tsv (same schema; only start/end
   change; contig may also change if a gene's best hit moves to a different
   plasmid backbone).

Usage
-----
    # First time (downloads gene refs):
    python3 data/setup/reannotate_plasmid_genes_blast.py \\
        --reference assets/HS11286_extended.fasta \\
        --out assets/plasmid_refs/plasmid_gene_coords_blast.tsv \\
        --fetch-refs

    # Subsequent runs (uses cached gene refs):
    python3 data/setup/reannotate_plasmid_genes_blast.py \\
        --reference assets/HS11286_extended.fasta \\
        --out assets/plasmid_refs/plasmid_gene_coords_blast.tsv
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ── Canonical AMR gene references (NCBI Reference Gene Catalog) ──────────────
# (gene_name, accession, expected_length_bp, default_absent_threshold)
GENE_REFS = [
    ("blaKPC-2",     "NG_049253.1",  882, 0.20),
    ("blaCTX-M-15",  "NG_048935.1",  876, 0.50),
    ("blaNDM-1",     "NG_049326.1",  813, 0.20),
    ("blaCTX-M-14",  "NG_048937.1",  876, 1.00),
    ("blaTEM-1",     "NG_050173.1",  861, 0.20),
    ("qnrB1",        "NG_050549.1",  657, 0.20),
    ("aac6-Ib-cr",   "NG_047330.1",  600, 0.10),
    ("blaOXA-48",    "NG_049731.1",  798, 0.20),
    ("blaNDM-5",     "NG_049327.1",  813, 0.20),
    ("blaOXA-181",   "NG_049838.1",  798, 0.20),
    ("blaCTX-M-65",  "NG_048988.1",  876, 0.20),
    ("blaCTX-M-27",  "NG_048948.1",  876, 0.20),
]

NCBI_EMAIL = "louise.cerdeira@gmail.com"   # required by NCBI E-utilities


def fetch_gene_references(out_fasta: Path) -> None:
    """Download canonical nucleotide sequences for each gene from NCBI Entrez."""
    try:
        from Bio import Entrez, SeqIO  # noqa: F401
    except ImportError:
        sys.exit("ERROR: Biopython required for --fetch-refs. "
                 "Run: pip install biopython")
    from Bio import Entrez

    Entrez.email = NCBI_EMAIL
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as fh:
        for gene, acc, _exp_len, _ in GENE_REFS:
            print(f"  Fetching {gene} ({acc}) …", flush=True)
            try:
                handle = Entrez.efetch(db="nucleotide", id=acc,
                                       rettype="fasta", retmode="text")
                fasta = handle.read()
                handle.close()
            except Exception as exc:                       # noqa: BLE001
                print(f"  WARNING: failed to fetch {gene} ({acc}): {exc}",
                      file=sys.stderr)
                continue
            # Replace the header with a gene-named one
            lines = fasta.splitlines()
            if lines and lines[0].startswith(">"):
                lines[0] = f">{gene}  ref={acc}"
            fh.write("\n".join(lines) + "\n")
    print(f"Wrote {out_fasta}")


def build_blast_db(reference: Path, dbdir: Path) -> Path:
    dbdir.mkdir(parents=True, exist_ok=True)
    dbprefix = dbdir / "ref"
    if not (dbdir / "ref.nhr").exists():
        print(f"Building BLAST db at {dbprefix} …")
        subprocess.run(
            ["makeblastdb", "-in", str(reference), "-dbtype", "nucl",
             "-out", str(dbprefix)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    return dbprefix


def blast_genes(query_fasta: Path, dbprefix: Path) -> dict[str, tuple]:
    """Return {gene_name: (contig, start, end, pident, length, evalue)}."""
    fmt = "6 qseqid sseqid sstart send pident length evalue bitscore"
    out = subprocess.run(
        ["blastn", "-query", str(query_fasta), "-db", str(dbprefix),
         "-outfmt", fmt, "-evalue", "1e-10", "-max_target_seqs", "5"],
        capture_output=True, text=True, check=True,
    )
    best: dict[str, tuple] = {}
    for line in out.stdout.strip().splitlines():
        qseqid, sseqid, sstart, send, pident, length, evalue, bitscore = line.split("\t")
        sstart, send = int(sstart), int(send)
        if sstart > send:                                  # antisense hit
            sstart, send = send, sstart
        pident = float(pident)
        bitscore = float(bitscore)
        prev = best.get(qseqid)
        if prev is None or bitscore > prev[-1]:
            best[qseqid] = (sseqid, sstart, send, pident, int(length),
                            float(evalue), bitscore)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path,
                    help="HS11286_extended.fasta")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output plasmid_gene_coords TSV")
    ap.add_argument("--gene-refs", default=None, type=Path,
                    help="Cached gene_references.fasta "
                         "(default: assets/plasmid_refs/gene_references.fasta)")
    ap.add_argument("--fetch-refs", action="store_true",
                    help="Download gene sequences from NCBI before BLAST.")
    ap.add_argument("--keep-blast-db", action="store_true",
                    help="Keep the BLAST DB after running (default: clean up).")
    args = ap.parse_args()

    refs_fasta = args.gene_refs or (args.reference.parent / "plasmid_refs"
                                    / "gene_references.fasta")

    if args.fetch_refs or not refs_fasta.exists():
        print(f"Fetching gene references → {refs_fasta}")
        fetch_gene_references(refs_fasta)

    if not refs_fasta.exists():
        sys.exit(f"ERROR: gene reference FASTA not found: {refs_fasta}")

    dbdir = Path(tempfile.mkdtemp(prefix="cnvrock_blast_"))
    try:
        dbprefix = build_blast_db(args.reference, dbdir)
        print("BLAST-ing genes against extended reference …")
        hits = blast_genes(refs_fasta, dbprefix)

        # Write output coords TSV (matches existing schema)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("gene\tcontig\tstart\tend\tabsent_threshold\n")
            for gene, _acc, exp_len, abs_thr in GENE_REFS:
                if gene not in hits:
                    print(f"  WARNING: no BLAST hit for {gene}",
                          file=sys.stderr)
                    continue
                contig, start, end, pident, hlen, evalue, score = hits[gene]
                length_ok = (0.85 * exp_len <= hlen <= 1.15 * exp_len)
                ident_ok = pident >= 95.0
                marker = "✓" if (length_ok and ident_ok) else "?"
                print(f"  {marker} {gene:<14s} → {contig}:{start:>9d}-{end:<9d}  "
                      f"pident={pident:5.1f}% len={hlen}/{exp_len} "
                      f"score={score:.0f}")
                fh.write(f"{gene}\t{contig}\t{start}\t{end}\t{abs_thr}\n")
        print(f"\nWrote {args.out}")
    finally:
        if not args.keep_blast_db:
            shutil.rmtree(dbdir, ignore_errors=True)


if __name__ == "__main__":
    main()
