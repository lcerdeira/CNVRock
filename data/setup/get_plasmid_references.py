"""
Download representative plasmid sequences for blaKPC-2, blaCTX-M-15, blaNDM-1
and build an extended reference FASTA for mapping.

Strategy
--------
blaKPC-2    — NC_016846.1 is already present in assets/HS11286.fasta (confirmed by
              AMRFinder: positions 20557-21435 on NC_016846.1). Extract it from the
              existing FASTA; no download needed.

blaCTX-M-15 — Fetch a representative K. pneumoniae plasmid from NCBI RefSeq via
blaNDM-1      Entrez esearch. The script picks the first RefSeq match, downloads it,
              runs AMRFinder to confirm the gene is present and record coordinates,
              then saves the sequence as assets/plasmid_refs/{gene}.fasta.

Outputs
-------
assets/plasmid_refs/
    blaKPC-2.fasta          — extracted NC_016846.1 from HS11286.fasta
    blaCTX-M-15.fasta       — downloaded plasmid
    blaNDM-1.fasta          — downloaded plasmid
    plasmid_gene_coords.tsv — contig, start, end per gene (for interval list and CNV caller)

assets/HS11286_extended.fasta
    Chromosome (NC_016845.1) + plasmid contigs concatenated.  Use this as the
    mapping reference for all new BAMs.

Usage (run from CNVRock/ repo root, requires biopython and amrfinder in PATH):
    python3 data/setup/get_plasmid_references.py
"""

import csv
import os
import subprocess
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

REPO_DIR   = Path("/home/lshlt19/CNVRock")
ASSETS     = REPO_DIR / "assets"
PLASMID_DIR = ASSETS / "plasmid_refs"

HS11286_FASTA    = ASSETS / "HS11286.fasta"
EXTENDED_FASTA   = ASSETS / "HS11286_extended.fasta"
COORDS_TSV       = PLASMID_DIR / "plasmid_gene_coords.tsv"

# AMRFinder gene patterns used to confirm presence and extract coordinates.
# Values are substrings matched against the "Gene symbol" column.
GENE_PATTERNS = {
    "blaKPC-2":    ["blaKPC"],
    "blaCTX-M-15": ["blaCTX-M"],
    "blaNDM-1":    ["blaNDM"],
}

# NCBI search query per gene: retrieves RefSeq plasmids in K. pneumoniae.
# The script picks the first hit; adjust the query if needed.
NCBI_QUERIES = {
    "blaCTX-M-15": (
        'blaCTX-M-15 AND "Klebsiella pneumoniae"[organism] AND plasmid[filter]'
    ),
    "blaNDM-1": (
        'blaNDM-1 AND "Klebsiella pneumoniae"[organism] AND plasmid[filter]'
    ),
}

Entrez.email = "louise.cerdeira@gmail.com"


# ---------------------------------------------------------------------------
# Step 1: extract NC_016846.1 (blaKPC-2) from existing HS11286.fasta
# ---------------------------------------------------------------------------

def extract_kpc_plasmid() -> dict:
    """Extract NC_016846.1 from HS11286.fasta and write to plasmid_refs/."""
    out_fa = PLASMID_DIR / "blaKPC-2.fasta"
    if out_fa.exists():
        print("blaKPC-2.fasta already exists — skipping extraction.")
    else:
        found = False
        with open(out_fa, "w") as fh:
            for rec in SeqIO.parse(str(HS11286_FASTA), "fasta"):
                if rec.id == "NC_016846.1":
                    SeqIO.write(rec, fh, "fasta")
                    found = True
                    print(f"Extracted NC_016846.1 ({len(rec.seq):,} bp) → {out_fa}")
                    break
        if not found:
            out_fa.unlink(missing_ok=True)
            raise RuntimeError(
                "NC_016846.1 not found in HS11286.fasta. "
                "Ensure HS11286.fasta is the full genome (chromosome + plasmids)."
            )

    # Known AMRFinder coordinates (confirmed from reference genome AMRFinder run)
    return {
        "gene":   "blaKPC-2",
        "contig": "NC_016846.1",
        "start":  20557,
        "end":    21435,
    }


# ---------------------------------------------------------------------------
# Step 2: download blaCTX-M-15 and blaNDM-1 plasmids from NCBI
# ---------------------------------------------------------------------------

def _entrez_fetch_plasmid(gene: str, query: str) -> tuple[str, str]:
    """Search NCBI for a plasmid, download the first RefSeq hit.

    Returns (accession, sequence_string).
    """
    print(f"\nSearching NCBI for {gene}:  {query}")
    time.sleep(0.5)  # be polite to NCBI

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=20, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        raise RuntimeError(f"No NCBI hits for query: {query}")

    # Prefer longer sequences (more likely to be complete plasmids).
    # Fetch summaries to check sequence lengths before downloading.
    time.sleep(0.5)
    handle   = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    best_id = max(summaries, key=lambda s: int(s.get("Length", 0)))["Id"]

    print(f"  Found {len(ids)} hit(s). Downloading longest match (ID {best_id}) …")
    time.sleep(0.5)

    handle = Entrez.efetch(db="nucleotide", id=best_id, rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()

    lines  = fasta_text.strip().split("\n")
    header = lines[0]
    acc    = header.split()[0].lstrip(">")
    seq    = "".join(lines[1:])

    return acc, fasta_text


def download_plasmid(gene: str) -> None:
    """Download and save a representative plasmid sequence for gene."""
    out_fa = PLASMID_DIR / f"{gene}.fasta"
    if out_fa.exists():
        print(f"{gene}.fasta already exists — skipping download.")
        return

    query       = NCBI_QUERIES[gene]
    acc, fasta  = _entrez_fetch_plasmid(gene, query)
    out_fa.write_text(fasta)
    size_kb = len(fasta) // 1000
    print(f"  Saved {acc} ({size_kb} kb) → {out_fa}")


# ---------------------------------------------------------------------------
# Step 3: run AMRFinder to confirm gene presence and get coordinates
# ---------------------------------------------------------------------------

def amrfinder_coords(gene: str, fasta_path: Path) -> dict:
    """Run AMRFinder on a plasmid FASTA and return gene coordinates."""
    tsv_path = PLASMID_DIR / f"{gene}.amrfinder.tsv"
    if not tsv_path.exists():
        print(f"  Running AMRFinder on {fasta_path.name} …", flush=True)
        result = subprocess.run(
            [
                "amrfinder",
                "-O", "Klebsiella_pneumoniae",
                "-n", str(fasta_path),
                "--output", str(tsv_path),
                "--quiet",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"AMRFinder failed for {gene}:\n{result.stderr[:400]}")

    patterns = GENE_PATTERNS[gene]
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("Gene symbol", "")
            if any(p in sym for p in patterns):
                contig = row.get("Contig id", "")
                start  = int(row.get("Start", 0))
                end    = int(row.get("Stop",  0))
                print(f"  AMRFinder: {sym} on {contig}:{start}-{end}")
                return {"gene": gene, "contig": contig, "start": start, "end": end}

    raise RuntimeError(
        f"AMRFinder did not find {gene} in {fasta_path}. "
        f"Try a different reference plasmid (edit NCBI_QUERIES in this script)."
    )


# ---------------------------------------------------------------------------
# Step 4: build extended reference FASTA
# ---------------------------------------------------------------------------

def build_extended_fasta(plasmid_fastas: list[Path]) -> None:
    """Concatenate chromosome + plasmid sequences into HS11286_extended.fasta."""
    if EXTENDED_FASTA.exists():
        print(f"\n{EXTENDED_FASTA.name} already exists — skipping.")
        return

    print(f"\nBuilding {EXTENDED_FASTA.name} …")
    contigs_written = []
    with open(EXTENDED_FASTA, "w") as out:
        # 1. Chromosome only from HS11286.fasta
        for rec in SeqIO.parse(str(HS11286_FASTA), "fasta"):
            if rec.id == "NC_016845.1":
                SeqIO.write(rec, out, "fasta")
                contigs_written.append(rec.id)
                break

        # 2. blaKPC-2 plasmid (NC_016846.1 already in HS11286 — keep same ID)
        kpc_fa = PLASMID_DIR / "blaKPC-2.fasta"
        for rec in SeqIO.parse(str(kpc_fa), "fasta"):
            SeqIO.write(rec, out, "fasta")
            contigs_written.append(rec.id)

        # 3. Downloaded plasmid sequences (CTX-M-15, NDM-1)
        for fa in plasmid_fastas:
            for rec in SeqIO.parse(str(fa), "fasta"):
                SeqIO.write(rec, out, "fasta")
                contigs_written.append(rec.id)
                break  # one contig per file

    print(f"Extended reference contigs: {contigs_written}")
    print(f"Saved → {EXTENDED_FASTA}")


# ---------------------------------------------------------------------------
# Step 5: write plasmid_gene_coords.tsv
# ---------------------------------------------------------------------------

def write_coords(coord_rows: list[dict]) -> None:
    """Save gene coordinates for use by the interval list builder and CNV caller."""
    with open(COORDS_TSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["gene", "contig", "start", "end"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(coord_rows)
    print(f"\nSaved plasmid gene coords ({len(coord_rows)} genes) → {COORDS_TSV}")
    for r in coord_rows:
        print(f"  {r['gene']:15s} {r['contig']}:{r['start']}-{r['end']}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    PLASMID_DIR.mkdir(parents=True, exist_ok=True)

    coord_rows = []

    # blaKPC-2 — already in HS11286.fasta
    print("=== blaKPC-2 (NC_016846.1) ===")
    coord_rows.append(extract_kpc_plasmid())

    # blaCTX-M-15 and blaNDM-1 — download from NCBI
    downloaded_fastas = []
    for gene in ["blaCTX-M-15", "blaNDM-1"]:
        print(f"\n=== {gene} ===")
        download_plasmid(gene)
        fa = PLASMID_DIR / f"{gene}.fasta"
        coords = amrfinder_coords(gene, fa)
        coord_rows.append(coords)
        downloaded_fastas.append(fa)

    write_coords(coord_rows)
    build_extended_fasta(downloaded_fastas)

    print("\nDone. Next steps:")
    print("  1. Index the extended reference:   hpc/build_extended_reference.sh")
    print("  2. Collect plasmid read counts:     hpc/collect_plasmid_readcounts.sh")
    print("  3. Build plasmid NPY store:         python3 models/experiments/21/plasmid_readcounts_to_npy.py")
    print("  4. For new BAMs: update REFERENCE in hpc/download_sra.sh to HS11286_extended.fasta")


if __name__ == "__main__":
    main()
