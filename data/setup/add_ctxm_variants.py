"""
Add blaCTX-M-14 (CTX-M-9 group) representative plasmid to HS11286_extended.fasta.

Context
-------
Exp 25 showed blaCTX-M-15 reference (MK552109.1) misses ~33 samples carrying
non-CTX-M-15 variants (FNR=0.27, FN PCN=0.00). The CTX-M-9 group (CTX-M-14,
CTX-M-27, CTX-M-9) accounts for most non-CTX-M-15 CTX-M carriers in KpSC.

Strategy
--------
1. Download a representative blaCTX-M-14 plasmid from NCBI.
2. Confirm the gene is present using blastn (blast_env).
3. Append the plasmid contig to HS11286_extended.fasta.
4. Update assets/plasmid_refs/plasmid_gene_coords.tsv with the new locus.
5. Re-index the extended reference and re-run remap_unmapped_to_plasmids.sh.

Usage (requires blast_env on PATH for blastn):
    python3 data/setup/add_ctxm_variants.py

After running:
    sbatch hpc/build_extended_reference.sh      # re-index BWA
    sbatch --array=... hpc/remap_unmapped_to_plasmids.sh  # remap with new contigs
    python3 data/setup/merge_plasmid_counts.py ...         # merge new counts
"""

import csv
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

REPO_DIR    = Path("/home/lshlt19/CNVRock")
ASSETS      = REPO_DIR / "assets"
PLASMID_DIR = ASSETS / "plasmid_refs"
EXTENDED_FA = ASSETS / "HS11286_extended.fasta"
COORDS_TSV  = PLASMID_DIR / "plasmid_gene_coords.tsv"

Entrez.email = "louise.cerdeira@gmail.com"

# New variants to add — NCBI search queries
# CTX-M-14 is representative of the CTX-M-9 group (CTX-M-9, -14, -27 are ~97% identical)
NEW_GENES = [
    {
        "gene":    "blaCTX-M-14",
        "query":   ('blaCTX-M-14 AND "Klebsiella pneumoniae"[organism] '
                    'AND plasmid[filter] AND refseq[filter]'),
        "pattern": "blaCTX-M",    # AMRFinder/BLAST gene symbol substring
        "absent_threshold": "1.00",  # same as CTX-M-15 (cross-mapping concern)
    },
]

BLAST_MIN_IDENTITY = 80.0
BLAST_MIN_COVERAGE = 80.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fetch_plasmid_ncbi(gene: str, query: str) -> tuple[str, str]:
    """Search NCBI nucleotide, download the longest RefSeq hit as FASTA.

    Returns (accession, fasta_text).
    """
    print(f"Searching NCBI for {gene}: {query}")
    time.sleep(0.5)

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=20, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        raise RuntimeError(f"No NCBI hits for: {query}")

    time.sleep(0.5)
    handle     = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries  = Entrez.read(handle)
    handle.close()
    best_id    = max(summaries, key=lambda s: int(s.get("Length", 0)))["Id"]

    print(f"  Downloading best match (id {best_id}, {len(ids)} candidates) …")
    time.sleep(0.5)

    handle     = Entrez.efetch(db="nucleotide", id=best_id, rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()

    lines = fasta_text.strip().split("\n")
    acc   = lines[0].split()[0].lstrip(">")
    print(f"  Downloaded: {acc} ({len(''.join(lines[1:])):,} bp)")
    return acc, fasta_text


def _blast_find_gene(gene_name: str, plasmid_fasta: str, query_gene_fa: str) -> dict | None:
    """Use blastn to locate the gene on the plasmid. Returns {contig, start, end} or None."""
    query_seq = ""
    with open(query_gene_fa) as f:
        for line in f:
            if not line.startswith(">"):
                query_seq += line.strip()
    query_len = len(query_seq)
    if query_len == 0:
        return None

    try:
        result = subprocess.run(
            ["blastn", "-query", query_gene_fa, "-subject", plasmid_fasta,
             "-outfmt", "6 sseqid pident length sstart send qlen",
             "-perc_identity", str(BLAST_MIN_IDENTITY),
             "-max_hsps", "1"],
            capture_output=True, text=True, check=True,
        )
        best = None
        for line in result.stdout.strip().splitlines():
            parts = line.split("\t")
            contig, pct_id, hsp_len, s_start, s_end, q_len = parts
            coverage = float(hsp_len) / float(q_len) * 100
            if float(pct_id) >= BLAST_MIN_IDENTITY and coverage >= BLAST_MIN_COVERAGE:
                s_lo = min(int(s_start), int(s_end))
                s_hi = max(int(s_start), int(s_end))
                best = {"contig": contig, "start": s_lo, "end": s_hi}
                print(f"  BLAST: {gene_name} found at {contig}:{s_lo}-{s_hi} "
                      f"(identity={pct_id}%, coverage={coverage:.0f}%)")
                break
        return best
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  WARNING: blastn failed: {e}", file=sys.stderr)
        return None


def _fetch_gene_sequence(gene_name: str) -> str:
    """Fetch the reference gene CDS from NCBI as a FASTA string (for BLAST query)."""
    # Search for the gene CDS sequence
    query = f'{gene_name}[gene] AND "Klebsiella pneumoniae"[organism] AND CDS[feature key]'
    print(f"  Fetching reference gene sequence for {gene_name} …")
    time.sleep(0.5)

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=5)
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        # Fallback: use a known NCBI accession for the gene
        fallback = {"blaCTX-M-14": "AY077516.1"}
        acc = fallback.get(gene_name)
        if acc is None:
            raise RuntimeError(f"Cannot find reference sequence for {gene_name}")
        ids = [acc]

    time.sleep(0.5)
    handle   = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
    fasta    = handle.read()
    handle.close()
    return fasta


def _contig_in_extended(contig_id: str) -> bool:
    """Check whether contig_id is already in HS11286_extended.fasta."""
    with open(EXTENDED_FA) as f:
        for line in f:
            if line.startswith(">") and contig_id in line:
                return True
    return False


def _load_coords() -> list[dict]:
    """Load existing plasmid_gene_coords.tsv rows."""
    rows = []
    with open(COORDS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def _save_coords(rows: list[dict]) -> None:
    fieldnames = ["gene", "contig", "start", "end", "absent_threshold"]
    with open(COORDS_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Updated {COORDS_TSV} ({len(rows)} genes)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    PLASMID_DIR.mkdir(parents=True, exist_ok=True)

    if not EXTENDED_FA.exists():
        print(f"ERROR: {EXTENDED_FA} not found. Run get_plasmid_references.py first.",
              file=sys.stderr)
        sys.exit(1)

    coords = _load_coords()
    existing_genes = {r["gene"] for r in coords}

    for gene_info in NEW_GENES:
        gene   = gene_info["gene"]
        query  = gene_info["query"]
        absent = gene_info["absent_threshold"]

        print(f"\n{'='*60}")
        print(f"Processing {gene}")

        if gene in existing_genes:
            print(f"  {gene} already in coords TSV — skipping.")
            continue

        # ── Download plasmid ──────────────────────────────────────────────────
        out_fa = PLASMID_DIR / f"{gene}.fasta"
        if out_fa.exists():
            print(f"  {gene}.fasta already downloaded.")
            # Read the accession from the existing file
            with open(out_fa) as f:
                acc = f.readline().split()[0].lstrip(">")
        else:
            acc, fasta_text = _fetch_plasmid_ncbi(gene, query)
            out_fa.write_text(fasta_text)
            print(f"  Saved → {out_fa}")

        # Check if contig already in extended reference
        if _contig_in_extended(acc):
            print(f"  {acc} already in HS11286_extended.fasta — skipping append.")
        else:
            # ── Find gene on plasmid via BLAST ────────────────────────────────
            print(f"  Locating {gene} on {acc} via BLAST …")
            gene_fasta_text = _fetch_gene_sequence(gene)
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                             delete=False) as tf:
                tf.write(gene_fasta_text)
                gene_query_path = tf.name
            try:
                coords_hit = _blast_find_gene(gene, str(out_fa), gene_query_path)
            finally:
                os.unlink(gene_query_path)

            if coords_hit is None:
                print(f"  ERROR: {gene} not found on {acc} by BLAST. "
                      f"Try a different plasmid.", file=sys.stderr)
                continue

            # ── Append plasmid contig to extended FASTA ───────────────────────
            print(f"  Appending {acc} to {EXTENDED_FA.name} …")
            with open(EXTENDED_FA, "a") as out, open(out_fa) as inp:
                for line in inp:
                    out.write(line)
            print(f"  Appended.")

            # ── Update coords ─────────────────────────────────────────────────
            coords.append({
                "gene":             gene,
                "contig":           coords_hit["contig"],
                "start":            coords_hit["start"],
                "end":              coords_hit["end"],
                "absent_threshold": absent,
            })
            _save_coords(coords)
            print(f"  Added {gene} → {coords_hit['contig']}:{coords_hit['start']}-{coords_hit['end']}")

    print(f"\n{'='*60}")
    print("Done. Next steps:")
    print("  1. Re-index extended reference:")
    print("       sbatch hpc/build_extended_reference.sh")
    print("  2. Re-run unmapped read remapping (new contigs will attract CTX-M-14 reads):")
    print("       N=$(wc -l < assets/kpsc_bam_accessions.txt)")
    print("       sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh")
    print("  3. Merge new counts:")
    print("       python3 data/setup/merge_plasmid_counts.py \\")
    print("           --counts-dir data/inputs/plasmid_remap_counts/ \\")
    print("           --store-path data/inputs/KpSC-plasmid-1000bp-npy/")
    print("  4. Run exp 26 wrap_up")


if __name__ == "__main__":
    main()
