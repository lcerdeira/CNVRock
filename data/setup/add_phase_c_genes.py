"""
Add Phase C plasmid resistance genes to HS11286_extended.fasta.

Genes added:
  blaTEM-1       — broad-spectrum beta-lactamase (TEM-1)
  blaOXA-48      — OXA-48 carbapenemase
  qnrB1          — plasmid-mediated quinolone resistance
  aac6-Ib-cr     — aminoglycoside acetyltransferase with ciprofloxacin activity

Strategy
--------
For each gene:
  1. Download a representative K. pneumoniae plasmid from NCBI RefSeq.
  2. Locate the gene on the plasmid via blastn.
  3. Append the plasmid contig to assets/HS11286_extended.fasta.
  4. Add a row to assets/plasmid_refs/plasmid_gene_coords.tsv.

Usage (requires blast_env on PATH for blastn):
    python3 data/setup/add_phase_c_genes.py

After running:
    sbatch hpc/build_extended_reference.sh
    N=$(wc -l < assets/kpsc_bam_accessions.txt)
    sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh
    python3 data/setup/merge_plasmid_counts.py \\
        --counts-dir data/inputs/plasmid_remap_counts/ \\
        --store-path data/inputs/KpSC-plasmid-1000bp-npy/
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

NEW_GENES = [
    {
        "gene":             "blaTEM-1",
        "query":            ('blaTEM-1 AND "Klebsiella pneumoniae"[organism] '
                             'AND plasmid[filter] AND refseq[filter]'),
        "pattern":          "blaTEM",
        "absent_threshold": "0.20",   # revisit if chromosomal TEM cross-mapping causes FPs
        "slen_range":       "800:1100",  # TEM CDS ~861 bp
    },
    {
        "gene":             "blaOXA-48",
        "query":            ('blaOXA-48 AND "Klebsiella pneumoniae"[organism] '
                             'AND plasmid[filter] AND 5000:500000[SLEN]'),
        "pattern":          "blaOXA",
        "absent_threshold": "0.20",
        "slen_range":       "700:950",   # OXA-48 CDS ~828 bp
    },
    {
        "gene":             "qnrB1",
        "query":            ('qnrB1 AND "Klebsiella pneumoniae"[organism] '
                             'AND plasmid[filter] AND refseq[filter]'),
        "pattern":          "qnrB",
        "absent_threshold": "0.20",
        "slen_range":       "580:750",   # qnrB CDS ~657 bp
    },
    {
        "gene":             "aac6-Ib-cr",
        "query":            ('"Klebsiella pneumoniae"[organism] AND aac AND Ib-cr '
                             'AND plasmid[filter] AND 10000:300000[SLEN]'),
        "pattern":          "aac",
        "absent_threshold": "0.20",
        "slen_range":       "450:650",   # aac(6')-Ib CDS ~543 bp
        # Internal name "aac6-Ib-cr" doesn't match NCBI gene name aac(6')-Ib-cr;
        # search without quotes so the parentheses aren't treated as field tags.
        "cds_query":        'aac AND Ib-cr AND aminoglycoside AND 450:650[SLEN]',
    },
]

BLAST_MIN_IDENTITY = 80.0
BLAST_MIN_COVERAGE = 80.0


# ---------------------------------------------------------------------------
# Helpers (identical to add_ctxm_variants.py)
# ---------------------------------------------------------------------------

def _fetch_plasmid_ncbi(gene: str, query: str) -> tuple[str, str]:
    print(f"Searching NCBI for {gene}: {query}")
    time.sleep(0.5)

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=20, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        raise RuntimeError(f"No NCBI hits for: {query}")

    time.sleep(0.5)
    handle    = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    # Prefer plasmid-sized sequences (<500 kb) to avoid picking chromosome contigs.
    plasmid_sized = [s for s in summaries if int(s.get("Length", 0)) < 500_000]
    candidates = plasmid_sized if plasmid_sized else summaries
    # NZ_ accessions are WGS assembly contigs, not standalone plasmids.
    # Prefer CP/AP/classical accessions (complete plasmid sequences).
    complete = [s for s in candidates
                if not s.get("AccessionVersion", "").startswith("NZ_")]
    candidates = complete if complete else candidates
    best_id = max(candidates, key=lambda s: int(s.get("Length", 0)))["Id"]

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
    query_seq = ""
    with open(query_gene_fa) as f:
        for line in f:
            if not line.startswith(">"):
                query_seq += line.strip()
    if not query_seq:
        return None

    try:
        result = subprocess.run(
            ["blastn", "-query", query_gene_fa, "-subject", plasmid_fasta,
             "-outfmt", "6 sseqid pident length sstart send qlen",
             "-perc_identity", str(BLAST_MIN_IDENTITY),
             "-max_hsps", "1"],
            capture_output=True, text=True, check=True,
        )
        for line in result.stdout.strip().splitlines():
            parts = line.split("\t")
            contig, pct_id, hsp_len, s_start, s_end, q_len = parts
            coverage = float(hsp_len) / float(q_len) * 100
            if float(pct_id) >= BLAST_MIN_IDENTITY and coverage >= BLAST_MIN_COVERAGE:
                s_lo = min(int(s_start), int(s_end))
                s_hi = max(int(s_start), int(s_end))
                print(f"  BLAST: {gene_name} found at {contig}:{s_lo}-{s_hi} "
                      f"(identity={pct_id}%, coverage={coverage:.0f}%)")
                return {"contig": contig, "start": s_lo, "end": s_hi}
        return None
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  WARNING: blastn failed: {e}", file=sys.stderr)
        return None


def _fetch_gene_sequence(gene_name: str, slen_range: str, cds_query: str | None = None) -> str:
    if cds_query:
        query = cds_query
    elif gene_name.startswith("bla"):
        query = f'"{gene_name}"[title] AND beta-lactamase[title] AND {slen_range}[SLEN]'
    else:
        query = f'"{gene_name}"[title] AND {slen_range}[SLEN]'

    print(f"  Fetching reference gene sequence for {gene_name} …")
    time.sleep(0.5)

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=5)
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        raise RuntimeError(
            f"Cannot find reference CDS sequence for {gene_name}. "
            f"Provide a fallback accession in FALLBACK_ACCS."
        )

    time.sleep(0.5)
    handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
    fasta  = handle.read()
    handle.close()
    return fasta


def _contig_in_extended(contig_id: str) -> bool:
    with open(EXTENDED_FA) as f:
        for line in f:
            if line.startswith(">") and contig_id in line:
                return True
    return False


def _load_coords() -> list[dict]:
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
        gene       = gene_info["gene"]
        query      = gene_info["query"]
        absent     = gene_info["absent_threshold"]
        slen_range = gene_info["slen_range"]

        print(f"\n{'='*60}")
        print(f"Processing {gene}")

        if gene in existing_genes:
            print(f"  {gene} already in coords TSV — skipping.")
            continue

        out_fa = PLASMID_DIR / f"{gene}.fasta"
        if out_fa.exists():
            print(f"  {gene}.fasta already downloaded.")
            with open(out_fa) as f:
                acc = f.readline().split()[0].lstrip(">")
        else:
            acc, fasta_text = _fetch_plasmid_ncbi(gene, query)
            out_fa.write_text(fasta_text)
            print(f"  Saved → {out_fa}")

        if _contig_in_extended(acc):
            print(f"  {acc} already in HS11286_extended.fasta — skipping append.")
        else:
            print(f"  Locating {gene} on {acc} via BLAST …")
            gene_fasta_text = _fetch_gene_sequence(gene, slen_range, gene_info.get("cds_query"))
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tf:
                tf.write(gene_fasta_text)
                gene_query_path = tf.name
            try:
                coords_hit = _blast_find_gene(gene, str(out_fa), gene_query_path)
            finally:
                os.unlink(gene_query_path)

            if coords_hit is None:
                print(f"  ERROR: {gene} not found on {acc} by BLAST. "
                      f"Deleting cached {out_fa.name} — re-run to try a different plasmid.",
                      file=sys.stderr)
                out_fa.unlink(missing_ok=True)
                continue

            print(f"  Appending {acc} to {EXTENDED_FA.name} …")
            with open(EXTENDED_FA, "a") as out, open(out_fa) as inp:
                for line in inp:
                    out.write(line)
            print(f"  Appended.")

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
    print("  2. Remap unmapped reads (all samples, all genes):")
    print("       N=$(wc -l < assets/kpsc_bam_accessions.txt)")
    print("       sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh")
    print("  3. Merge new counts:")
    print("       python3 data/setup/merge_plasmid_counts.py \\")
    print("           --counts-dir data/inputs/plasmid_remap_counts/ \\")
    print("           --store-path data/inputs/KpSC-plasmid-1000bp-npy/")
    print("  4. Create models/experiments/27/config.yaml and run exp 27")


if __name__ == "__main__":
    main()
