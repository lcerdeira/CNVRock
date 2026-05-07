"""
Add plasmid resistance gene variants to HS11286_extended.fasta using PLSDB.

Replaces manual per-gene NCBI searches (add_phase_c_genes.py, add_ctxm_variants.py)
with a systematic query against the PLSDB 2025 plasmid database, which has
72,360 plasmids pre-annotated with AMRFinderPlus 3.12.8 + RGI 6.0.3.

PLSDB API base: https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api/
Python package: conda install -c ccb-sb plsdbapi  (or pip install plsdbapi)

Strategy per gene variant
--------------------------
1. Query PLSDB filter_nuccore(AMR_genes=[variant]) → candidate plasmid accessions
2. Filter to Klebsiella host (TAXONOMY_genus = "Klebsiella") where possible
3. Prefer: complete topology > longest sequence > Klebsiella host
4. Download FASTA for the best representative
5. BLAST gene CDS onto the plasmid to get exact coordinates
6. Append contig to assets/HS11286_extended.fasta
7. Add row to assets/plasmid_refs/plasmid_gene_coords.tsv

Idempotent: genes already in coords TSV or contigs already in extended FASTA
are skipped.

Usage (HPC, with blast_env activated):
    python3 data/setup/add_plsdb_genes.py [--dry-run] [--genes blaCTX-M-65,blaCTX-M-14]

After running:
    sbatch hpc/build_extended_reference.sh
    N=$(wc -l < assets/kpsc_bam_accessions.txt)
    sbatch --array=1-${N}%50 hpc/remap_unmapped_to_plasmids.sh
    python3 data/setup/merge_plasmid_counts.py \\
        --counts-dir data/inputs/plasmid_remap_counts/ \\
        --store-path data/inputs/KpSC-plasmid-1000bp-npy/
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import requests

# ── Try Biopython for CDS lookup fallback ───────────────────────────────────
try:
    from Bio import Entrez, SeqIO
    Entrez.email = "louise.cerdeira@gmail.com"
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# ── Paths ────────────────────────────────────────────────────────────────────
REPO_DIR    = Path(__file__).resolve().parents[2]
ASSETS      = REPO_DIR / "assets"
PLASMID_DIR = ASSETS / "plasmid_refs"
EXTENDED_FA = ASSETS / "HS11286_extended.fasta"
COORDS_TSV  = PLASMID_DIR / "plasmid_gene_coords.tsv"

# ── PLSDB API ────────────────────────────────────────────────────────────────
PLSDB_BASE  = "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api"
PLSDB_NCBI  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SESSION      = requests.Session()
SESSION.headers.update({"Accept": "application/json", "User-Agent": "CNVRock/1.0"})

# ── BLAST thresholds ─────────────────────────────────────────────────────────
BLAST_MIN_IDENTITY = 80.0
BLAST_MIN_COVERAGE = 80.0

# ════════════════════════════════════════════════════════════════════════════
# Target gene variants
#
# Each entry:
#   gene             — internal name used in plasmid_gene_coords.tsv
#   amr_gene         — exact AMRFinderPlus gene symbol used in PLSDB query
#   family_pattern   — substring matched against BLAST subject (for CDS lookup)
#   absent_threshold — PCN below which gene is called absent
#   cds_acc          — optional: known NCBI CDS accession (avoids CDS search)
#   priority         — 1=highest; entries processed in priority order
# ════════════════════════════════════════════════════════════════════════════

GENE_TARGETS = [
    # ── CTX-M variants ── highest priority: these are the ST11/other FN root cause
    {
        "gene":             "blaCTX-M-65",
        "amr_gene":         "blaCTX-M-65",
        "family_pattern":   "CTX-M",
        "absent_threshold": 0.20,
        "priority":         1,
    },
    {
        "gene":             "blaCTX-M-27",
        "amr_gene":         "blaCTX-M-27",
        "family_pattern":   "CTX-M",
        "absent_threshold": 0.20,
        "priority":         1,
    },
    {
        "gene":             "blaCTX-M-3",
        "amr_gene":         "blaCTX-M-3",
        "family_pattern":   "CTX-M",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    {
        "gene":             "blaCTX-M-55",
        "amr_gene":         "blaCTX-M-55",
        "family_pattern":   "CTX-M",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    {
        "gene":             "blaCTX-M-2",
        "amr_gene":         "blaCTX-M-2",
        "family_pattern":   "CTX-M",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    # ── NDM variants ──
    {
        "gene":             "blaNDM-5",
        "amr_gene":         "blaNDM-5",
        "family_pattern":   "NDM",
        "absent_threshold": 0.20,
        "priority":         1,
    },
    {
        "gene":             "blaNDM-7",
        "amr_gene":         "blaNDM-7",
        "family_pattern":   "NDM",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    {
        "gene":             "blaNDM-4",
        "amr_gene":         "blaNDM-4",
        "family_pattern":   "NDM",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    # ── OXA-48-like variants ──
    # GT already covers the OXA-48 family, but reads from OXA-181/OXA-232
    # carriers currently map to the OXA-48 contig at low identity — adding
    # dedicated references increases sensitivity.
    {
        "gene":             "blaOXA-181",
        "amr_gene":         "blaOXA-181",
        "family_pattern":   "OXA",
        "absent_threshold": 0.20,
        "priority":         1,
    },
    {
        "gene":             "blaOXA-232",
        "amr_gene":         "blaOXA-232",
        "family_pattern":   "OXA",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    # ── KPC variants ──
    {
        "gene":             "blaKPC-3",
        "amr_gene":         "blaKPC-3",
        "family_pattern":   "KPC",
        "absent_threshold": 0.20,
        "priority":         2,
    },
    # ── Additional clinically important plasmid genes ──
    {
        "gene":             "mcr-1",
        "amr_gene":         "mcr-1",
        "family_pattern":   "mcr",
        "absent_threshold": 0.20,
        "priority":         3,
    },
    {
        "gene":             "blaSHV-12",
        "amr_gene":         "blaSHV-12",
        "family_pattern":   "SHV",
        "absent_threshold": 0.20,
        "priority":         3,
    },
]


# ════════════════════════════════════════════════════════════════════════════
# PLSDB API helpers
#
# Actual response format (confirmed from live API):
#   filter_nuccore?AMR_genes=X   → {"NUCCORE_ACC": ["acc1", "acc2", ...]}
#   filter_taxonomy?TAXONOMY_genus=X → {"NUCCORE_ACC": ["acc1", ...]}
#   summary/?NUCCORE_ACC=X       → {"Metadata_annotations": {"NUCCORE": {...}}}
#     NUCCORE fields: NUCCORE_ACC, NUCCORE_Topology, NUCCORE_Completeness,
#                     Length, NUCCORE_Description, AMR_genes
# ════════════════════════════════════════════════════════════════════════════

def _plsdb_filter_nuccore(amr_gene: str, timeout: int = 90) -> list[str]:
    """Return list of NUCCORE_ACC strings for plasmids annotated with amr_gene."""
    try:
        r = SESSION.get(f"{PLSDB_BASE}/filter_nuccore",
                        params={"AMR_genes": amr_gene}, timeout=timeout)
        r.raise_for_status()
        return r.json().get("NUCCORE_ACC", [])
    except Exception as e:
        print(f"  WARNING: PLSDB filter_nuccore failed for {amr_gene}: {e}",
              file=sys.stderr)
        return []


def _plsdb_filter_taxonomy(genus: str = "Klebsiella", timeout: int = 90) -> set[str]:
    """Return set of NUCCORE_ACC strings hosted in <genus>."""
    try:
        r = SESSION.get(f"{PLSDB_BASE}/filter_taxonomy",
                        params={"TAXONOMY_genus": genus}, timeout=timeout)
        r.raise_for_status()
        return set(r.json().get("NUCCORE_ACC", []))
    except Exception as e:
        print(f"  WARNING: PLSDB filter_taxonomy failed: {e}", file=sys.stderr)
        return set()


def _plsdb_summary(acc: str, timeout: int = 30) -> dict:
    """Fetch NUCCORE metadata dict for a single accession."""
    try:
        r = SESSION.get(f"{PLSDB_BASE}/summary/",
                        params={"NUCCORE_ACC": acc}, timeout=timeout)
        r.raise_for_status()
        return r.json().get("Metadata_annotations", {}).get("NUCCORE", {})
    except Exception:
        return {}


def _pick_best_representative(candidates: list[str],
                               klebsiella_accs: set[str],
                               max_inspect: int = 20) -> str | None:
    """
    From a list of candidate accessions, choose the best representative:
      1. Prefer Klebsiella host (intersection with taxonomy filter)
      2. Prefer complete + circular topology
      3. Prefer longer sequence
    Fetches summary metadata for up to max_inspect candidates to avoid
    excessive API calls.
    """
    if not candidates:
        return None

    # Prioritise Klebsiella hits; fall back to all candidates
    pool = [a for a in candidates if a in klebsiella_accs] or candidates
    # Prefer complete (CP/AP prefixes) over NZ_ WGS before even calling summary
    complete = [a for a in pool if not a.startswith("NZ_")]
    pool = complete[:max_inspect] if complete else pool[:max_inspect]

    best_acc, best_score = None, -1.0
    for acc in pool:
        meta  = _plsdb_summary(acc)
        time.sleep(0.2)
        desc  = meta.get("NUCCORE_Description", "")
        score = (
            2.0 if "Klebsiella" in desc else 0.0
        ) + (
            1.0 if meta.get("NUCCORE_Topology") == "circular" else 0.0
        ) + (
            1.0 if meta.get("NUCCORE_Completeness") == "complete" else 0.0
        ) + (
            meta.get("Length", 0) / 1_000_000
        )
        if score > best_score:
            best_score, best_acc = score, acc
    return best_acc


def _download_fasta_ncbi(acc: str, out_path: Path) -> bool:
    """
    Download a plasmid FASTA from NCBI using eutils REST API.
    Falls back to Biopython Entrez if requests fails.
    """
    ncbi_efetch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nucleotide&id={acc}&rettype=fasta&retmode=text"
        "&email=louise.cerdeira@gmail.com"
    )
    try:
        time.sleep(0.4)
        r = SESSION.get(ncbi_efetch, timeout=60)
        r.raise_for_status()
        fa = r.text
        if fa.startswith(">"):
            out_path.write_text(fa)
            return True
        print(f"  WARNING: NCBI returned unexpected content for {acc}",
              file=sys.stderr)
    except Exception as e:
        print(f"  WARNING: NCBI eutils request failed for {acc}: {e}",
              file=sys.stderr)

    # Biopython fallback
    if HAS_BIOPYTHON:
        try:
            time.sleep(0.4)
            h  = Entrez.efetch(db="nucleotide", id=acc,
                               rettype="fasta", retmode="text")
            fa = h.read(); h.close()
            if fa.startswith(">"):
                out_path.write_text(fa)
                return True
        except Exception as e:
            print(f"  WARNING: Biopython efetch fallback failed for {acc}: {e}",
                  file=sys.stderr)
    return False


# ════════════════════════════════════════════════════════════════════════════
# NCBI CDS lookup (Biopython fallback for gene coordinates)
# ════════════════════════════════════════════════════════════════════════════

def _ncbi_efetch(acc: str) -> str | None:
    """Fetch FASTA from NCBI eutils (no Biopython required)."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nucleotide&id={acc}&rettype=fasta&retmode=text"
        "&email=louise.cerdeira@gmail.com"
    )
    try:
        time.sleep(0.4)
        r = SESSION.get(url, timeout=60)
        r.raise_for_status()
        fa = r.text
        return fa if fa.startswith(">") else None
    except Exception:
        return None


def _ncbi_esearch(term: str) -> list[str]:
    """Search NCBI nucleotide and return a list of GI/accession IDs."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=nucleotide&term={requests.utils.quote(term)}&retmax=5&retmode=json"
        "&email=louise.cerdeira@gmail.com"
    )
    try:
        time.sleep(0.4)
        r = SESSION.get(url, timeout=30)
        r.raise_for_status()
        return r.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        return []


def _fetch_cds_fasta(gene: str, cds_acc: str | None = None) -> str | None:
    """
    Fetch a short CDS FASTA to use as BLAST query.
    Priority: explicit cds_acc > NCBI eutils search by gene name.
    Uses requests directly; Biopython is optional.
    """
    if cds_acc:
        fa = _ncbi_efetch(cds_acc)
        if fa:
            print(f"  CDS fetched via {cds_acc}")
            return fa
        print(f"  WARNING: direct CDS fetch failed ({cds_acc})", file=sys.stderr)

    # Fallback: Biopython Entrez search
    if HAS_BIOPYTHON:
        queries = [
            f'"{gene}"[gene] AND 500:1500[SLEN]',
            f'"{gene}"[title] AND beta-lactamase AND 500:1500[SLEN]',
        ]
        for q in queries:
            try:
                time.sleep(0.4)
                h = Entrez.esearch(db="nucleotide", term=q, retmax=5)
                rec = Entrez.read(h); h.close()
                ids = rec.get("IdList", [])
                if ids:
                    time.sleep(0.4)
                    h = Entrez.efetch(db="nucleotide", id=ids[0],
                                      rettype="fasta", retmode="text")
                    fa = h.read(); h.close()
                    if fa.startswith(">"):
                        print(f"  CDS fetched via NCBI search (query: {q[:60]}…)")
                        return fa
            except Exception as e:
                print(f"  WARNING: CDS search failed: {e}", file=sys.stderr)
        return None

    # Pure-requests fallback: eutils esearch
    queries = [
        f'"{gene}"[gene] AND 500:1500[SLEN]',
        f'"{gene}"[title] AND 500:1500[SLEN]',
    ]
    for q in queries:
        ids = _ncbi_esearch(q)
        if ids:
            fa = _ncbi_efetch(ids[0])
            if fa:
                print(f"  CDS fetched via eutils search (query: {q[:60]}…)")
                return fa
    return None


# ════════════════════════════════════════════════════════════════════════════
# BLAST helpers
# ════════════════════════════════════════════════════════════════════════════

def _blast_locate_gene(plasmid_fa: Path, cds_fasta: str) -> dict | None:
    """
    Run blastn of cds_fasta against plasmid_fa.
    Returns {"contig": ..., "start": int, "end": int} or None.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa",
                                     delete=False) as tf:
        tf.write(cds_fasta)
        query_path = tf.name
    try:
        result = subprocess.run(
            ["blastn",
             "-query",   query_path,
             "-subject", str(plasmid_fa),
             "-outfmt",  "6 sseqid pident length sstart send qlen",
             "-perc_identity", str(BLAST_MIN_IDENTITY),
             "-max_hsps", "1"],
            capture_output=True, text=True, check=True,
        )
        for line in result.stdout.strip().splitlines():
            contig, pct_id, hsp_len, s_start, s_end, q_len = line.split("\t")
            cov = float(hsp_len) / float(q_len) * 100
            if float(pct_id) >= BLAST_MIN_IDENTITY and cov >= BLAST_MIN_COVERAGE:
                s_lo = min(int(s_start), int(s_end))
                s_hi = max(int(s_start), int(s_end))
                print(f"  BLAST hit: {contig}:{s_lo}-{s_hi} "
                      f"(id={pct_id}%, cov={cov:.0f}%)")
                return {"contig": contig, "start": s_lo, "end": s_hi}
        print("  WARNING: no BLAST hit above thresholds", file=sys.stderr)
        return None
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  WARNING: blastn failed: {e}", file=sys.stderr)
        return None
    finally:
        os.unlink(query_path)


# ════════════════════════════════════════════════════════════════════════════
# Extended FASTA / coords TSV helpers
# ════════════════════════════════════════════════════════════════════════════

def _contig_in_extended(contig_id: str) -> bool:
    if not EXTENDED_FA.exists():
        return False
    with open(EXTENDED_FA) as f:
        for line in f:
            if line.startswith(">") and contig_id in line:
                return True
    return False


def _load_coords() -> list[dict]:
    if not COORDS_TSV.exists():
        return []
    with open(COORDS_TSV) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _save_coords(rows: list[dict]) -> None:
    fieldnames = ["gene", "contig", "start", "end", "absent_threshold"]
    with open(COORDS_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  Coords TSV updated → {len(rows)} genes total")


def _append_to_extended(plasmid_fa: Path) -> None:
    with open(EXTENDED_FA, "a") as out, open(plasmid_fa) as inp:
        for line in inp:
            out.write(line)
    print(f"  Appended {plasmid_fa.name} → {EXTENDED_FA.name}")


# ════════════════════════════════════════════════════════════════════════════
# Per-gene processing
# ════════════════════════════════════════════════════════════════════════════

def process_gene(gene_info: dict,
                 coords: list[dict],
                 klebsiella_accs: set[str],
                 dry_run: bool = False) -> bool:
    """
    Process one gene variant. Returns True if gene was added.
    """
    gene      = gene_info["gene"]
    amr_gene  = gene_info["amr_gene"]
    absent_t  = gene_info["absent_threshold"]
    cds_acc   = gene_info.get("cds_acc")

    existing_genes = {r["gene"] for r in coords}
    if gene in existing_genes:
        print(f"  {gene}: already in coords TSV — skip")
        return False

    # ── 1. Query PLSDB ────────────────────────────────────────────────────
    print(f"  Querying PLSDB for {amr_gene} …")
    records = _plsdb_filter_nuccore(amr_gene)
    if not records:
        print(f"  WARNING: no PLSDB records for {amr_gene}", file=sys.stderr)
        return False
    print(f"  PLSDB returned {len(records)} candidates")

    # ── 2. Pick best representative ───────────────────────────────────────
    best_acc = _pick_best_representative(records, klebsiella_accs)
    if not best_acc:
        print(f"  WARNING: could not identify a best accession", file=sys.stderr)
        return False
    print(f"  Best representative: {best_acc}")

    if dry_run:
        print(f"  [dry-run] would download {best_acc} and add {gene}")
        return False

    # ── 3. Download plasmid FASTA ─────────────────────────────────────────
    out_fa = PLASMID_DIR / f"{gene}.fasta"
    if out_fa.exists():
        print(f"  {out_fa.name} already downloaded")
    else:
        print(f"  Downloading FASTA for {best_acc} …")
        ok = _download_fasta_ncbi(best_acc, out_fa)
        if not ok or not out_fa.exists():
            print(f"  WARNING: FASTA download failed for {best_acc}",
                  file=sys.stderr)
            return False
        print(f"  Saved → {out_fa}")

    # Re-read accession from file header (PLSDB may return versioned acc)
    with open(out_fa) as f:
        first_line = f.readline()
    file_acc = first_line.split()[0].lstrip(">").strip()

    if _contig_in_extended(file_acc):
        print(f"  {file_acc} already in extended FASTA — skip append")
        return False

    # ── 4. Fetch CDS for BLAST ────────────────────────────────────────────
    print(f"  Fetching CDS sequence for {gene} …")
    cds_fasta = _fetch_cds_fasta(gene, cds_acc)
    if not cds_fasta:
        print(f"  WARNING: could not obtain CDS FASTA for {gene}",
              file=sys.stderr)
        return False

    # ── 5. BLAST to locate gene on plasmid ────────────────────────────────
    hit = _blast_locate_gene(out_fa, cds_fasta)
    if not hit:
        print(f"  ERROR: {gene} not found on {file_acc} by BLAST — "
              f"removing cached {out_fa.name} for retry", file=sys.stderr)
        out_fa.unlink(missing_ok=True)
        return False

    # ── 6. Append & update coords ─────────────────────────────────────────
    _append_to_extended(out_fa)
    coords.append({
        "gene":             gene,
        "contig":           hit["contig"],
        "start":            hit["start"],
        "end":              hit["end"],
        "absent_threshold": absent_t,
    })
    _save_coords(coords)
    print(f"  Added {gene} → {hit['contig']}:{hit['start']}-{hit['end']}")
    return True


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Add PLSDB-sourced plasmid resistance gene variants to "
                    "HS11286_extended.fasta and plasmid_gene_coords.tsv")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without downloading or modifying files")
    parser.add_argument("--genes", default="",
                        help="Comma-separated subset of gene names to process "
                             "(default: all targets in priority order)")
    parser.add_argument("--max-priority", type=int, default=3,
                        help="Only process genes with priority <= this value (default: 3)")
    args = parser.parse_args()

    if not EXTENDED_FA.exists():
        print(f"ERROR: {EXTENDED_FA} not found. Run get_plasmid_references.py first.",
              file=sys.stderr)
        sys.exit(1)

    PLASMID_DIR.mkdir(parents=True, exist_ok=True)
    coords = _load_coords()

    # ── Filter targets ────────────────────────────────────────────────────
    targets = sorted(GENE_TARGETS, key=lambda g: g["priority"])
    if args.genes:
        keep = set(args.genes.split(","))
        targets = [t for t in targets if t["gene"] in keep]
    targets = [t for t in targets if t["priority"] <= args.max_priority]

    if not targets:
        print("No targets to process.")
        return

    print(f"Targets to process: {len(targets)}")
    for t in targets:
        print(f"  priority={t['priority']}  {t['gene']}")

    # ── Pre-fetch Klebsiella accessions once (used for all genes) ─────────
    print("\nFetching Klebsiella plasmid accessions from PLSDB taxonomy filter …")
    klebsiella_accs = _plsdb_filter_taxonomy("Klebsiella")
    print(f"  {len(klebsiella_accs):,} Klebsiella-associated accessions in PLSDB")

    # ── Process each gene ─────────────────────────────────────────────────
    added, skipped, failed = [], [], []
    for gene_info in targets:
        gene = gene_info["gene"]
        print(f"\n{'='*60}")
        print(f"Processing {gene}  (priority {gene_info['priority']})")
        try:
            ok = process_gene(gene_info, coords, klebsiella_accs, dry_run=args.dry_run)
            (added if ok else skipped).append(gene)
        except Exception as e:
            print(f"  ERROR: unexpected failure for {gene}: {e}", file=sys.stderr)
            failed.append(gene)
        time.sleep(1)   # polite pause between API calls

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Added  : {len(added)}   {added}")
    print(f"  Skipped: {len(skipped)}  (already present or no PLSDB hit)")
    print(f"  Failed : {len(failed)}  {failed}")

    if added and not args.dry_run:
        print(f"""
Next steps on HPC:
  1. Re-index extended reference:
       sbatch hpc/build_extended_reference.sh

  2. Remap unmapped reads (all samples):
       N=$(wc -l < assets/kpsc_bam_accessions.txt)
       sbatch --array=1-${{N}}%50 hpc/remap_unmapped_to_plasmids.sh

  3. Merge new plasmid counts:
       python3 data/setup/merge_plasmid_counts.py \\
           --counts-dir data/inputs/plasmid_remap_counts/ \\
           --store-path data/inputs/KpSC-plasmid-1000bp-npy/

  4. Create new experiment config (copy previous, update out_dir):
       cp -r models/experiments/30 models/experiments/31
       # edit models/experiments/31/config.yaml
""")


if __name__ == "__main__":
    main()
