#!/usr/bin/env python3
"""
Re-annotate plasmid_gene_coords.tsv from GenBank annotations of each contig in
HS11286_extended.fasta.

Complements reannotate_plasmid_genes_blast.py: GenBank annotations are
curated by submitters and explicitly say "this CDS is blaKPC-2 from position
X to Y on accession Z". BLAST can drift to paralogous hits; GenBank is the
authoritative source where it exists.

Outputs the union/disagreement between the two approaches so we can pick the
right coordinates manually if they disagree.

Workflow
--------
1. For each contig in HS11286_extended.fasta (from .fai or --contigs CSV),
   fetch its GenBank record via NCBI Entrez (cached to
   assets/plasmid_refs/genbank_records/).
2. Parse the GenBank record's CDS / gene features.
3. Match each feature's `gene` / `product` / `note` qualifiers against the
   AMR gene panel (substring + alias match).
4. Write the per-contig hits to a TSV.
5. If --merge-blast is given, also reads the BLAST output TSV and writes a
   `_combined` TSV that flags agreement / disagreement.

Usage
-----
    # Fetch + parse GenBank for every contig in the extended reference:
    python3 data/setup/reannotate_plasmid_genes_genbank.py \\
        --reference assets/HS11286_extended.fasta \\
        --out       assets/plasmid_refs/plasmid_gene_coords_genbank.tsv \\
        --fetch-genbank

    # Compare with BLAST result:
    python3 data/setup/reannotate_plasmid_genes_genbank.py \\
        --reference     assets/HS11286_extended.fasta \\
        --out           assets/plasmid_refs/plasmid_gene_coords_genbank.tsv \\
        --merge-blast   assets/plasmid_refs/plasmid_gene_coords_blast.tsv \\
        --combined-out  assets/plasmid_refs/plasmid_gene_coords_combined.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# ── Gene panel (with alias regexes for matching messy GenBank annotations) ───
# Each entry: (canonical_name, regex_pattern, default_absent_threshold)
GENE_PANEL = [
    ("blaKPC-2",    r"(bla[_]?KPC[-_]?2|kpc[-_]?2|KPC-?2)",          0.20),
    ("blaCTX-M-15", r"(bla[_]?CTX[-_]?M[-_]?15|CTX-?M-?15)",          0.50),
    ("blaCTX-M-14", r"(bla[_]?CTX[-_]?M[-_]?14|CTX-?M-?14)",          1.00),
    ("blaCTX-M-65", r"(bla[_]?CTX[-_]?M[-_]?65|CTX-?M-?65)",          0.20),
    ("blaCTX-M-27", r"(bla[_]?CTX[-_]?M[-_]?27|CTX-?M-?27)",          0.20),
    ("blaNDM-1",    r"(bla[_]?NDM[-_]?1|NDM-?1)\b",                   0.20),
    ("blaNDM-5",    r"(bla[_]?NDM[-_]?5|NDM-?5)",                     0.20),
    ("blaOXA-48",   r"(bla[_]?OXA[-_]?48|OXA-?48)",                   0.20),
    ("blaOXA-181",  r"(bla[_]?OXA[-_]?181|OXA-?181)",                 0.20),
    ("blaTEM-1",    r"(bla[_]?TEM[-_]?1|TEM-?1)",                     0.20),
    ("qnrB1",       r"(qnrB[-_]?1)\b",                                0.20),
    ("aac6-Ib-cr",  r"aac\(?6.?\)?[-_]?Ib[-_]?cr",                    0.10),
]

NCBI_EMAIL = "louise.cerdeira@gmail.com"


def contigs_from_fai(reference: Path) -> list[str]:
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        sys.exit(f"ERROR: missing {fai}. Run `samtools faidx {reference}` first.")
    out = []
    with open(fai) as fh:
        for line in fh:
            out.append(line.split("\t")[0])
    return out


def fetch_genbank(accession: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{accession}.gb"
    if path.exists() and path.stat().st_size > 0:
        return path
    try:
        from Bio import Entrez
    except ImportError:
        sys.exit("ERROR: pip install biopython")
    Entrez.email = NCBI_EMAIL
    print(f"  Fetching {accession} from NCBI …", flush=True)
    handle = Entrez.efetch(db="nucleotide", id=accession,
                           rettype="gb", retmode="text")
    text = handle.read()
    handle.close()
    path.write_text(text)
    return path


def parse_features(gb_path: Path) -> list[dict]:
    """Return [{type, gene, product, note, start, end, strand}, …] for
    every CDS / gene feature in the GenBank record."""
    try:
        from Bio import SeqIO
    except ImportError:
        sys.exit("ERROR: pip install biopython")
    rec = SeqIO.read(gb_path, "genbank")
    rows = []
    for feat in rec.features:
        if feat.type not in {"CDS", "gene", "misc_feature"}:
            continue
        q = feat.qualifiers
        rows.append({
            "type":    feat.type,
            "gene":    (q.get("gene", [""])  [0] or ""),
            "product": (q.get("product", [""])[0] or ""),
            "note":    (q.get("note", [""])   [0] or ""),
            "start":   int(feat.location.start) + 1,   # GenBank → 1-based
            "end":     int(feat.location.end),
            "strand":  feat.location.strand or 1,
        })
    return rows


def match_panel(features: list[dict]) -> dict[str, dict]:
    """Match each gene in the panel to the best feature in the contig."""
    hits = {}
    for gene, pat, _ in GENE_PANEL:
        regex = re.compile(pat, re.IGNORECASE)
        for feat in features:
            blob = " ".join([feat["gene"], feat["product"], feat["note"]])
            if regex.search(blob):
                # Prefer CDS over gene-only feature
                prev = hits.get(gene)
                if prev is None or (feat["type"] == "CDS"
                                    and prev["type"] != "CDS"):
                    hits[gene] = feat
    return hits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cache-dir", default=None, type=Path,
                    help="Where to store .gb files "
                         "(default: assets/plasmid_refs/genbank_records/)")
    ap.add_argument("--fetch-genbank", action="store_true",
                    help="Force re-download even if cached.")
    ap.add_argument("--merge-blast", default=None, type=Path,
                    help="BLAST-derived coords TSV to compare against.")
    ap.add_argument("--combined-out", default=None, type=Path,
                    help="Where to write the combined GenBank+BLAST TSV.")
    args = ap.parse_args()

    cache_dir = args.cache_dir or (args.reference.parent / "plasmid_refs"
                                   / "genbank_records")

    contigs = contigs_from_fai(args.reference)
    # Drop the main KpSC chromosome — it's huge and we don't expect plasmid
    # AMR genes on it (chromosomal blaSHV is handled elsewhere).
    contigs = [c for c in contigs if c != "NC_016845.1"]

    print(f"Parsing GenBank annotations for {len(contigs)} contigs …")
    per_contig_hits: dict[str, dict[str, dict]] = {}
    for acc in contigs:
        try:
            gb = fetch_genbank(acc, cache_dir)
        except Exception as exc:                          # noqa: BLE001
            print(f"  {acc}: fetch failed — {exc}", file=sys.stderr)
            continue
        features = parse_features(gb)
        hits = match_panel(features)
        if hits:
            print(f"  {acc}: {', '.join(hits.keys())}")
            per_contig_hits[acc] = hits

    # ── Write the GenBank-derived TSV ─────────────────────────────────────────
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("gene\tcontig\tstart\tend\tabsent_threshold\tsource\n")
        for gene, _pat, abs_thr in GENE_PANEL:
            best = None
            for acc, hits in per_contig_hits.items():
                if gene in hits:
                    feat = hits[gene]
                    # Prefer CDS over other types
                    if best is None or (feat["type"] == "CDS"
                                        and best[1]["type"] != "CDS"):
                        best = (acc, feat)
            if best is None:
                fh.write(f"{gene}\t-\t-\t-\t{abs_thr}\tNO_GENBANK_HIT\n")
                continue
            acc, feat = best
            fh.write(f"{gene}\t{acc}\t{feat['start']}\t{feat['end']}\t"
                     f"{abs_thr}\tgenbank:{feat['type']}\n")
    print(f"\nWrote {args.out}")

    # ── Optional merge with BLAST ─────────────────────────────────────────────
    if args.merge_blast:
        if not args.combined_out:
            sys.exit("--combined-out required with --merge-blast")
        merge_with_blast(args.out, args.merge_blast, args.combined_out)


def merge_with_blast(gb_path: Path, blast_path: Path, out: Path) -> None:
    import pandas as pd
    gb = pd.read_csv(gb_path, sep="\t")
    bl = pd.read_csv(blast_path, sep="\t").rename(
        columns={"contig": "contig_blast", "start": "start_blast",
                 "end": "end_blast"}
    )[["gene", "contig_blast", "start_blast", "end_blast"]]
    merged = gb.merge(bl, on="gene", how="outer")
    merged["agreement"] = "?"
    for i, row in merged.iterrows():
        if row["contig"] == "-" or pd.isna(row.get("contig_blast")):
            merged.at[i, "agreement"] = "missing_one"
        elif row["contig"] != row["contig_blast"]:
            merged.at[i, "agreement"] = "different_contig"
        else:
            overlap = max(0,
                          min(int(row["end"]), int(row["end_blast"])) -
                          max(int(row["start"]), int(row["start_blast"])))
            length = int(row["end"]) - int(row["start"])
            merged.at[i, "agreement"] = "agree" if overlap > 0.5 * length else "diff_pos"
    merged.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}\n")
    print("agreement summary:")
    print(merged["agreement"].value_counts().to_string())


if __name__ == "__main__":
    main()
