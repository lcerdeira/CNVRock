"""
Generate Kleborate-based ground truth for KpSC chromosomal gene status.

Why Kleborate instead of AMRFinder for ompK35/ompK36/ramR:
  AMRFinder 2023 marks these genes "absent" (copy count 0) when the assembly
  sequence is divergent from its reference — even when the gene is physically
  present.  This means AMRFinder GT incorrectly labels many functional-gene
  samples as "deleted", producing CRR ~0.85 for the "deleted" category and
  making the evaluation look like FNR=1.0 (all deletions missed) when the model
  is actually correct.

  Kleborate uses curated KpSC-specific gene sequences and explicitly reports
  truncations vs full-length genes, which is the clinically meaningful
  distinction our model tries to detect.

Kleborate interpretation (ompK35, ompK36):
  - Full-length (e.g. "OmpK35", "OmpK36")     → cn = 1  (functional)
  - Truncated / premature stop / absent ("-")  → cn = 0  (loss of function)
  - Unknown / not reported                     → cn = -1 (uncallable)

ramR:
  Kleborate standard output does not report ramR directly.  This script uses
  blastn (BLAST+) against the assemblies with ramR extracted from HS11286.
  Requires:  blastn in PATH  AND  --assemblies-dir supplied.
  If BLAST is unavailable, ramR column is filled with -1 (uncallable).

Output format (compatible with kpsc_amrfinder_ground_truth.tsv):
    sample_id  ompK35  ompK36  ramR
    (integer: 1=functional, 0=absent/truncated, -1=uncallable)

Usage:
  # Step 1 — run Kleborate (v2 or v3) on assemblies:
  #   kleborate -a assemblies/*.fasta --resistance -o assets/kleborate_out.tsv  # v2
  #   kleborate -a assemblies/*.fasta --preset kpsc  -o assets/kleborate_out.tsv  # v3
  #
  # Step 2 — parse and optionally merge with AMRFinder GT:
  python data/setup/get_kleborate_gt.py \\
      --kleborate-tsv  assets/kleborate_out.tsv \\
      --amrfinder-gt   assets/kpsc_amrfinder_ground_truth.tsv \\
      --reference      assets/HS11286.fasta \\
      --assemblies-dir data/assemblies/ \\
      --out            assets/kpsc_kleborate_ground_truth.tsv

  # Without BLAST (ramR will be -1):
  python data/setup/get_kleborate_gt.py \\
      --kleborate-tsv  assets/kleborate_out.tsv \\
      --amrfinder-gt   assets/kpsc_amrfinder_ground_truth.tsv \\
      --out            assets/kpsc_kleborate_ground_truth.tsv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Gene coordinates in HS11286 (NC_016845.1, RefSeq annotation)
# ---------------------------------------------------------------------------

# ramR: KPHS_06060, 465 bp, plus strand
RAMR_CONTIG = "NC_016845.1"
RAMR_START  = 648627   # 1-based, inclusive
RAMR_END    = 649091
RAMR_BLAST_MIN_IDENTITY  = 80.0   # %
RAMR_BLAST_MIN_COVERAGE  = 60.0   # % of query length covered by HSP

# blaSHV-11: KPHS_25220, 861 bp, minus strand
BLASHV_CONTIG = "NC_016845.1"
BLASHV_START  = 2549403  # 1-based, inclusive
BLASHV_END    = 2550263
BLASHV_BLAST_MIN_IDENTITY = 80.0   # % — covers most blaSHV alleles (high family conservation)
BLASHV_BLAST_MIN_COVERAGE = 80.0   # % of query — require near-full-length hit per copy


# ---------------------------------------------------------------------------
# Kleborate parsing
# ---------------------------------------------------------------------------

def _parse_kleborate_v2(rows: list[dict]) -> dict[str, dict]:
    """Parse Kleborate ≤2.x output.  Strain column = 'strain'."""
    result = {}
    for row in rows:
        sid = row.get("strain", "").strip()
        if not sid:
            continue
        result[sid] = {
            "ompK35": _v2_porin_call(row.get("OmpK35", ""), "OmpK35"),
            "ompK36": _v2_porin_call(row.get("OmpK36", ""), "OmpK36"),
        }
    return result


def _v2_porin_call(value: str, gene_prefix: str) -> int:
    """Interpret a Kleborate v2 OmpK35/OmpK36 column value.

    Returns: 1 = functional, 0 = absent or truncated, -1 = not reported.
    """
    v = value.strip()
    if v == "-" or v == "":
        return 0   # absent
    if v == gene_prefix:
        return 1   # full-length functional
    # Truncations look like "OmpK35-6aa", "OmpK35*", "OmpK36-trunc" etc.
    if v.startswith(gene_prefix):
        return 0   # present but truncated → loss of function
    return -1      # unrecognised format


def _parse_kleborate_v3(rows: list[dict]) -> dict[str, dict]:
    """Parse Kleborate ≥3.x output.

    In v3 the KPSC preset combines both porins in a single column:
        klebsiella_pneumo_complex__amr__Omp_mutations

    Example values:
        '-'                                         → both functional
        'OmpK35:p.Glu42fs'                          → OmpK35 disrupted
        'OmpK36:p.Ala183fs'                         → OmpK36 disrupted
        'OmpK36:p.134_135insGlyAsp;OmpK35:p.Glu42fs' → both disrupted

    Any entry with 'OmpK35:' → ompK35=0 (disrupted).
    Any entry with 'OmpK36:' → ompK36=0 (disrupted).
    '-' or gene not mentioned → gene functional (cn=1).
    """
    # v3 KPSC column name for Omp mutations
    OMP_COL = "klebsiella_pneumo_complex__amr__Omp_mutations"

    result = {}
    for row in rows:
        sid = row.get("strain", "").strip()
        if not sid:
            continue

        omp_val = row.get(OMP_COL, "").strip()
        ompk35, ompk36 = _v3_omp_mutations_call(omp_val)
        result[sid] = {"ompK35": ompk35, "ompK36": ompk36}
    return result


def _v3_omp_mutations_call(omp_value: str) -> tuple[int, int]:
    """Parse Kleborate v3 Omp_mutations field → (ompK35_cn, ompK36_cn).

    Returns 1=functional, 0=disrupted for each gene.
    Logic: if the gene name appears in the field, it has a mutation → disrupted.
    If the field is '-' or the gene is not mentioned, it is functional.
    """
    v = omp_value.strip()
    if v == "-" or v == "":
        return 1, 1   # no mutations → both functional
    ompk35 = 0 if "OmpK35:" in v else 1
    ompk36 = 0 if "OmpK36:" in v else 1
    return ompk35, ompk36


def parse_kleborate(tsv_path: str) -> dict[str, dict]:
    """Auto-detect Kleborate version from header and parse accordingly."""
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise RuntimeError(f"Kleborate output is empty: {tsv_path}")

    header = set(rows[0].keys())
    # v2: separate 'OmpK35' / 'OmpK36' columns (capital K)
    # v3: combined 'klebsiella_pneumo_complex__amr__Omp_mutations' column
    if "OmpK35" in header or "OmpK36" in header:
        print(f"Detected Kleborate v2 format ({tsv_path})", flush=True)
        return _parse_kleborate_v2(rows)
    elif any("Omp_mutations" in h for h in header):
        print(f"Detected Kleborate v3 KPSC preset format ({tsv_path})", flush=True)
        return _parse_kleborate_v3(rows)
    else:
        print(f"Detected Kleborate v3 format ({tsv_path}) — attempting v3 parse", flush=True)
        return _parse_kleborate_v3(rows)


# ---------------------------------------------------------------------------
# BLAST-based ramR detection
# ---------------------------------------------------------------------------

def _extract_ramr_sequence(reference_fasta: str) -> str:
    """Extract ramR nucleotide sequence from HS11286 reference using samtools faidx."""
    region = f"{RAMR_CONTIG}:{RAMR_START}-{RAMR_END}"
    try:
        result = subprocess.run(
            ["samtools", "faidx", reference_fasta, region],
            capture_output=True, text=True, check=True,
        )
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: samtools faidx failed ({e}). Will not run BLAST for ramR.", file=sys.stderr)
        return ""


def _blast_ramr(query_fasta: str, assemblies_dir: str) -> dict[str, int]:
    """Run blastn of ramR query against all assemblies.  Returns {sample_id: cn}."""
    query_len = RAMR_END - RAMR_START + 1
    ramr_calls = {}

    fasta_files = [
        f for f in os.listdir(assemblies_dir)
        if f.endswith((".fasta", ".fa", ".fna"))
    ]
    if not fasta_files:
        print(f"WARNING: no FASTA files found in {assemblies_dir}", file=sys.stderr)
        return ramr_calls

    print(f"BLASTing ramR against {len(fasta_files)} assemblies in {assemblies_dir} ...", flush=True)

    for fname in fasta_files:
        sample_id = fname.rsplit(".", 1)[0]
        assembly  = os.path.join(assemblies_dir, fname)
        try:
            result = subprocess.run(
                [
                    "blastn",
                    "-query", query_fasta,
                    "-subject", assembly,
                    "-outfmt", "6 pident length qstart qend",
                    "-perc_identity", str(RAMR_BLAST_MIN_IDENTITY),
                    "-max_hsps", "1",
                ],
                capture_output=True, text=True, check=True,
            )
            best = None
            for line in result.stdout.strip().splitlines():
                pct_id, hsp_len, q_start, q_end = map(float, line.split("\t"))
                coverage = (q_end - q_start + 1) / query_len * 100.0
                if pct_id >= RAMR_BLAST_MIN_IDENTITY and coverage >= RAMR_BLAST_MIN_COVERAGE:
                    best = 1
                    break
            ramr_calls[sample_id] = best if best is not None else 0
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  WARNING: BLAST failed for {sample_id}: {e}", file=sys.stderr)
            ramr_calls[sample_id] = -1

    n_present = sum(1 for v in ramr_calls.values() if v == 1)
    n_absent  = sum(1 for v in ramr_calls.values() if v == 0)
    print(f"ramR BLAST: present={n_present}  absent={n_absent}  failed={len(ramr_calls)-n_present-n_absent}", flush=True)
    return ramr_calls


# ---------------------------------------------------------------------------
# BLAST-based blaSHV copy counting
# ---------------------------------------------------------------------------

def _extract_blashv_sequence(reference_fasta: str) -> str:
    """Extract blaSHV-11 nucleotide sequence from HS11286 reference using samtools faidx."""
    region = f"{BLASHV_CONTIG}:{BLASHV_START}-{BLASHV_END}"
    try:
        result = subprocess.run(
            ["samtools", "faidx", reference_fasta, region],
            capture_output=True, text=True, check=True,
        )
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: samtools faidx failed for blaSHV ({e}).", file=sys.stderr)
        return ""


def _blast_blashv(query_fasta: str, assemblies_dir: str) -> dict[str, int]:
    """Run blastn of blaSHV-11 query against all assemblies.

    Returns {sample_id: copy_count} where copy_count is the number of
    non-overlapping HSPs meeting identity and coverage thresholds.
    cn=1 → single chromosomal copy (normal).
    cn≥2 → amplification.
    cn=0 → blaSHV absent (rare in K. pneumoniae).
    """
    query_len = BLASHV_END - BLASHV_START + 1
    blashv_calls: dict[str, int] = {}

    fasta_files = [
        f for f in os.listdir(assemblies_dir)
        if f.endswith((".fasta", ".fa", ".fna"))
    ]
    if not fasta_files:
        print(f"WARNING: no FASTA files found in {assemblies_dir}", file=sys.stderr)
        return blashv_calls

    print(f"BLASTing blaSHV against {len(fasta_files)} assemblies in {assemblies_dir} ...", flush=True)

    for fname in fasta_files:
        sample_id = fname.rsplit(".", 1)[0]
        assembly  = os.path.join(assemblies_dir, fname)
        try:
            result = subprocess.run(
                [
                    "blastn",
                    "-query",  query_fasta,
                    "-subject", assembly,
                    "-outfmt", "6 pident length qstart qend sstart send",
                    "-perc_identity", str(BLASHV_BLAST_MIN_IDENTITY),
                    # Allow up to 20 HSPs — amplified strains may have many copies
                    "-max_hsps", "20",
                ],
                capture_output=True, text=True, check=True,
            )
            # Count non-overlapping hits on the subject (assembly) side
            # Each copy of blaSHV occupies a distinct genomic region → distinct sstart/send
            hits: list[tuple[int, int]] = []  # (sstart, send) of accepted HSPs
            for line in result.stdout.strip().splitlines():
                pct_id, hsp_len, q_start, q_end, s_start, s_end = map(float, line.split("\t"))
                coverage = (q_end - q_start + 1) / query_len * 100.0
                if pct_id < BLASHV_BLAST_MIN_IDENTITY or coverage < BLASHV_BLAST_MIN_COVERAGE:
                    continue
                s_lo, s_hi = int(min(s_start, s_end)), int(max(s_start, s_end))
                # Suppress if >50% overlaps an already-counted hit
                overlaps = False
                for (h_lo, h_hi) in hits:
                    overlap = max(0, min(s_hi, h_hi) - max(s_lo, h_lo))
                    if overlap / (s_hi - s_lo + 1) > 0.5:
                        overlaps = True
                        break
                if not overlaps:
                    hits.append((s_lo, s_hi))
            blashv_calls[sample_id] = len(hits)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  WARNING: BLAST failed for {sample_id}: {e}", file=sys.stderr)
            blashv_calls[sample_id] = -1

    n1 = sum(1 for v in blashv_calls.values() if v == 1)
    n2 = sum(1 for v in blashv_calls.values() if v >= 2)
    n0 = sum(1 for v in blashv_calls.values() if v == 0)
    nf = sum(1 for v in blashv_calls.values() if v == -1)
    print(f"blaSHV BLAST: single={n1}  amplified(≥2)={n2}  absent={n0}  failed={nf}", flush=True)
    return blashv_calls


# ---------------------------------------------------------------------------
# Merge with AMRFinder GT
# ---------------------------------------------------------------------------

def load_amrfinder_gt(tsv_path: str) -> dict[str, dict]:
    """Load AMRFinder GT as {sample_id: {gene: cn, ...}}."""
    gt = {}
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sid = row.pop("sample_id")
            gt[sid] = {k: (int(v) if v not in ("", "NA", "nan") else -1)
                       for k, v in row.items()}
    return gt


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kleborate-tsv",  required=True,
                        help="Kleborate output TSV (v2 or v3).")
    parser.add_argument("--amrfinder-gt",   default=None,
                        help="AMRFinder GT TSV to merge blaSHV/blaKPC/blaCTX-M/blaNDM "
                             "columns into the output.  Optional.")
    parser.add_argument("--reference",      default=None,
                        help="HS11286.fasta path — needed for samtools faidx to extract "
                             "ramR and blaSHV sequences.  Skip if --assemblies-dir not supplied.")
    parser.add_argument("--assemblies-dir", default=None,
                        help="Directory of per-sample assembly FASTA files for BLAST-based "
                             "ramR and blaSHV copy-number calling.  "
                             "Files must be named <sample_id>.fasta/.fa/.fna.")
    parser.add_argument("--out",            required=True,
                        help="Output TSV path.")
    args = parser.parse_args()

    # ── Parse Kleborate ───────────────────────────────────────────────────────
    kleborate = parse_kleborate(args.kleborate_tsv)
    print(f"Parsed {len(kleborate):,} samples from Kleborate TSV.", flush=True)

    # ── ramR + blaSHV via BLAST (optional) ───────────────────────────────────
    ramr_calls:   dict[str, int] = {}
    blashv_calls: dict[str, int] = {}
    if args.assemblies_dir and args.reference:
        # ramR — presence/absence
        ramr_seq = _extract_ramr_sequence(args.reference)
        if ramr_seq:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tf:
                tf.write(ramr_seq)
                query_path = tf.name
            try:
                ramr_calls = _blast_ramr(query_path, args.assemblies_dir)
            finally:
                os.unlink(query_path)
        # blaSHV — copy number
        blashv_seq = _extract_blashv_sequence(args.reference)
        if blashv_seq:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tf:
                tf.write(blashv_seq)
                query_path = tf.name
            try:
                blashv_calls = _blast_blashv(query_path, args.assemblies_dir)
            finally:
                os.unlink(query_path)
    elif args.assemblies_dir:
        print("WARNING: --reference not supplied; skipping BLAST for ramR and blaSHV.", file=sys.stderr)
    else:
        print("--assemblies-dir not supplied; ramR and blaSHV will use AMRFinder GT.", flush=True)

    # ── Load AMRFinder GT for non-overridden genes ────────────────────────────
    # Columns we override with BLAST or Kleborate:
    #   ompK35, ompK36 — Kleborate (sequence-level)
    #   ramR           — BLAST (if assemblies-dir supplied)
    #   blaSHV         — BLAST copy count (if assemblies-dir supplied), else AMRFinder
    BLAST_OVERRIDE_COLS = {"ompK35", "ompK36", "ramR"}
    if blashv_calls:
        BLAST_OVERRIDE_COLS.add("blaSHV")

    amr_gt: dict[str, dict] = {}
    amr_extra_cols: list[str] = []
    if args.amrfinder_gt:
        amr_gt = load_amrfinder_gt(args.amrfinder_gt)
        sample_amr = next(iter(amr_gt.values()), {})
        amr_extra_cols = [c for c in sample_amr if c not in BLAST_OVERRIDE_COLS]
        # If blaSHV not in BLAST override, carry it from AMRFinder
        if "blaSHV" not in BLAST_OVERRIDE_COLS and "blaSHV" in sample_amr:
            amr_extra_cols = ["blaSHV"] + [c for c in amr_extra_cols if c != "blaSHV"]
        print(f"Loaded AMRFinder GT ({len(amr_gt):,} samples); "
              f"carrying forward: {amr_extra_cols}", flush=True)

    # ── Build merged output ───────────────────────────────────────────────────
    all_sample_ids = sorted(set(kleborate) | set(amr_gt) | set(blashv_calls))
    # blaSHV goes first (primary chromosomal target), then porins, ramR, then rest
    blast_cols = ["blaSHV"] if blashv_calls else []
    out_cols = ["sample_id"] + blast_cols + ["ompK35", "ompK36", "ramR"] + \
               [c for c in amr_extra_cols if c not in blast_cols]

    n_written = 0
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for sid in all_sample_ids:
            kleb = kleborate.get(sid, {})
            amr  = amr_gt.get(sid, {})
            row: dict = {
                "sample_id": sid,
                "ompK35":    kleb.get("ompK35", -1),
                "ompK36":    kleb.get("ompK36", -1),
                "ramR":      ramr_calls.get(sid, -1),
            }
            if blashv_calls:
                row["blaSHV"] = blashv_calls.get(sid, -1)
            for col in amr_extra_cols:
                if col not in row:
                    row[col] = amr.get(col, -1)
            writer.writerow(row)
            n_written += 1

    print(f"\nSaved Kleborate GT ({n_written:,} samples) → {args.out}", flush=True)

    # ── Diagnostics ───────────────────────────────────────────────────────────
    ompk35_counts = {"functional": 0, "truncated": 0, "uncallable": 0}
    ompk36_counts = {"functional": 0, "truncated": 0, "uncallable": 0}
    for d in kleborate.values():
        for gene, counts_dict in [("ompK35", ompk35_counts), ("ompK36", ompk36_counts)]:
            v = d.get(gene, -1)
            if v == 1:
                counts_dict["functional"] += 1
            elif v == 0:
                counts_dict["truncated"] += 1
            else:
                counts_dict["uncallable"] += 1
    print(f"  ompK35: {ompk35_counts}")
    print(f"  ompK36: {ompk36_counts}")
    if ramr_calls:
        n_p = sum(1 for v in ramr_calls.values() if v == 1)
        n_a = sum(1 for v in ramr_calls.values() if v == 0)
        n_u = sum(1 for v in ramr_calls.values() if v == -1)
        print(f"  ramR: functional={n_p}  absent={n_a}  uncallable={n_u}")
    if blashv_calls:
        n1 = sum(1 for v in blashv_calls.values() if v == 1)
        n2 = sum(1 for v in blashv_calls.values() if v >= 2)
        n0 = sum(1 for v in blashv_calls.values() if v == 0)
        print(f"  blaSHV copies: 1={n1}  ≥2(amp)={n2}  0(absent)={n0}")


if __name__ == "__main__":
    main()
