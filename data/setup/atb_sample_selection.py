"""
AllTheBacteria KpSC sample selection utility.

Downloads AllTheBacteria metadata, filters for KpSC species, applies quality
thresholds, and outputs:
  - assets/kpsc_sra_accessions.txt   — one SRA run accession per line
  - assets/kpsc_sample_metadata.tsv  — sample_id, Species, ST, N50, coverage ...
  - assets/kpsc_amrfinder_ground_truth.tsv  (--summarise-amrfinder mode only)

Usage
-----
# Step 1: Download metadata and write sample list + metadata TSV
python data/setup/atb_sample_selection.py \
    --atb-metadata atb_metadata.tsv \
    --out-accessions assets/kpsc_sra_accessions.txt \
    --out-metadata   assets/kpsc_sample_metadata.tsv \
    [--min-coverage 30] [--max-samples 10000]

# Step 2: After running AMRFinder+ on each assembly, aggregate results
python data/setup/atb_sample_selection.py \
    --summarise-amrfinder gt/ \
    --out-gt assets/kpsc_amrfinder_ground_truth.tsv

AllTheBacteria metadata
-----------------------
The metadata TSV is distributed with AllTheBacteria releases:
  https://github.com/AllTheBacteria/AllTheBacteria

Download a pre-built metadata file or pull from Zenodo (see the ATB GitHub).
Expected columns (column names may vary by ATB release):
  sample_id / run_accession  — SRA run accession (ERR.../SRR...)
  species                    — e.g. "Klebsiella pneumoniae"
  genus                      — "Klebsiella"
  n50                        — assembly N50 in bp
  completeness               — CheckM2 completeness (0–100)
  contamination              — CheckM2 contamination (0–100)
  coverage                   — estimated read coverage depth
  ST                         — sequence type from Kleborate (optional)

KpSC species filter (based on Kleborate species calls):
  Klebsiella pneumoniae
  Klebsiella quasipneumoniae subsp. quasipneumoniae
  Klebsiella quasipneumoniae subsp. similipneumoniae
  Klebsiella variicola subsp. variicola
  Klebsiella variicola subsp. tropica
  Klebsiella africana
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


KPSC_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. similipneumoniae",
    "Klebsiella variicola",
    "Klebsiella variicola subsp. variicola",
    "Klebsiella variicola subsp. tropica",
    "Klebsiella africana",
}

# AMRFinder genes of interest — must match GENES_OF_INTEREST in 06_gene_cnv_caller.py
AMRFINDER_GENES = ["blaSHV", "ompK35", "ompK36", "ramR"]


# ---------------------------------------------------------------------------
# Step 1: Filter ATB metadata → sample list + metadata TSV
# ---------------------------------------------------------------------------

def select_samples(
    atb_metadata_path: str,
    out_accessions: str,
    out_metadata: str,
    min_coverage: float = 30.0,
    min_completeness: float = 90.0,
    max_contamination: float = 5.0,
    min_n50: int = 50_000,
    max_samples: int | None = None,
) -> pd.DataFrame:
    print(f"Loading ATB metadata: {atb_metadata_path}", flush=True)
    df = pd.read_csv(atb_metadata_path, sep="\t", low_memory=False)
    print(f"  Total rows: {len(df):,}", flush=True)

    # Normalise column names — ATB releases use different names across versions
    col_map = {
        "run_accession": "sample_id",
        "accession":     "sample_id",
        "species_name":  "Species",
        "species":       "Species",
        "genome_n50":    "N50",
        "n50":           "N50",
        "checkm2_completeness": "completeness",
        "checkm2_contamination": "contamination",
        "depth":         "coverage",
        "mean_coverage": "coverage",
        "sequence_type": "ST",
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    if "sample_id" not in df.columns:
        raise ValueError(
            "Could not find sample/run accession column in metadata. "
            f"Columns present: {list(df.columns)}"
        )
    if "Species" not in df.columns:
        raise ValueError(
            "Could not find species column in metadata. "
            f"Columns present: {list(df.columns)}"
        )

    # Species filter
    kpsc_mask = df["Species"].isin(KPSC_SPECIES)
    df = df[kpsc_mask].copy()
    print(f"  After KpSC species filter: {len(df):,}", flush=True)

    # Quality filters
    if "completeness" in df.columns:
        df = df[df["completeness"] >= min_completeness]
        print(f"  After completeness ≥ {min_completeness}: {len(df):,}", flush=True)
    if "contamination" in df.columns:
        df = df[df["contamination"] <= max_contamination]
        print(f"  After contamination ≤ {max_contamination}: {len(df):,}", flush=True)
    if "N50" in df.columns:
        df = df[df["N50"] >= min_n50]
        print(f"  After N50 ≥ {min_n50:,}: {len(df):,}", flush=True)
    if "coverage" in df.columns:
        df = df[df["coverage"] >= min_coverage]
        print(f"  After coverage ≥ {min_coverage}×: {len(df):,}", flush=True)

    if max_samples is not None and len(df) > max_samples:
        df = df.sample(n=max_samples, random_state=42)
        print(f"  Subsampled to {max_samples:,} samples.", flush=True)

    # Write outputs
    Path(out_accessions).parent.mkdir(parents=True, exist_ok=True)
    df["sample_id"].to_csv(out_accessions, index=False, header=False)
    print(f"Wrote {len(df):,} SRA accessions → {out_accessions}", flush=True)

    meta_cols = ["sample_id", "Species"]
    for col in ["ST", "N50", "coverage", "completeness", "contamination"]:
        if col in df.columns:
            meta_cols.append(col)
    df[meta_cols].to_csv(out_metadata, sep="\t", index=False)
    print(f"Wrote sample metadata ({len(df):,} rows) → {out_metadata}", flush=True)

    # Summary by species
    print("\nSpecies breakdown:")
    for sp, n in df["Species"].value_counts().items():
        print(f"  {sp:<50} {n:>6,}")

    return df


# ---------------------------------------------------------------------------
# Step 2: Aggregate AMRFinder+ results per sample → ground truth TSV
# ---------------------------------------------------------------------------

def summarise_amrfinder(amrfinder_dir: str, out_gt: str) -> None:
    """Read per-sample AMRFinder+ TSV files and aggregate gene copy counts.

    Each file is expected to be the direct output of:
        amrfinder --plus -O Klebsiella_pneumoniae -n assembly.fasta --output sample.tsv

    The function counts how many times each gene of interest appears per sample
    (copy count). It outputs a TSV with columns:
        sample_id, blaSHV, ompK35, ompK36, ramR

    For blaSHV: count ≥ 2 → amplification event.
    For ompK35/ompK36/ramR: count = 0 → deletion/disruption event.
    """
    amrfinder_dir = Path(amrfinder_dir)
    tsv_files = list(amrfinder_dir.glob("*.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No .tsv files found in {amrfinder_dir}")

    print(f"Processing {len(tsv_files):,} AMRFinder+ outputs…", flush=True)

    rows = []
    for tsv_path in tsv_files:
        sample_id = tsv_path.stem
        try:
            amr = pd.read_csv(tsv_path, sep="\t", low_memory=False)
        except Exception as exc:
            print(f"  WARN: skipping {tsv_path.name}: {exc}", flush=True)
            continue

        # AMRFinder+ gene name column — usually "Gene symbol" or "Name"
        gene_col = next(
            (c for c in amr.columns if "gene" in c.lower() or c == "Name"),
            None
        )
        if gene_col is None:
            print(f"  WARN: no gene name column in {tsv_path.name}, skipping.", flush=True)
            continue

        row = {"sample_id": sample_id}
        for gene in AMRFINDER_GENES:
            # Count rows where gene name starts with the target (handles blaSHV-11, blaSHV-12 etc.)
            count = int((amr[gene_col].str.startswith(gene, na=False)).sum())
            row[gene] = count
        rows.append(row)

    gt = pd.DataFrame(rows, columns=["sample_id"] + AMRFINDER_GENES)
    gt.to_csv(out_gt, sep="\t", index=False)
    print(f"Ground truth written ({len(gt):,} samples) → {out_gt}", flush=True)

    # Summary
    print("\nGround truth summary:")
    for gene in AMRFINDER_GENES:
        n_amp = (gt[gene] >= 2).sum()
        n_del = (gt[gene] == 0).sum()
        n_ok  = (gt[gene] == 1).sum()
        print(f"  {gene:<10}  n=1: {n_ok:>6,}  n=0 (del): {n_del:>6,}  n≥2 (amp): {n_amp:>6,}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="AllTheBacteria KpSC sample selection and ground truth aggregation"
    )
    subparsers = parser.add_subparsers(dest="command")

    sel = subparsers.add_parser("select", help="Filter ATB metadata and write sample list")
    sel.add_argument("--atb-metadata",    required=True, help="Path to AllTheBacteria metadata TSV")
    sel.add_argument("--out-accessions",  default="assets/kpsc_sra_accessions.txt")
    sel.add_argument("--out-metadata",    default="assets/kpsc_sample_metadata.tsv")
    sel.add_argument("--min-coverage",    type=float, default=30.0)
    sel.add_argument("--min-completeness", type=float, default=90.0)
    sel.add_argument("--max-contamination", type=float, default=5.0)
    sel.add_argument("--min-n50",         type=int,   default=50_000)
    sel.add_argument("--max-samples",     type=int,   default=None)

    agg = subparsers.add_parser("summarise-amrfinder",
                                 help="Aggregate AMRFinder+ TSVs into ground truth")
    agg.add_argument("amrfinder_dir", help="Directory containing per-sample AMRFinder+ TSVs")
    agg.add_argument("--out-gt", default="assets/kpsc_amrfinder_ground_truth.tsv")

    args = parser.parse_args()

    if args.command == "select":
        select_samples(
            atb_metadata_path = args.atb_metadata,
            out_accessions    = args.out_accessions,
            out_metadata      = args.out_metadata,
            min_coverage      = args.min_coverage,
            min_completeness  = args.min_completeness,
            max_contamination = args.max_contamination,
            min_n50           = args.min_n50,
            max_samples       = args.max_samples,
        )
    elif args.command == "summarise-amrfinder":
        summarise_amrfinder(
            amrfinder_dir = args.amrfinder_dir,
            out_gt        = args.out_gt,
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
