"""
Download CABBAGE phenotypic AMR data and join to ATB sample accessions.

CABBAGE (Comprehensive Assessment of Bacterial-Based AMR prediction from
Genotypes) is the largest public AST dataset, hosted at EBI AMR portal:
  https://www.ebi.ac.uk/amr/

Release files: https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-12/

Data schema (phenotype.parquet):
  BioSample_ID, SRA_accession, assembly_ID
  organism, genus, species
  antibiotic_name, ast_standard (CLSI/EUCAST)
  measurement, measurement_sign, measurement_units  — MIC value
  resistance_phenotype    — susceptible / intermediate / resistant
  Updated_phenotype_CLSI, Updated_phenotype_EUCAST  — normalised calls
  collection_year, ISO_country_code, isolation_source, host, ...

Output (assets/cabbage_kpsc_phenotypes.tsv):
  sample_id (= SRA run accession), antibiotic_name, ast_standard,
  measurement, measurement_sign, measurement_units, resistance_phenotype,
  updated_phenotype, collection_year, country, isolation_source

Usage:
  python data/setup/get_cabbage_phenotypes.py [--release 2025-12]

The output joins to CNVRock results via the sample_id column.
To find overlap with current cohort:
  python data/setup/get_cabbage_phenotypes.py --report-overlap

NOTE: 22,774 of ~26,265 CABBAGE KpSC samples appear in ATB (AMRFinder+
available). The remaining ~3,491 are not in ATB — their assemblies must be
downloaded separately if needed (see --report-overlap for the list).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd

REPO_DIR  = Path(__file__).resolve().parents[2]
ASSETS    = REPO_DIR / "assets"
CACHE_DIR = ASSETS / "cabbage_cache"

EBI_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases"

# KpSC organism names as they appear in CABBAGE
KPSC_ORGANISMS = {
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae",
    "Klebsiella variicola",
    "Klebsiella africana",
    "Klebsiella pneumoniae subsp. pneumoniae",
    "Klebsiella pneumoniae subsp. ozaenae",
    "Klebsiella pneumoniae subsp. rhinoscleromatis",
}

# Antibiotics of clinical relevance for KpSC AMR — used for focused output
# (all antibiotics are kept in the TSV; this just adds a priority flag)
PRIORITY_ANTIBIOTICS = {
    # Carbapenems
    "imipenem", "meropenem", "ertapenem", "doripenem",
    # 3rd/4th-gen cephalosporins
    "ceftriaxone", "cefotaxime", "ceftazidime", "cefepime",
    # Quinolones
    "ciprofloxacin", "levofloxacin",
    # Aminoglycosides
    "amikacin", "gentamicin", "tobramycin",
    # Tigecycline, colistin
    "tigecycline", "colistin", "polymyxin b",
    # Beta-lactam/beta-lactamase inhibitors
    "piperacillin-tazobactam", "amoxicillin-clavulanic acid", "ceftazidime-avibactam",
    "meropenem-vaborbactam",
    # Aztreonam
    "aztreonam",
    # Trimethoprim
    "trimethoprim", "trimethoprim-sulfamethoxazole",
}


def _latest_release(ftp_base: str = EBI_FTP_BASE) -> str:
    """Return the most recent release directory name."""
    import urllib.request
    import re
    with urllib.request.urlopen(f"{ftp_base}/") as r:
        html = r.read().decode()
    releases = sorted(re.findall(r'href="(\d{4}-\d{2})/"', html))
    if not releases:
        raise RuntimeError("Could not find any CABBAGE releases")
    return releases[-1]


def _download_cached(url: str, cache_path: Path, label: str = "") -> Path:
    """Download url to cache_path if not already present."""
    if cache_path.exists():
        print(f"  Using cached {cache_path.name}", flush=True)
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Downloading {label or url} …", flush=True)
    urlretrieve(url, cache_path)
    size_mb = cache_path.stat().st_size / 1024**2
    print(f"  → {cache_path}  ({size_mb:.1f} MB)", flush=True)
    return cache_path


def load_cabbage_phenotypes(release: str) -> pd.DataFrame:
    """Download (once) and return the full CABBAGE phenotype table."""
    url = f"{EBI_FTP_BASE}/{release}/phenotype.parquet"
    cache = CACHE_DIR / f"phenotype_{release}.parquet"
    path = _download_cached(url, cache, f"CABBAGE phenotype {release}")
    print(f"  Loading phenotype table …", flush=True)
    df = pd.read_parquet(path)
    print(f"  {len(df):,} rows × {len(df.columns)} columns", flush=True)
    return df


def filter_kpsc(df: pd.DataFrame) -> pd.DataFrame:
    """Retain only KpSC rows."""
    mask = df["organism"].isin(KPSC_ORGANISMS) | df["species"].isin(KPSC_ORGANISMS)
    out = df[mask].copy()
    print(f"  KpSC rows: {len(out):,} (from {len(df):,} total)", flush=True)
    return out


def normalise(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select and normalise output columns.

    Join key for CNVRock: sample_id = SRA_accession.
    Multiple rows per sample (one per antibiotic test).
    """
    # Pick the best available breakpoint interpretation
    df["updated_phenotype"] = df["Updated_phenotype_EUCAST"].where(
        df["Updated_phenotype_EUCAST"].notna(),
        df["Updated_phenotype_CLSI"],
    )

    # Flag priority antibiotics
    df["priority_antibiotic"] = (
        df["antibiotic_name"].str.lower().isin(PRIORITY_ANTIBIOTICS)
    )

    out = df[[
        "SRA_accession",
        "BioSample_ID",
        "assembly_ID",
        "organism",
        "antibiotic_name",
        "priority_antibiotic",
        "ast_standard",
        "laboratory_typing_method",
        "measurement",
        "measurement_sign",
        "measurement_units",
        "resistance_phenotype",
        "updated_phenotype",
        "collection_year",
        "ISO_country_code",
        "country",
        "isolation_source",
        "isolation_source_category",
        "host",
        "geographical_region",
    ]].rename(columns={"SRA_accession": "sample_id"})

    out["sample_id"] = out["sample_id"].astype(str).str.strip()
    return out


def report_overlap(kpsc: pd.DataFrame, atb_sra: pd.Series) -> None:
    """Print overlap statistics between CABBAGE and ATB."""
    cab_ids = set(kpsc["sample_id"].dropna())
    atb_ids = set(atb_sra.dropna())

    in_both   = cab_ids & atb_ids
    cab_only  = cab_ids - atb_ids
    atb_only  = atb_ids - cab_ids

    print(f"\nCABBAGE KpSC unique samples   : {len(cab_ids):>7,}")
    print(f"ATB KpSC samples              : {len(atb_ids):>7,}")
    print(f"In both (genotype+phenotype)  : {len(in_both):>7,}")
    print(f"CABBAGE only (need assembly)  : {len(cab_only):>7,}")
    print(f"ATB only (no phenotype data)  : {len(atb_only):>7,}")

    if cab_only:
        cab_only_path = ASSETS / "cabbage_only_accessions.txt"
        pd.Series(sorted(cab_only)).to_csv(cab_only_path, index=False, header=False)
        print(f"\nCABBAGE-only accessions → {cab_only_path}")
        print("  These need assemblies downloaded separately for Kleborate GT.")

    # Antibiotic coverage in the overlap
    overlap_df = kpsc[kpsc["sample_id"].isin(in_both)]
    print("\nTop antibiotics in overlap (unique sample counts):")
    ab_counts = (overlap_df.groupby("antibiotic_name")["sample_id"]
                 .nunique().sort_values(ascending=False).head(15))
    for ab, n in ab_counts.items():
        flag = " ***" if ab.lower() in PRIORITY_ANTIBIOTICS else ""
        print(f"  {ab:<45} {n:>6,}{flag}")


def main():
    parser = argparse.ArgumentParser(
        description="Download CABBAGE KpSC phenotypes and join to ATB accessions"
    )
    parser.add_argument("--release", default=None,
                        help="CABBAGE release (e.g. 2025-12). Default: latest.")
    parser.add_argument("--report-overlap", action="store_true",
                        help="Print overlap stats with current ATB cohort and exit")
    parser.add_argument("--out",
                        default=str(ASSETS / "cabbage_kpsc_phenotypes.tsv"),
                        help="Output TSV path")
    parser.add_argument("--priority-only", action="store_true",
                        help="Output only rows for priority antibiotics")
    args = parser.parse_args()

    # ── Resolve release ───────────────────────────────────────────────────────
    if args.release is None:
        print("Resolving latest CABBAGE release …", flush=True)
        try:
            args.release = _latest_release()
        except Exception as e:
            args.release = "2025-12"
            print(f"  WARNING: could not auto-detect release ({e}); using {args.release}",
                  file=sys.stderr)
    print(f"Using CABBAGE release: {args.release}", flush=True)

    # ── Load and filter ───────────────────────────────────────────────────────
    raw   = load_cabbage_phenotypes(args.release)
    kpsc  = filter_kpsc(raw)
    kpsc  = normalise(kpsc)

    # ── Overlap report ────────────────────────────────────────────────────────
    if args.report_overlap:
        # Load ATB run accessions from expansion list or existing cohort
        atb_path = ASSETS / "kpsc_expansion_sra_accessions.txt"
        if not atb_path.exists():
            atb_path = ASSETS / "kpsc_bam_accessions.txt"
        if atb_path.exists():
            atb_sra = pd.read_csv(atb_path, header=None, names=["run_accession"])["run_accession"]
        else:
            atb_sra = pd.Series(dtype=str)
            print("WARNING: no ATB accession file found; overlap computed against empty set",
                  file=sys.stderr)
        report_overlap(kpsc, atb_sra)
        # Fall through to also write the output file

    # ── Optional filter to priority antibiotics ───────────────────────────────
    if args.priority_only:
        kpsc = kpsc[kpsc["priority_antibiotic"]]
        print(f"After priority antibiotic filter: {len(kpsc):,} rows", flush=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    unique_samples = kpsc["sample_id"].nunique()
    unique_ab      = kpsc["antibiotic_name"].nunique()
    print(f"\nKpSC phenotype records : {len(kpsc):,}")
    print(f"Unique samples         : {unique_samples:,}")
    print(f"Unique antibiotics     : {unique_ab:,}")

    if "measurement" in kpsc.columns:
        has_mic = kpsc["measurement"].notna().sum()
        print(f"Records with MIC value : {has_mic:,} ({has_mic/len(kpsc)*100:.0f}%)")

    print("\nResistance phenotype breakdown:")
    print(kpsc["resistance_phenotype"].value_counts().to_string())

    print("\nOrganism breakdown:")
    print(kpsc["organism"].value_counts().to_string())

    # ── Write output ──────────────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kpsc.to_csv(out_path, sep="\t", index=False)
    print(f"\nCABBAGE KpSC phenotypes → {out_path}", flush=True)
    print(f"Join to CNVRock results on: sample_id ↔ sample_id", flush=True)
    print(f"\nNext: to link PCN to MIC, run:")
    print(f"  python data/setup/join_pcn_phenotypes.py \\")
    print(f"      --pcn  data/results/<exp>/plasmid_gene_calls.tsv \\")
    print(f"      --pheno {out_path}")


if __name__ == "__main__":
    main()
