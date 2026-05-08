"""
Phase 2 cohort expansion — generate download lists for all HQ KpSC samples
in AllTheBacteria that are not yet in the current cohort.

Outputs (written to assets/):
  kpsc_expansion_metadata.tsv      — sample_id, species, ST, N50, assembly_url
  kpsc_expansion_sra_accessions.txt — SRA run accessions for BAM download (download_sra.sh)
  kpsc_expansion_run_ftp_map.tsv   — run_accession, ftp_r1, ftp_r2  (needed by download_sra.sh)
  kpsc_expansion_assembly_urls.tsv — sample_accession, run_accession, aws_url (for Kleborate)

Prerequisites (run locally, ATB CLI installed):
  atb fetch --tables assembly.parquet,mlst.parquet --yes

Then transfer outputs to HPC:
  rsync assets/kpsc_expansion_* pryor:CNVRock/assets/

Usage:
  python data/setup/expand_cohort.py [--dry-run] [--max-samples N]

HPC next steps (after transfer):
  # Append new accessions to existing list:
  cat assets/kpsc_bam_accessions.txt assets/kpsc_expansion_sra_accessions.txt \\
      | sort -u > assets/kpsc_all_sra_accessions.txt

  # Download FASTQs + align:
  N=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)
  sbatch --array=1-${N}%50 hpc/download_sra.sh  # uses kpsc_expansion_sra_accessions.txt

  # Download assemblies for Kleborate:
  sbatch --array=1-${N}%50 hpc/download_assemblies_s3.sh

  # Extract read counts (chromosomal):
  sbatch --array=1-${N}%100 hpc/collect_readcounts.sh
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import requests

# ── KpSC species included in the expansion ───────────────────────────────────
KPSC_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae",
    "Klebsiella variicola",
    "Klebsiella africana",
}

REPO_DIR = Path(__file__).resolve().parents[2]
ASSETS   = REPO_DIR / "assets"

# ATB CLI data directory (Mac default)
ATB_DATA = Path.home() / "Library" / "Application Support" / "atb" / "data"

ENA_API  = "https://www.ebi.ac.uk/ena/portal/api/filereport"


# ── ATB parquet loader ────────────────────────────────────────────────────────

def load_atb_metadata() -> pd.DataFrame:
    """Read assembly + MLST parquet files exported by ATB CLI."""
    assembly_pq = ATB_DATA / "assembly.parquet"
    mlst_pq     = ATB_DATA / "mlst.parquet"

    if not assembly_pq.exists():
        print(f"ERROR: {assembly_pq} not found. Run: atb fetch --tables assembly.parquet --yes",
              file=sys.stderr)
        sys.exit(1)

    print(f"Loading ATB assembly metadata …", flush=True)
    asm = pd.read_parquet(assembly_pq)
    # Normalise column names to lowercase
    asm.columns = [c.lower().replace(" ", "_") for c in asm.columns]

    print(f"  {len(asm):,} total assemblies", flush=True)

    # Filter to KpSC species
    species_col = next((c for c in asm.columns if "species" in c or "scientific" in c), None)
    hq_col      = next((c for c in asm.columns if "hq" in c or "filter" in c), None)
    run_col     = next((c for c in asm.columns if "run" in c and "accession" in c), None)
    sample_col  = next((c for c in asm.columns if c in ("sample_accession", "name",
                                                          "sample_accession")), None)
    aws_col     = next((c for c in asm.columns if "aws" in c or "s3" in c or "url" in c), None)
    n50_col     = next((c for c in asm.columns if "n50" in c), None)
    comp_col    = next((c for c in asm.columns if "complet" in c), None)
    cont_col    = next((c for c in asm.columns if "contam" in c), None)

    for col, label in [(species_col, "species"), (run_col, "run_accession")]:
        if col is None:
            print(f"ERROR: could not find {label} column. Columns: {list(asm.columns)[:20]}",
                  file=sys.stderr)
            sys.exit(1)

    kpsc_mask = asm[species_col].isin(KPSC_SPECIES)
    if hq_col:
        kpsc_mask = kpsc_mask & (asm[hq_col] == "PASS")
    df = asm[kpsc_mask].copy()
    print(f"  {len(df):,} HQ KpSC assemblies after species + HQ filter", flush=True)

    # Rename to standard names
    rename = {species_col: "species", run_col: "run_accession"}
    if sample_col: rename[sample_col] = "sample_accession"
    if aws_col:    rename[aws_col]    = "aws_url"
    if n50_col:    rename[n50_col]    = "N50"
    if comp_col:   rename[comp_col]   = "completeness"
    if cont_col:   rename[cont_col]   = "contamination"
    df = df.rename(columns=rename)

    # Join MLST
    if mlst_pq.exists():
        mlst = pd.read_parquet(mlst_pq)
        mlst.columns = [c.lower().replace(" ", "_") for c in mlst.columns]
        st_col  = next((c for c in mlst.columns if c in ("mlst_st", "st", "sequence_type")), None)
        sch_col = next((c for c in mlst.columns if "scheme" in c), None)
        sid_col = next((c for c in mlst.columns if "sample" in c and "accession" in c), None)
        if st_col and sid_col:
            # Keep only klebsiella-scheme rows where available
            if sch_col:
                kl = mlst[mlst[sch_col].str.lower().str.contains("klebsiella", na=False)]
                if len(kl) > 0:
                    mlst = kl
            mlst = mlst[[sid_col, st_col]].rename(
                columns={sid_col: "sample_accession", st_col: "ST"}
            ).drop_duplicates("sample_accession")
            if "sample_accession" in df.columns:
                df = df.merge(mlst, on="sample_accession", how="left")
                print(f"  MLST joined: {df['ST'].notna().sum():,} samples have ST",
                      flush=True)
    else:
        df["ST"] = None

    return df


# ── ENA FTP URL builder ───────────────────────────────────────────────────────

def _ena_ftp_url(acc: str) -> tuple[str, str]:
    """
    Compute ENA FTP FASTQ URLs from a run accession.
    Pattern: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{acc[:6]}/{sub}/{acc}/{acc}_{1,2}.fastq.gz
    where sub = '' for 6-digit, '00N' for longer accessions (N = last digit).
    """
    base = "ftp.sra.ebi.ac.uk/vol1/fastq"
    prefix = acc[:6]
    if len(acc) <= 9:
        sub = ""
    else:
        sub = f"00{acc[-1]}"
    path = f"{base}/{prefix}/{sub}/{acc}" if sub else f"{base}/{prefix}/{acc}"
    return f"ftp://{path}/{acc}_1.fastq.gz", f"ftp://{path}/{acc}_2.fastq.gz"


def build_ftp_map(run_accessions: list[str],
                  batch_size: int = 500,
                  use_api: bool = True) -> pd.DataFrame:
    """
    Build run_accession → (ftp_r1, ftp_r2) map.
    Uses ENA Portal API for accuracy; falls back to computed URLs.
    """
    rows = []
    total = len(run_accessions)

    if use_api:
        print(f"Fetching FTP URLs from ENA Portal API ({total:,} accessions, "
              f"batch_size={batch_size}) …", flush=True)
        for i in range(0, total, batch_size):
            batch = run_accessions[i:i + batch_size]
            try:
                r = requests.post(
                    ENA_API,
                    data={
                        "accession": ",".join(batch),
                        "result":    "read_run",
                        "fields":    "run_accession,fastq_ftp",
                        "format":    "tsv",
                        "limit":     0,
                    },
                    timeout=60,
                )
                r.raise_for_status()
                lines = r.text.strip().splitlines()
                if len(lines) < 2:
                    raise ValueError("empty ENA response")
                for line in lines[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        acc = parts[0].strip()
                        ftps = [x.strip() for x in parts[1].split(";") if x.strip()]
                        r1 = next((f for f in ftps if "_1.fastq" in f or f.endswith(".fastq.gz")
                                   and "_2.fastq" not in f), "")
                        r2 = next((f for f in ftps if "_2.fastq" in f), "")
                        rows.append({"run_accession": acc, "ftp_r1": r1, "ftp_r2": r2})
                print(f"  batch {i//batch_size + 1}/{(total-1)//batch_size + 1}: "
                      f"{len(lines)-1} rows", flush=True)
            except Exception as e:
                print(f"  WARNING: ENA API failed for batch {i//batch_size+1}: {e}; "
                      f"using computed URLs", file=sys.stderr)
                for acc in batch:
                    r1, r2 = _ena_ftp_url(acc)
                    rows.append({"run_accession": acc, "ftp_r1": r1, "ftp_r2": r2})
    else:
        print(f"Computing FTP URLs for {total:,} accessions …", flush=True)
        for acc in run_accessions:
            r1, r2 = _ena_ftp_url(acc)
            rows.append({"run_accession": acc, "ftp_r1": r1, "ftp_r2": r2})

    return pd.DataFrame(rows)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate download lists for ATB KpSC cohort expansion"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Print summary without writing files")
    parser.add_argument("--max-samples", type=int, default=None,
                        help="Cap total new samples (default: all)")
    parser.add_argument("--skip-ftp-map", action="store_true",
                        help="Skip ENA FTP map generation (faster; use if only assemblies needed)")
    parser.add_argument("--out-dir", default=str(ASSETS),
                        help=f"Output directory (default: {ASSETS})")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load ATB metadata ────────────────────────────────────────────────────
    atb = load_atb_metadata()

    # ── Load existing cohort ─────────────────────────────────────────────────
    existing_meta = ASSETS / "kpsc_sample_metadata.tsv"
    if existing_meta.exists():
        existing = pd.read_csv(existing_meta, sep="\t")
        existing_ids = set(existing["sample_id"])
        print(f"\nExisting cohort: {len(existing_ids):,} samples", flush=True)
    else:
        existing_ids = set()
        print("No existing cohort found — treating all ATB samples as new", flush=True)

    # ── Identify new samples ─────────────────────────────────────────────────
    new = atb[~atb["run_accession"].isin(existing_ids)].copy()
    print(f"New samples available: {len(new):,}", flush=True)

    # Species breakdown
    print("\nNew samples by species:")
    for sp, n in new["species"].value_counts().items():
        print(f"  {sp:<50} {n:>7,}")

    # ST breakdown for K. pneumoniae
    if "ST" in new.columns:
        kp_new = new[new["species"] == "Klebsiella pneumoniae"]
        print(f"\nTop 20 STs in new K. pneumoniae ({len(kp_new):,} samples):")
        for st, n in kp_new["ST"].astype(str).value_counts().head(20).items():
            print(f"  ST{st}: {n:,}")

    # ── Optional cap ─────────────────────────────────────────────────────────
    if args.max_samples and len(new) > args.max_samples:
        # Stratified sample: proportional to species
        new = new.groupby("species", group_keys=False).apply(
            lambda g: g.sample(
                min(len(g), max(1, int(args.max_samples * len(g) / len(new)))),
                random_state=42,
            )
        ).head(args.max_samples)
        print(f"\nSubsampled to {len(new):,} samples (--max-samples {args.max_samples})")

    if args.dry_run:
        print("\n[dry-run] would write:")
        print(f"  {out_dir}/kpsc_expansion_metadata.tsv        ({len(new):,} rows)")
        print(f"  {out_dir}/kpsc_expansion_sra_accessions.txt  ({len(new):,} accessions)")
        print(f"  {out_dir}/kpsc_expansion_assembly_urls.tsv   ({len(new):,} rows)")
        if not args.skip_ftp_map:
            print(f"  {out_dir}/kpsc_expansion_run_ftp_map.tsv   ({len(new):,} rows)")
        return

    # ── Write metadata TSV ───────────────────────────────────────────────────
    meta_cols = ["run_accession", "species"]
    for c in ["ST", "sample_accession", "N50", "completeness", "contamination"]:
        if c in new.columns:
            meta_cols.append(c)
    meta_out = new[meta_cols].rename(columns={"run_accession": "sample_id"})
    meta_path = out_dir / "kpsc_expansion_metadata.tsv"
    meta_out.to_csv(meta_path, sep="\t", index=False)
    print(f"\nMetadata → {meta_path}  ({len(meta_out):,} rows)", flush=True)

    # ── Write SRA accession list ─────────────────────────────────────────────
    acc_path = out_dir / "kpsc_expansion_sra_accessions.txt"
    new["run_accession"].to_csv(acc_path, index=False, header=False)
    print(f"SRA accessions → {acc_path}", flush=True)

    # ── Write assembly URL list ──────────────────────────────────────────────
    if "aws_url" in new.columns:
        url_cols = ["run_accession"]
        if "sample_accession" in new.columns:
            url_cols.insert(0, "sample_accession")
        url_cols.append("aws_url")
        url_out = new[url_cols]
        url_path = out_dir / "kpsc_expansion_assembly_urls.tsv"
        url_out.to_csv(url_path, sep="\t", index=False)
        print(f"Assembly S3 URLs → {url_path}", flush=True)
    else:
        print("WARNING: no aws_url column found — assembly URL list not written",
              file=sys.stderr)

    # ── Build and write FTP map ──────────────────────────────────────────────
    if not args.skip_ftp_map:
        run_accs = new["run_accession"].dropna().tolist()
        ftp_map  = build_ftp_map(run_accs)
        ftp_path = out_dir / "kpsc_expansion_run_ftp_map.tsv"
        ftp_map.to_csv(ftp_path, sep="\t", index=False)
        print(f"FTP map → {ftp_path}  ({len(ftp_map):,} rows)", flush=True)

    print("\nDone. Transfer to HPC:")
    print(f"  rsync {out_dir}/kpsc_expansion_* pryor:CNVRock/assets/")
    print("\nHPC next steps:")
    print("  N=$(wc -l < assets/kpsc_expansion_sra_accessions.txt)")
    print("  sbatch --array=1-${N}%50 hpc/download_sra.sh  "
          "# uses kpsc_expansion_sra_accessions.txt")
    print("  sbatch --array=1-${N}%50 hpc/download_assemblies_s3.sh  "
          "# for Kleborate")


if __name__ == "__main__":
    main()
