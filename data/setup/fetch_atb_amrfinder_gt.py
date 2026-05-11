"""
fetch_atb_amrfinder_gt.py
─────────────────────────
Pull pre-computed AMRFinderPlus gene calls from the AllTheBacteria (ATB)
parquet for:
  • The 88,135 Phase 2 expansion samples  (run accession → BioSample via
    kpsc_expansion_metadata.tsv)
  • The 545 Phase 1 training samples      (run accession → BioSample via
    ATB assembly.parquet)

Output
------
assets/amrfinder_gt_expansion.tsv   — expansion cohort GT (sample_id + genes)
assets/amrfinder_gt_training.tsv    — training cohort GT  (sample_id + genes)

Columns (both files):
  sample_id        — run accession (DRR*/ERR*/SRR*)
  biosample        — BioSample accession (SAMD*/SAMN*/SAME*)
  blaSHV           — 1/0 any blaSHV variant present
  blaKPC           — 1/0 any blaKPC variant present
  blaNDM           — 1/0 any blaNDM variant present
  blaCTX-M         — 1/0 any blaCTX-M variant present
  qnrB             — 1/0 any qnrB variant present
  blaOXA-48-like   — 1/0 any OXA-48-group variant (Subclass==CARBAPENEM)
  aac6-Ib-cr       — 1/0 aac(6')-Ib-cr or aac(6')-Ib-cr5 present

Gene-detection logic uses `Element symbol` from AMRFinderPlus; OXA-48-like
is additionally filtered by Subclass==CARBAPENEM to exclude chromosomal OXA
genes (OXA-1, OXA-9, etc.).
"""

import os
import pandas as pd

ATB_DIR   = os.path.expanduser("~/Library/Application Support/atb/data/")
REPO_DIR  = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ASSETS    = os.path.join(REPO_DIR, "assets")

EXP_META  = os.path.join(ASSETS, "kpsc_expansion_metadata.tsv")
TRAIN_META = os.path.join(ASSETS, "kpsc_sample_metadata.tsv")

OUT_EXP   = os.path.join(ASSETS, "amrfinder_gt_expansion.tsv")
OUT_TRAIN = os.path.join(ASSETS, "amrfinder_gt_training.tsv")

# ── Gene patterns (applied to 'Element symbol') ───────────────────────────────
GENE_PATTERNS = {
    "blaSHV":        r"^blaSHV",
    "blaKPC":        r"^blaKPC",
    "blaNDM":        r"^blaNDM",
    "blaCTX-M":      r"^blaCTX-M",
    "qnrB":          r"^qnrB",
    # OXA-48-like: additionally filtered by Subclass==CARBAPENEM below
    "blaOXA-48-like": r"^blaOXA",
    # -cr variants confer fluoroquinolone resistance; include cr5 (same phenotype)
    "aac6-Ib-cr":    r"^aac\(6'\)-Ib-cr",
}


def load_amrfinder_for_samples(biosample_set: set) -> pd.DataFrame:
    """
    Load AMRFinderPlus rows for a given set of BioSample accessions.
    Returns a DataFrame with columns: Name, Element symbol, Subclass.
    """
    print(f"  Loading AMRFinderPlus parquet for {len(biosample_set):,} BioSamples …")
    amr = pd.read_parquet(
        os.path.join(ATB_DIR, "amrfinderplus.parquet"),
        columns=["Name", "Element symbol", "Subclass"],
    )
    amr = amr[amr["Name"].isin(biosample_set)].copy()
    print(f"  → {len(amr):,} AMRFinder rows for {amr['Name'].nunique():,} samples")
    return amr


def pivot_to_presence(amr: pd.DataFrame) -> pd.DataFrame:
    """
    Convert long AMRFinder rows to a wide BioSample × gene presence matrix.
    Each gene column is 1 if any matching row exists for that sample, else 0.
    """
    results = {"biosample": sorted(amr["Name"].unique())}

    for gene, pattern in GENE_PATTERNS.items():
        mask = amr["Element symbol"].str.match(pattern, na=False)
        if gene == "blaOXA-48-like":
            mask &= (amr["Subclass"] == "CARBAPENEM")
        hits = amr[mask].groupby("Name").size().gt(0).astype(int)
        results[gene] = (
            pd.Series(results["biosample"])
            .map(hits)
            .fillna(0)
            .astype(int)
            .values
        )

    return pd.DataFrame(results).rename(columns={"biosample": "biosample"})


def process_cohort(
    run_to_biosample: dict,
    label: str,
    out_path: str,
) -> pd.DataFrame:
    """
    For a cohort (dict: run_accession → biosample), fetch AMRFinder GT and
    write TSV.  Returns the final DataFrame.
    """
    print(f"\n{'='*60}")
    print(f"Processing {label} ({len(run_to_biosample):,} samples)")
    print(f"{'='*60}")

    biosample_set = set(run_to_biosample.values())
    amr = load_amrfinder_for_samples(biosample_set)

    missing = biosample_set - set(amr["Name"].unique())
    if missing:
        print(f"  ⚠  {len(missing)} BioSamples not in AMRFinder (no gene calls):")
        for s in sorted(missing)[:10]:
            print(f"       {s}")
        if len(missing) > 10:
            print(f"       … and {len(missing)-10} more")

    wide = pivot_to_presence(amr)

    # Add run accession (sample_id) by reversing the mapping
    biosample_to_run = {v: k for k, v in run_to_biosample.items()}
    wide.insert(0, "sample_id", wide["biosample"].map(biosample_to_run))

    # Samples with no AMRFinder rows → add rows of 0s
    covered = set(wide["biosample"])
    missing_rows = [
        {"sample_id": run, "biosample": bs, **{g: 0 for g in GENE_PATTERNS}}
        for bs, run in biosample_to_run.items()
        if bs not in covered
    ]
    if missing_rows:
        wide = pd.concat([wide, pd.DataFrame(missing_rows)], ignore_index=True)

    wide = wide.sort_values("sample_id").reset_index(drop=True)
    wide.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Written: {out_path}  ({len(wide):,} rows)")

    # Summary
    print("\n  Gene prevalence:")
    for gene in GENE_PATTERNS:
        n = wide[gene].sum()
        pct = 100 * n / len(wide)
        print(f"    {gene:<20s}: {n:6,} / {len(wide):,}  ({pct:.1f}%)")

    return wide


def main():
    # ── Expansion cohort ──────────────────────────────────────────────────────
    exp_meta = pd.read_csv(EXP_META, sep="\t")
    exp_map  = dict(zip(exp_meta["sample_id"], exp_meta["sample_accession"]))
    process_cohort(exp_map, "Phase 2 expansion (88k)", OUT_EXP)

    # ── Training cohort ───────────────────────────────────────────────────────
    train_meta = pd.read_csv(TRAIN_META, sep="\t")
    asm = pd.read_parquet(
        os.path.join(ATB_DIR, "assembly.parquet"),
        columns=["sample_accession", "run_accession"],
    )
    asm["run_accession"] = asm["run_accession"].str.split(",").str[0].str.strip()
    asm = asm[asm["run_accession"].isin(train_meta["sample_id"])].drop_duplicates("run_accession")
    train_map = dict(zip(asm["run_accession"], asm["sample_accession"]))
    process_cohort(train_map, "Phase 1 training (545)", OUT_TRAIN)

    print("\n✓  Done.")


if __name__ == "__main__":
    main()
