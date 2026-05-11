#!/usr/bin/env python3
"""
Join CNVRock plasmid gene calls (PCN + binary call) with CABBAGE AMR phenotypes.

CABBAGE links BioSample IDs to AMR phenotypes (MIC, R/S/I).
Our cohort uses run accessions (DRR/ERR/SRR); the expansion metadata
provides the bridge: run_accession → sample_accession (BioSample).

Usage
-----
# Expansion cohort (needs kpsc_expansion_metadata.tsv for BioSample bridge):
python3 join_pcn_phenotypes.py \
    --calls  data/results/29_kpsc_phase_c_v3/plasmid_gene_calls.tsv \
    --meta   assets/kpsc_expansion_metadata.tsv \
    --cabbage assets/cabbage_kpsc_phenotypes.tsv \
    --out    data/results/pcn_phenotype_joined.tsv

Outputs
-------
data/results/pcn_phenotype_joined.tsv
  One row per (sample, gene, antibiotic) combination.
  Columns: run_acc, biosample, gene, pcn, call, antibiotic_name,
           priority_antibiotic, ast_standard, measurement,
           measurement_sign, measurement_units, updated_phenotype,
           collection_year, country, isolation_source_category

data/results/pcn_phenotype_summary.tsv
  Per (gene × antibiotic) summary: n_R, n_S, n_I,
  pcn_median_R, pcn_median_S, Mann-Whitney U p-value (R vs S).
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

# ── Gene → antibiotic relevance map ──────────────────────────────────────────
# Used to filter to biologically meaningful pairs (avoids noise from unrelated
# gene–antibiotic combinations in the output).
GENE_ANTIBIOTIC = {
    "blaKPC-2":     ["meropenem", "imipenem", "ertapenem", "doripenem",
                     "ceftazidime-avibactam"],
    "blaNDM-1":     ["meropenem", "imipenem", "ertapenem", "doripenem"],
    "blaOXA-48":    ["meropenem", "imipenem", "ertapenem", "doripenem"],
    "blaCTX-M-15":  ["cefotaxime", "ceftriaxone", "ceftazidime", "cefepime",
                     "aztreonam"],
    "blaCTX-M-14":  ["cefotaxime", "ceftriaxone", "ceftazidime", "cefepime"],
    "qnrB1":        ["ciprofloxacin", "levofloxacin"],
    "aac6-Ib-cr":   ["amikacin", "tobramycin", "gentamicin", "ciprofloxacin"],
    "blaTEM-1":     ["amoxicillin-clavulanic acid", "piperacillin-tazobactam"],
    "blaSHV":       ["cefotaxime", "ceftriaxone", "ceftazidime"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls",   required=True,
                   help="plasmid_gene_calls.tsv from CNVRock evaluation")
    p.add_argument("--meta",    required=True,
                   help="kpsc_expansion_metadata.tsv (run_acc → BioSample bridge)")
    p.add_argument("--cabbage", required=True,
                   help="cabbage_kpsc_phenotypes.tsv")
    p.add_argument("--out",     default="data/results/pcn_phenotype_joined.tsv",
                   help="Output joined table [%(default)s]")
    p.add_argument("--all-pairs", action="store_true",
                   help="Include all gene–antibiotic pairs (not just relevant ones)")
    p.add_argument("--priority-only", action="store_true", default=True,
                   help="Restrict to CABBAGE priority antibiotics [default: True]")
    return p.parse_args()


def load_calls(path: str) -> pd.DataFrame:
    """Load gene calls TSV; melt to long format with (sample, gene, pcn, call)."""
    df = pd.read_csv(path, sep="\t")
    pcn_cols  = [c for c in df.columns if c.startswith("pcn_")]
    call_cols = [c for c in df.columns if not c.startswith("pcn_") and c != "sample_id"]

    # Melt calls
    calls_long = df[["sample_id"] + call_cols].melt(
        id_vars="sample_id", var_name="gene", value_name="call"
    )
    # Melt PCN (strip "pcn_" prefix to match gene names)
    pcn_long = df[["sample_id"] + pcn_cols].melt(
        id_vars="sample_id", var_name="gene", value_name="pcn"
    )
    pcn_long["gene"] = pcn_long["gene"].str.removeprefix("pcn_")

    merged = calls_long.merge(pcn_long, on=["sample_id", "gene"])
    return merged


def main():
    args = parse_args()

    # ── Load inputs ──────────────────────────────────────────────────────────
    print("Loading gene calls …", flush=True)
    calls = load_calls(args.calls)
    print(f"  {calls.sample_id.nunique()} samples, {calls.gene.nunique()} genes, "
          f"{len(calls):,} rows")

    print("Loading expansion metadata (run → BioSample bridge) …", flush=True)
    meta = pd.read_csv(args.meta, sep="\t")[["sample_id", "sample_accession"]] \
             .rename(columns={"sample_id": "run_acc",
                              "sample_accession": "biosample"})
    print(f"  {len(meta):,} samples in metadata")

    print("Loading CABBAGE phenotypes …", flush=True)
    cab = pd.read_csv(args.cabbage, sep="\t")
    if args.priority_only:
        cab = cab[cab.priority_antibiotic]
        print(f"  {len(cab):,} rows after filtering to priority antibiotics")
    print(f"  {cab.sample_id.nunique()} unique CABBAGE samples, "
          f"{cab.antibiotic_name.nunique()} antibiotics")

    # ── Bridge run accessions → BioSample ───────────────────────────────────
    # calls.sample_id = run accession (DRR/ERR/SRR)
    # meta bridges run_acc → biosample (SAMN/SAMD/SAME)
    # CABBAGE BioSample_ID uses the same BioSample namespace
    calls_with_bio = calls.merge(meta, left_on="sample_id", right_on="run_acc",
                                  how="inner")
    n_bridged = calls_with_bio.run_acc.nunique()
    print(f"\nRun accessions with BioSample bridge: {n_bridged} "
          f"(of {calls.sample_id.nunique()} in calls)")

    # ── Join to CABBAGE ──────────────────────────────────────────────────────
    joined = calls_with_bio.merge(
        cab, left_on="biosample", right_on="BioSample_ID", how="inner"
    )
    print(f"After CABBAGE join: {joined.run_acc.nunique()} samples, "
          f"{len(joined):,} rows")

    # ── Filter to relevant gene–antibiotic pairs ─────────────────────────────
    if not args.all_pairs:
        relevant_rows = []
        for gene, abxs in GENE_ANTIBIOTIC.items():
            mask = (joined.gene == gene) & (joined.antibiotic_name.isin(abxs))
            relevant_rows.append(joined[mask])
        joined = pd.concat(relevant_rows, ignore_index=True)
        print(f"After gene–antibiotic relevance filter: {len(joined):,} rows, "
              f"{joined.run_acc.nunique()} samples")

    # ── Select and rename output columns ────────────────────────────────────
    keep = ["run_acc", "biosample", "gene", "pcn", "call",
            "antibiotic_name", "priority_antibiotic", "ast_standard",
            "measurement", "measurement_sign", "measurement_units",
            "updated_phenotype", "collection_year",
            "ISO_country_code", "country",
            "isolation_source_category"]
    keep = [c for c in keep if c in joined.columns]
    out_df = joined[keep].copy()
    out_df = out_df.sort_values(["gene", "antibiotic_name", "run_acc"])

    # ── Save joined table ────────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nJoined table → {out_path}  ({len(out_df):,} rows)")

    # ── Summary table: per (gene × antibiotic) stats ─────────────────────────
    summary_rows = []
    for (gene, abx), grp in out_df.groupby(["gene", "antibiotic_name"]):
        r = grp[grp.updated_phenotype == "resistant"]["pcn"]
        s = grp[grp.updated_phenotype == "susceptible"]["pcn"]
        i = grp[grp.updated_phenotype == "intermediate"]["pcn"]
        n_total = len(grp.run_acc.unique())

        if len(r) >= 5 and len(s) >= 5:
            u_stat, p_val = stats.mannwhitneyu(r, s, alternative="greater")
        else:
            u_stat, p_val = np.nan, np.nan

        summary_rows.append({
            "gene":          gene,
            "antibiotic":    abx,
            "n_samples":     n_total,
            "n_R":           len(r),
            "n_S":           len(s),
            "n_I":           len(i),
            "pcn_median_R":  round(r.median(), 3) if len(r) else np.nan,
            "pcn_median_S":  round(s.median(), 3) if len(s) else np.nan,
            "mwu_p_R_gt_S":  round(p_val, 4) if not np.isnan(p_val) else np.nan,
        })

    if not summary_rows:
        print("No data for summary (0 rows after join — check that --meta covers "
              "the same samples as --calls).")
        return

    summary = pd.DataFrame(summary_rows).sort_values(["gene", "antibiotic"])
    summary_path = out_path.parent / (out_path.stem.replace("joined", "summary") + ".tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Summary table  → {summary_path}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
