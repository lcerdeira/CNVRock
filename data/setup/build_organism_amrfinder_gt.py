#!/usr/bin/env python3
"""
Build AMRFinder presence ground-truth for an arbitrary organism cohort.

Generalises build_abaumannii_amrfinder_gt.py. Reads the AllTheBacteria
pre-computed AMRFinderPlus parquet, restricts to the cohort BioSamples, and
emits one binary presence column per evaluable gene.

Only *acquired* / variably-present genes are written as GT. Intrinsic
chromosomal genes (norA in S. aureus; the ampC/acrAB efflux system in
E. coli) are present in essentially every strain, so AMRFinderPlus does not
catalogue them as presence calls — the parquet scan confirms this directly
(norA: 0 rows; acrB: 820 across 2.6M samples). Those genes have no binary
GT and are evaluated by their copy-ratio distribution, exactly as the
intrinsic efflux genes were for A. baumannii.

    python3 data/setup/build_organism_amrfinder_gt.py --organism saureus
    python3 data/setup/build_organism_amrfinder_gt.py --organism ecoli
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ATB = Path(os.path.expanduser(
    "~/Library/Application Support/atb/data/amrfinderplus.parquet"))

# Per organism: the manifest, the output path, and the evaluable-gene
# predicates. Each predicate maps (element_symbol, subclass) -> bool. Column
# names must match the call_ids in the organism's gene_coords.tsv.
ORGANISMS = {
    "saureus": dict(
        manifest="assets/saureus_subset.tsv",
        out="assets/amrfinder_gt_saureus.tsv",
        genes={
            "mecA": lambda s, sub: s == "mecA",
            "blaZ": lambda s, sub: s == "blaZ",
            "ermC": lambda s, sub: s == "erm(C)",
            "ermB": lambda s, sub: s == "erm(B)",
            "tetK": lambda s, sub: s == "tet(K)",
            "mupA": lambda s, sub: s in ("mupA", "mupA-like"),
            # norA is intrinsic -> not a presence GT column (see module docstring)
        },
    ),
    "ecoli": dict(
        manifest="assets/ecoli_subset.tsv",
        out="assets/amrfinder_gt_ecoli.tsv",
        genes={
            "blaCTX-M":      lambda s, sub: s.startswith("blaCTX-M"),
            "blaTEM":        lambda s, sub: s.startswith("blaTEM"),
            "blaSHV":        lambda s, sub: s.startswith("blaSHV"),
            "blaKPC":        lambda s, sub: s.startswith("blaKPC"),
            "blaNDM":        lambda s, sub: s.startswith("blaNDM"),
            "blaOXA-48-like": lambda s, sub: any(
                s.startswith(p) for p in
                ("blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-244")),
            "qnrB":          lambda s, sub: s.startswith("qnrB"),
            "aac6-Ib-cr":    lambda s, sub: s.replace("'", "").startswith("aac(6)-Ib"),
            "sul1":          lambda s, sub: s == "sul1",
            "sul2":          lambda s, sub: s == "sul2",
            "dfrA12":        lambda s, sub: s == "dfrA12",
            "dfrA14":        lambda s, sub: s == "dfrA14",
            "mcr-1":         lambda s, sub: s.startswith("mcr-1"),
            # ampC / acrB are intrinsic -> evaluated by CRR, not presence
        },
    ),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, choices=sorted(ORGANISMS))
    args = ap.parse_args()
    cfg = ORGANISMS[args.organism]
    out = REPO / cfg["out"]

    if not ATB.exists():
        raise SystemExit(f"ATB parquet not found: {ATB}")

    manifest = pd.read_csv(REPO / cfg["manifest"], sep="\t", dtype=str)
    bs_to_acc = dict(zip(manifest["biosample"], manifest["accession"]))
    bs_set = set(bs_to_acc)
    print(f"{args.organism} cohort: {len(bs_set):,} biosamples")

    print("Loading ATB AMRFinderPlus parquet …")
    amr = pd.read_parquet(ATB, columns=["Name", "Element symbol", "Subclass"])
    sub = amr[amr["Name"].isin(bs_set)].copy()
    sub["Element symbol"] = sub["Element symbol"].astype(str)
    sub["Subclass"] = sub["Subclass"].astype(str)
    print(f"ATB rows for cohort: {len(sub):,} across "
          f"{sub['Name'].nunique():,} biosamples")

    results: dict[str, object] = {"biosample": sorted(bs_set)}
    bs_index = pd.Series(results["biosample"])
    for gene, pred in cfg["genes"].items():
        mask = sub.apply(lambda r: pred(r["Element symbol"], r["Subclass"]),
                         axis=1)
        hits = sub[mask].groupby("Name").size().gt(0).astype(int)
        results[gene] = bs_index.map(hits).fillna(0).astype(int).values

    df = pd.DataFrame(results)
    df.insert(0, "sample_id", df["biosample"].map(bs_to_acc))
    n_no_rows = len(bs_set) - sub["Name"].nunique()
    if n_no_rows:
        print(f"  {n_no_rows} biosamples had no AMRFinder rows (all genes -> 0)")

    df = df.sort_values("sample_id").reset_index(drop=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(df):,} rows)")
    print("\nGene prevalence:")
    for gene in cfg["genes"]:
        n = int(df[gene].sum())
        print(f"  {gene:<16} {n:>6} / {len(df)}  ({100*n/len(df):.1f} %)")


if __name__ == "__main__":
    main()
