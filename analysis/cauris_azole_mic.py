#!/usr/bin/env python3
"""
Fetch azole MIC data for C. auris isolates in our cohort.

Sources:
  1. NCBI BioSample antibiogram attributes (fluconazole/voriconazole)
  2. Supplementary tables from key publications (Kappel et al. 2024,
     PRJEB36822; Rybak/TAC1B 2020; Colombia PRJNA1003896)
  3. ENA sample attributes (isolation_source, host, geo_loc, phenotype)

ERG11 mutation status is separate (see cauris_mutation_gt.py).
This script builds the MIC table so the cross-tab
  ERG11 amplification (CNVRock CRR) × ERG11 mutation × azole MIC
can be computed.

Run on HPC:
  python3 analysis/cauris_azole_mic.py

Output:
  data/results/cauris_mutation_gt/cauris_azole_mic.tsv
  data/results/cauris_mutation_gt/cauris_mic_bioprojects.tsv
"""
from __future__ import annotations
import io, json, time
from pathlib import Path
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import pandas as pd

REPO = Path("/home/lshlt19/CNVRock")
MANIFEST = REPO / "assets/cauris_subset_500_mut.tsv"
OUT_DIR  = REPO / "data/results/cauris_mutation_gt"
OUT_MIC  = OUT_DIR / "cauris_azole_mic.tsv"
OUT_META = OUT_DIR / "cauris_bioproject_meta.tsv"

# BioProjects known to have ERG11 mutations + azole MIC data (from literature)
# Reference: Kappel/Rhodes 2024 (npj, PRJEB36822), Rybak 2020 (PRJNA493622),
#            Colombia 2024 (PRJNA1003896), original Lockhart (PRJNA328792)
TARGET_BIOPROJECTS = [
    "PRJEB36822",   # UK 207 isolates, Kappel et al. 2024 — all clades, Y132F/K143R/F126L
    "PRJNA493622",  # USA 318 isolates, multiple introductions — Clades I/III/IV
    "PRJNA1003896", # Colombia 99 isolates, Clade IV, Y132F
    "PRJNA328792",  # Original Lockhart 2017 — all 4 clades
    "PRJEB20230",   # UK Rhodes 2018 — 27 Clade I, all Y132F
]

AZOLE_ATTRS = ["fluconazole", "voriconazole", "anidulafungin",
               "caspofungin", "micafungin", "amphotericin"]


def ncbi_biosamples_for_project(bioproject: str, retmax: int = 500) -> list[str]:
    """Return NCBI BioSample UIDs for a BioProject accession."""
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=biosample&term={bioproject}[BioProject]"
           f"&retmax={retmax}&retmode=json")
    try:
        resp = json.loads(urlopen(url, timeout=30).read())
        return resp["esearchresult"]["idlist"]
    except Exception as e:
        print(f"  esearch {bioproject}: {e}")
        return []


def ncbi_biosample_attrs(uids: list[str]) -> list[dict]:
    """Fetch BioSample XML for a list of UIDs, return attribute dicts."""
    if not uids:
        return []
    chunk = 100
    results = []
    for i in range(0, len(uids), chunk):
        batch = uids[i:i+chunk]
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               f"?db=biosample&id={','.join(batch)}&rettype=xml")
        try:
            xml_bytes = urlopen(url, timeout=60).read()
            root = ET.fromstring(xml_bytes)
            for bs in root.findall(".//BioSample"):
                acc = bs.get("accession", "?")
                attrs = {a.get("attribute_name","").lower(): a.text
                         for a in bs.findall(".//Attribute")}
                results.append({"biosample": acc, **attrs})
            time.sleep(0.4)
        except Exception as e:
            print(f"  efetch batch {i}: {e}")
    return results


def ena_run_table(bioproject: str) -> pd.DataFrame | None:
    """Fetch ENA run table for a project (accession + sample metadata)."""
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={bioproject}&result=read_run"
           f"&fields=run_accession,sample_accession,study_accession,"
           f"country,collection_date,host,isolation_source"
           f"&format=tsv&limit=0")
    try:
        raw = urlopen(url, timeout=60).read()
        df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str)
        print(f"  ENA {bioproject}: {len(df)} runs")
        return df
    except Exception as e:
        print(f"  ENA {bioproject}: {e}")
        return None


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mf = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    our_bs  = set(mf["biosample"].dropna())
    our_acc = set(mf["accession"].dropna())
    print(f"Our 500-sample manifest: {len(mf)} samples")

    all_rows: list[dict] = []
    proj_meta: list[dict] = []

    for proj in TARGET_BIOPROJECTS:
        print(f"\n── {proj} ──────────────────────────────")

        # 1. ENA run table (fast, gives run→sample mapping + geo)
        ena = ena_run_table(proj)
        if ena is not None:
            # check overlap with our manifest
            overlap_acc = set(ena.get("run_accession", pd.Series()).dropna()) & our_acc
            overlap_bs  = set(ena.get("sample_accession", pd.Series()).dropna()) & our_bs
            print(f"  Overlap with our manifest: {len(overlap_acc)} runs, "
                  f"{len(overlap_bs)} biosamples")
            proj_meta.append({"bioproject": proj, "n_runs": len(ena),
                               "overlap_acc": len(overlap_acc),
                               "overlap_bs": len(overlap_bs)})

        # 2. NCBI BioSample attributes (MIC data if available)
        uids = ncbi_biosamples_for_project(proj)
        print(f"  NCBI UIDs: {len(uids)}")
        if uids:
            attrs_list = ncbi_biosample_attrs(uids)
            for attrs in attrs_list:
                mic_data = {k: v for k, v in attrs.items()
                            if any(a in k for a in AZOLE_ATTRS)}
                if mic_data or attrs.get("biosample") in our_bs:
                    row = {"biosample": attrs.get("biosample","?"),
                           "bioproject": proj}
                    row.update(mic_data)
                    row["in_our_cohort"] = int(attrs.get("biosample","?") in our_bs)
                    all_rows.append(row)
            n_mic = sum(1 for r in attrs_list
                        if any(a in k for k in r for a in AZOLE_ATTRS))
            print(f"  BioSamples with any azole MIC: {n_mic}")

    # Save results
    mic_df = pd.DataFrame(all_rows)
    mic_df.to_csv(OUT_MIC, sep="\t", index=False)
    print(f"\nWrote {OUT_MIC} ({len(mic_df)} rows)")

    meta_df = pd.DataFrame(proj_meta)
    meta_df.to_csv(OUT_META, sep="\t", index=False)
    print(f"Wrote {OUT_META}")

    # Summary
    if not mic_df.empty:
        in_cohort = mic_df[mic_df.get("in_our_cohort",pd.Series(0)) == 1]
        print(f"\nSamples in our cohort with MIC data: {len(in_cohort)}")
        print(f"Total overlapping samples found: {len(mic_df)}")

        for col in [c for c in mic_df.columns if any(a in c for a in AZOLE_ATTRS)]:
            n = mic_df[col].notna().sum()
            if n > 0:
                print(f"  {col}: {n} values")


if __name__ == "__main__":
    main()
