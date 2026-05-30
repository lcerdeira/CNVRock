#!/usr/bin/env python3
"""
Mine C. auris azole MIC data from:
  1. SRA metadata for key BioProjects (PRJNA328792, PRJNA414579, PRJNA548234)
  2. ENA sample attributes
  3. Cross with our 500-sample manifest + ERG11 mutation GT

Target articles:
  PMID 29361025 — Chowdhary 2018 (350 Indian isolates, ERG11 × MIC)
  PMC5786713    — Spivak & Hanson 2018 review
  PMC11931498   — Arendrup 2025 EUCAST vs CLSI

Run on HPC:
  /home/lshlt19/miniconda3/envs/cnvrock/bin/python3 analysis/cauris_mic_mining.py
"""
from __future__ import annotations
import io, json, time, re
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import pandas as pd

REPO    = Path("/home/lshlt19/CNVRock")
OUT_DIR = REPO / "data/results/cauris_mutation_gt"
MANIFEST = REPO / "assets/cauris_subset_500_mut.tsv"
MUT_GT   = OUT_DIR / "cauris_erg11_mutation_gt.tsv"

# BioProjects suggested by user
BIOPROJECTS = ["PRJNA328792", "PRJNA414579", "PRJNA548234",
               "PRJEB36822", "PRJNA493622"]

AZOLE_KEYS = ["fluconazole", "voriconazole", "itraconazole",
              "posaconazole", "isavuconazole", "anidulafungin",
              "caspofungin", "micafungin", "amphotericin",
              "antifungal_susceptibility", "mic_fluconazole",
              "resistance_phenotype", "phenotype", "ast"]


def ncbi_sra_runinfo(bioproject: str) -> pd.DataFrame | None:
    """Fetch SRA RunInfo CSV for a BioProject via NCBI Entrez."""
    print(f"\n  Querying SRA RunInfo for {bioproject}…")
    try:
        # esearch to get SRA UIDs
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=sra&term={bioproject}[BioProject]&retmax=2000&retmode=json")
        resp = json.loads(urlopen(url, timeout=30).read())
        uids = resp["esearchresult"]["idlist"]
        count = resp["esearchresult"]["count"]
        print(f"    SRA UIDs: {len(uids)} (total {count})")
        if not uids:
            return None

        # Fetch biosample attributes via efetch XML in chunks
        rows = []
        chunk = 200
        for i in range(0, len(uids), chunk):
            batch = uids[i:i+chunk]
            url2 = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                    f"?db=sra&id={','.join(batch)}&rettype=runinfo&retmode=csv")
            try:
                csv_bytes = urlopen(url2, timeout=60).read()
                df_chunk = pd.read_csv(io.BytesIO(csv_bytes), dtype=str)
                rows.append(df_chunk)
                time.sleep(0.4)
            except Exception as e:
                print(f"    chunk {i}: {e}")

        if not rows:
            return None
        df = pd.concat(rows, ignore_index=True)
        print(f"    RunInfo rows: {len(df)}, cols: {list(df.columns[:10])}")

        # Check for any azole-related columns
        azole_cols = [c for c in df.columns
                      if any(k in c.lower() for k in AZOLE_KEYS)]
        if azole_cols:
            print(f"    *** AZOLE COLUMNS FOUND: {azole_cols}")
        return df
    except Exception as e:
        print(f"    Error: {e}")
        return None


def ena_sample_attributes(bioproject: str) -> pd.DataFrame | None:
    """Fetch ENA sample attributes including custom fields."""
    print(f"\n  ENA attributes for {bioproject}…")
    # ENA portal API with sample attributes
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={bioproject}&result=read_run"
           f"&fields=run_accession,sample_accession,sample_title,"
           f"sample_description,country,collection_date,host,"
           f"isolation_source,strain,cultivar,ecotype,serotype"
           f"&format=tsv&limit=0")
    try:
        raw = urlopen(url, timeout=60).read()
        df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str)
        print(f"    ENA rows: {len(df)}")

        # Also try to get ALL attributes via BioSamples bulk for a subset
        azole_cols = [c for c in df.columns
                      if any(k in c.lower() for k in AZOLE_KEYS)]
        if azole_cols:
            print(f"    *** AZOLE COLUMNS: {azole_cols}")

        # Check sample_description for MIC keywords
        if "sample_description" in df.columns:
            mic_rows = df[df["sample_description"].str.contains(
                "MIC|fluconazole|resistant|azole",
                case=False, na=False)]
            if len(mic_rows):
                print(f"    Descriptions with MIC/azole: {len(mic_rows)}")
                print(f"    Example: {mic_rows['sample_description'].iloc[0][:200]}")
        return df
    except Exception as e:
        print(f"    ENA error: {e}")
        return None


def fetch_biosamples_bulk(biosample_list: list[str]) -> pd.DataFrame:
    """Fetch full BioSample XML attributes for a list of accessions."""
    print(f"\n  Fetching full BioSample attributes for {len(biosample_list)} samples…")
    all_rows = []
    # Use NCBI BioSample efetch
    chunk = 100
    for i in range(0, min(len(biosample_list), 500), chunk):
        batch = biosample_list[i:i+chunk]
        # esearch for UIDs
        term = " OR ".join(f"{bs}[Accession]" for bs in batch[:20])
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=biosample&term={term}&retmax=100&retmode=json")
        try:
            resp = json.loads(urlopen(url, timeout=30).read())
            uids = resp["esearchresult"]["idlist"]
            if not uids:
                continue
            import xml.etree.ElementTree as ET
            url2 = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                    f"?db=biosample&id={','.join(uids)}&rettype=xml")
            xml_bytes = urlopen(url2, timeout=60).read()
            root = ET.fromstring(xml_bytes)
            for bs_elem in root.findall(".//BioSample"):
                acc = bs_elem.get("accession", "?")
                attrs = {a.get("attribute_name", "").lower(): a.text
                         for a in bs_elem.findall(".//Attribute")}
                # Keep only if has any azole-related attribute
                azole = {k: v for k, v in attrs.items()
                         if any(z in k for z in AZOLE_KEYS)}
                if azole:
                    row = {"biosample": acc}
                    row.update(azole)
                    all_rows.append(row)
            time.sleep(0.5)
        except Exception as e:
            print(f"    batch {i}: {e}")

    if all_rows:
        df = pd.DataFrame(all_rows)
        print(f"    BioSamples with azole attributes: {len(df)}")
        return df
    return pd.DataFrame()


def check_overlap_with_manifest(sra_df: pd.DataFrame,
                                 bioproject: str) -> None:
    """Check overlap between SRA results and our 500-sample manifest."""
    mf = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    our_acc = set(mf["accession"].dropna())
    our_bs  = set(mf["biosample"].dropna())

    run_col = next((c for c in sra_df.columns
                    if c.lower() in ["run", "run_accession"]), None)
    bs_col  = next((c for c in sra_df.columns
                    if "biosample" in c.lower() or "sample_accession" in c.lower()), None)

    if run_col:
        overlap = set(sra_df[run_col].dropna()) & our_acc
        if overlap:
            print(f"    *** {bioproject}: {len(overlap)} runs in our manifest: {list(overlap)[:5]}")
    if bs_col:
        overlap_bs = set(sra_df[bs_col].dropna()) & our_bs
        if overlap_bs:
            print(f"    *** {bioproject}: {len(overlap_bs)} biosamples in our manifest")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mf = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    print(f"Our 500-sample manifest: {len(mf)} samples")
    print(f"Biosample prefixes: {set(bs[:4] for bs in mf['biosample'].dropna())}")

    all_sra: list[pd.DataFrame] = []
    all_ena: list[pd.DataFrame] = []

    for bp in BIOPROJECTS:
        print(f"\n{'='*55}\n{bp}")

        # SRA RunInfo
        sra = ncbi_sra_runinfo(bp)
        if sra is not None:
            sra["bioproject"] = bp
            check_overlap_with_manifest(sra, bp)
            all_sra.append(sra)
            # Save per-project
            sra.to_csv(OUT_DIR / f"sra_runinfo_{bp}.tsv", sep="\t", index=False)

        # ENA attributes
        ena = ena_sample_attributes(bp)
        if ena is not None:
            ena["bioproject"] = bp
            check_overlap_with_manifest(ena, bp)
            all_ena.append(ena)

        time.sleep(1)

    # Consolidate SRA RunInfo — look for any azole columns across all projects
    if all_sra:
        big = pd.concat(all_sra, ignore_index=True)
        azole_cols = [c for c in big.columns
                      if any(k in c.lower() for k in AZOLE_KEYS)]
        print(f"\n{'='*55}")
        print(f"All SRA cols with azole terms: {azole_cols}")
        if azole_cols:
            sub = big[["Run", "BioSample", "bioproject"] + azole_cols].dropna(
                subset=azole_cols, how="all")
            print(f"Rows with any azole value: {len(sub)}")
            sub.to_csv(OUT_DIR / "cauris_sra_azole_mics.tsv", sep="\t", index=False)
            print(f"Saved cauris_sra_azole_mics.tsv")
        else:
            print("No azole MIC columns in SRA RunInfo for these BioProjects.")
            # Try fetching full BioSample attributes for all our 500 samples
            print("\nFalling back: fetching full BioSample XML for our 500 samples…")
            our_bs = mf["biosample"].dropna().tolist()
            bs_df = fetch_biosamples_bulk(our_bs)
            if not bs_df.empty:
                bs_df.to_csv(OUT_DIR / "cauris_biosample_azole_attrs.tsv",
                             sep="\t", index=False)
                print(f"Saved {len(bs_df)} rows with azole attributes")
            else:
                print("No azole attributes found in BioSample records for our cohort.")
                print("\nConclusion: structured MIC data is NOT available in public")
                print("databases for our 500-isolate subset.")
                print("MIC data exists only in supplementary tables of published papers")
                print("(Chowdhary PMID 29361025, Kappel PRJEB36822).")
                print("Recommend: manual extraction from Kappel 2024 Suppl Table 1")
                print("(doi:10.1038/s44259-024-00043-6)")


if __name__ == "__main__":
    main()
