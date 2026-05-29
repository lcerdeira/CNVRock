#!/usr/bin/env python3
"""
Long-read validation of CNVRock chromosomal CNV calls (KpSC blaSHV).

Reviewer question: what is the read-depth "ground truth"? Short-read
assembly collapses tandem duplications, so it cannot serve. Long-read /
hybrid assemblies DO resolve them. Of the 10K KpSC cohort, 415 biosamples
also have an ONT/PacBio run and ~131 have an NCBI assembly — an orthogonal
ground truth for the blaSHV-amplification headline.

Method:
  1. For 10K-cohort biosamples that have a long-read run, fetch the NCBI
     assembly (datasets).
  2. BLAST the blaSHV CDS against each assembly; count distinct high-identity
     hits (>=95% identity, >=90% query coverage) = assembly blaSHV copy
     number (a long-read assembly resolves tandem copies short reads miss).
  3. Join to CNVRock's per-sample crr_blaSHV (gene_calls.tsv) via
     biosample -> run accession.
  4. Report: does CNVRock CRR rise with assembly copy number?
     (Spearman; CRR distribution by assembly copy 1 vs >=2.)

Run on HPC (datasets + blast in blast_env / cnvrock):
  python3 analysis/longread_validation.py
"""
from __future__ import annotations
import io, os, subprocess, tempfile, glob
from pathlib import Path
from urllib.request import urlopen
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML  # noqa: F401 (kept for reference)

REPO = Path("/home/lshlt19/CNVRock")
META = REPO / "assets/kpsc_expansion_metadata_runlevel.tsv"
MANIFEST10K = REPO / "assets/kpsc_expansion_subset_10k.tsv"
HS11286 = REPO / "assets/HS11286_extended.fasta"
GENE_CALLS = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT = REPO / "data/results/longread_validation"
BLASTN = "blastn"  # from blast_env on PATH
MAKEBLASTDB = "makeblastdb"
SHV_CONTIG, SHV_START, SHV_END = "NC_016845.1", 2549403, 2550263


def extract_shv_cds():
    for rec in SeqIO.parse(str(HS11286), "fasta"):
        if rec.id == SHV_CONTIG:
            return str(rec.seq[SHV_START - 1:SHV_END])
    raise RuntimeError("blaSHV contig not found")


def overlap_biosamples():
    meta = pd.read_csv(META, sep="\t", dtype=str)
    s10k = set(pd.read_csv(MANIFEST10K, sep="\t", dtype=str)["accession"])
    our = meta[meta["sample_id"].isin(s10k)][["sample_id", "sample_accession"]].dropna()
    url = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run"
           "&query=tax_eq(573)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22"
           "%20OR%20instrument_platform%3D%22PACBIO_SMRT%22)"
           "&fields=sample_accession&format=tsv&limit=0")
    lr = pd.read_csv(io.BytesIO(urlopen(url, timeout=600).read()), sep="\t", dtype=str)
    lr_bs = set(lr["sample_accession"].dropna())
    return our[our["sample_accession"].isin(lr_bs)]


def shv_copies_in_assembly(asm_fa, shv_query_fa, scratch):
    db = os.path.join(scratch, "asm")
    subprocess.run([MAKEBLASTDB, "-in", asm_fa, "-dbtype", "nucl", "-out", db],
                   capture_output=True, check=True)
    res = subprocess.run(
        [BLASTN, "-query", shv_query_fa, "-db", db, "-outfmt",
         "6 pident length qlen sstart send", "-perc_identity", "95"],
        capture_output=True, text=True, check=True)
    copies = 0
    for ln in res.stdout.splitlines():
        f = ln.split("\t")
        pident, length, qlen = float(f[0]), int(f[1]), int(f[2])
        if pident >= 95 and length / qlen >= 0.90:
            copies += 1
    return copies


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    shv = extract_shv_cds()
    shv_fa = OUT / "blaSHV_cds.fasta"
    shv_fa.write_text(f">blaSHV\n{shv}\n")

    ov = overlap_biosamples()
    print(f"overlap (10K cohort × long-read): {len(ov)} biosamples")
    calls = pd.read_csv(GENE_CALLS, sep="\t")[["sample_id", "crr_blaSHV", "blaSHV"]]
    ov = ov.merge(calls, on="sample_id", how="inner")
    print(f"  with a CNVRock blaSHV call: {len(ov)}")

    rows = []
    for _, r in ov.iterrows():
        bs = r["sample_accession"]
        with tempfile.TemporaryDirectory() as scratch:
            try:
                subprocess.run(
                    ["datasets", "download", "genome", "accession", bs,
                     "--assembly-source", "all", "--include", "genome",
                     "--filename", f"{scratch}/a.zip"],
                    capture_output=True, timeout=120, check=True)
                subprocess.run(["unzip", "-o", f"{scratch}/a.zip", "-d", scratch],
                               capture_output=True, check=True)
                fnas = glob.glob(f"{scratch}/**/*.fna", recursive=True)
                if not fnas:
                    continue
                copies = shv_copies_in_assembly(fnas[0], str(shv_fa), scratch)
            except Exception:
                continue
        rows.append({"sample_id": r["sample_id"], "biosample": bs,
                     "asm_shv_copies": copies,
                     "cnvrock_crr_blaSHV": r["crr_blaSHV"]})
        if len(rows) % 20 == 0:
            print(f"  processed {len(rows)} assemblies…", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "longread_shv_validation.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT/'longread_shv_validation.tsv'} ({len(df)} isolates)")
    if len(df) >= 10:
        from scipy import stats
        d = df.dropna(subset=["cnvrock_crr_blaSHV", "asm_shv_copies"])
        rho = stats.spearmanr(d["asm_shv_copies"], d["cnvrock_crr_blaSHV"]).statistic
        multi = d[d["asm_shv_copies"] >= 2]
        single = d[d["asm_shv_copies"] == 1]
        print(f"Spearman(asm copies, CNVRock CRR) = {rho:.3f}  (n={len(d)})")
        print(f"  CRR median — assembly single-copy: {single['cnvrock_crr_blaSHV'].median():.2f}"
              f"  multi-copy: {multi['cnvrock_crr_blaSHV'].median():.2f}")
        print(f"  assembly multi-copy isolates: {len(multi)}")


if __name__ == "__main__":
    main()
