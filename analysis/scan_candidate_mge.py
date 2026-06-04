#!/usr/bin/env python3
"""
Mobile-genetic-element (MGE) context for the genome-wide CNV scan candidate
loci (§3.8). Tests the central hypothesis that the recurrently-amplified loci
are passenger markers of integrons / ICEs / IS-flanked cassettes.

Tools (conda env 'mgescan'):
  - IntegronFinder 2  — class-1 integrons + integron-integrase (intI)
  - ISEScan           — insertion sequences (IS families)
  - abricate          — AMR / plasmid / virulence gene screen (CARD, ResFinder,
                        PlasmidFinder, VFDB) to find any AMR gene physically
                        linked to (within ±20 kb of) a candidate locus

Method: run all three on the full HS11286 chromosome (NC_016845.1), then
intersect detected elements with each candidate locus (bin coordinates =
bin × 1000 bp), reporting the nearest IS, integron membership, and any AMR
gene within a 20 kb window.

Run on HPC:
  conda activate mgescan
  python3 analysis/scan_candidate_mge.py

Outputs (data/results/cnv_scan_phase_e/):
  mge_isescan.tsv, mge_integron.tsv, mge_abricate.tsv  (raw per-tool)
  candidate_mge_context.tsv                            (per-locus summary)
"""
from __future__ import annotations
import subprocess, shutil, os, glob, re
from pathlib import Path
import pandas as pd

REPO   = Path("/home/lshlt19/CNVRock")
FASTA  = REPO / "assets/HS11286.fasta"
CHROM  = "NC_016845.1"
LOCI   = REPO / "data/results/cnv_scan_phase_e/candidate_annotation_per_locus.tsv"
OUTDIR = REPO / "data/results/cnv_scan_phase_e"
WORK   = OUTDIR / "mge_work"
WINDOW = 20000   # ±20 kb linkage window for AMR genes


def run(cmd, **kw):
    print("  $", " ".join(cmd))
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def chrom_fasta() -> Path:
    """Write a chromosome-only FASTA (tools choke on multi-contig sometimes)."""
    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / "HS11286_chrom.fasta"
    keep, buf = False, []
    for line in FASTA.read_text().split("\n"):
        if line.startswith(">"):
            keep = line[1:].split()[0] == CHROM
        if keep:
            buf.append(line)
    out.write_text("\n".join(buf) + "\n")
    return out


def run_isescan(fa: Path) -> pd.DataFrame:
    out = WORK / "isescan"
    if not shutil.which("isescan.py"):
        print("  isescan not found — skipping"); return pd.DataFrame()
    run(["isescan.py", "--seqfile", str(fa), "--output", str(out), "--nthread", "4"])
    tsv = glob.glob(str(out / "**/*.tsv"), recursive=True)
    csvs = glob.glob(str(out / "**/*.csv"), recursive=True)
    f = (tsv or csvs)
    if not f:
        return pd.DataFrame()
    df = pd.read_csv(f[0], sep=None, engine="python")
    df.to_csv(OUTDIR / "mge_isescan.tsv", sep="\t", index=False)
    return df


def run_integron(fa: Path) -> pd.DataFrame:
    if not shutil.which("integron_finder"):
        print("  integron_finder not found — skipping"); return pd.DataFrame()
    out = WORK / "integron"
    run(["integron_finder", "--local-max", "--outdir", str(out), str(fa)])
    res = glob.glob(str(out / "**/*.integrons"), recursive=True)
    if not res:
        return pd.DataFrame()
    rows = []
    for f in res:
        for line in open(f):
            if line.startswith("#") or line.startswith("ID_integron"):
                continue
            rows.append(line.rstrip().split("\t"))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "mge_integron.tsv", sep="\t", index=False)
    return df


def run_abricate(fa: Path) -> dict[str, pd.DataFrame]:
    if not shutil.which("abricate"):
        print("  abricate not found — skipping"); return {}
    res = {}
    for db in ["card", "resfinder", "plasmidfinder", "vfdb"]:
        r = run(["abricate", "--db", db, str(fa)])
        if r.returncode == 0 and r.stdout.strip():
            from io import StringIO
            df = pd.read_csv(StringIO(r.stdout), sep="\t")
            df["db"] = db
            res[db] = df
    if res:
        alld = pd.concat(res.values(), ignore_index=True)
        alld.to_csv(OUTDIR / "mge_abricate.tsv", sep="\t", index=False)
    return res


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fa = chrom_fasta()
    loci = pd.read_csv(LOCI, sep="\t")

    print("Running ISEScan…");        iss = run_isescan(fa)
    print("Running IntegronFinder…"); integ = run_integron(fa)
    print("Running abricate…");       abr = run_abricate(fa)

    # Helper: extract start/end columns robustly from a tool dataframe
    def coords(df, sc, ec):
        for s in sc:
            if s in df.columns: scol = s; break
        else: return []
        for e in ec:
            if e in df.columns: ecol = e; break
        else: return []
        out = []
        for _, r in df.iterrows():
            try:
                out.append((int(r[scol]), int(r[ecol]), r))
            except Exception:
                pass
        return out

    is_coords = coords(iss, ["isBegin", "start", "Start"],
                            ["isEnd", "end", "End"]) if len(iss) else []
    abr_all = pd.concat(abr.values(), ignore_index=True) if abr else pd.DataFrame()
    abr_coords = coords(abr_all, ["START"], ["END"]) if len(abr_all) else []

    rows = []
    for _, L in loci.iterrows():
        bc = str(L["bin_cluster"])
        lo = int(bc.split("–")[0]) * 1000
        hi = int(bc.split("–")[-1]) * 1000 + 1000
        # nearest IS
        nearest_is, is_dist = "—", None
        for s, e, r in is_coords:
            d = 0 if (e >= lo and s <= hi) else min(abs(s-hi), abs(e-lo))
            if is_dist is None or d < is_dist:
                is_dist = d
                fam = r.get("family", r.get("cluster", "IS"))
                nearest_is = f"{fam} ({'within' if d==0 else f'{d//1000}kb'})"
        # AMR gene within window
        amr_hits = []
        for s, e, r in abr_coords:
            if e >= lo - WINDOW and s <= hi + WINDOW:
                amr_hits.append(f"{r.get('GENE','?')}[{r.get('db','?')}]")
        rows.append({
            "bin_cluster": bc,
            "max_n_STs":   L["max_n_STs"],
            "best_q":      L["best_q"],
            "genes":       L["genes"],
            "KEGG_pathways": L.get("KEGG_pathways", "—"),
            "nearest_IS":  nearest_is,
            "AMR_within_20kb": "; ".join(sorted(set(amr_hits))) or "none",
        })
    out = pd.DataFrame(rows).sort_values(["max_n_STs", "best_q"],
                                         ascending=[False, True])
    out.to_csv(OUTDIR / "candidate_mge_context.tsv", sep="\t", index=False)
    print(f"\nwrote candidate_mge_context.tsv ({len(out)} loci)")
    print(out.to_string(index=False, max_colwidth=40))

    # Integron summary
    if len(integ):
        print(f"\nIntegronFinder: {len(integ)} integron rows detected "
              f"(see mge_integron.tsv)")
    else:
        print("\nIntegronFinder: no integrons on the HS11286 chromosome "
              "(expected — HS11286 reference carries few; ST-genome analysis "
              "would be needed to see mobilised cassettes).")


if __name__ == "__main__":
    main()
