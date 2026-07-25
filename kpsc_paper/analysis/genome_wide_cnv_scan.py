#!/usr/bin/env python3
"""
Genome-wide CNV–phenotype scan (Phase E, exploratory).

For every 1 kb bin in NC_016845.1 we compute:

  1. Per-bin copy ratio  r = x / x̂   from the chrom MQ20 NPY store and the
     VAE reconstructions (data/results/33_kpsc_expansion_10k/reconstructions.npy).
  2. Per-bin prevalence-of-amplification  (fraction of samples with r > 1.75).
  3. Per-(bin × antibiotic) Mann-Whitney U test on  r | R  vs  r | S
     using CABBAGE phenotypes joined via the run → BioSample bridge.
     FDR-corrected across 5,334 bins per drug with Benjamini–Hochberg.
  4. Bin → gene lift via HS11286.gff3 (chromosome NC_016845.1).
  5. Cross-ST recurrence: number of distinct STs in which the bin shows
     amplification (controls for clonal expansion vs convergent selection).

Output: tables ranking candidate loci by (-log10 p_FDR) × recurrence × prevalence.

Run on HPC:
    sbatch hpc/run_genome_wide_scan.sh
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_observed_depth(store_dir: Path, sample_ids: np.ndarray) -> np.ndarray:
    """Load the chrom NPY store and reorder rows to match sample_ids.

    The store layout is a single tensor `counts.npy` (n_samples × n_bins)
    plus `sample_ids.npy` giving the row order, plus `contigs.npy` and
    `medians.npy`. We slice to NC_016845.1 (the chromosome) bins only.
    """
    counts_path = store_dir / "counts.npy"
    sids_path   = store_dir / "sample_ids.npy"
    contigs_path = store_dir / "contigs.npy"
    if not counts_path.exists():
        raise FileNotFoundError(f"missing {counts_path}")
    counts = np.load(counts_path).astype(np.float32)
    store_sids = np.array([str(s) for s in np.load(sids_path, allow_pickle=True)])
    if contigs_path.exists():
        contigs = np.load(contigs_path, allow_pickle=True)
        contigs = np.array([str(c) for c in contigs])
        chrom_mask = (contigs == "NC_016845.1")
        if chrom_mask.any():
            counts = counts[:, chrom_mask]
            print(f"      sliced to NC_016845.1: {chrom_mask.sum()} bins")
    # reorder rows to match requested sample_ids
    sid_to_row = {s: i for i, s in enumerate(store_sids)}
    rows = [sid_to_row.get(str(s), -1) for s in sample_ids]
    out = np.full((len(sample_ids), counts.shape[1]), np.nan, dtype=np.float32)
    found = np.array([r for r in rows if r >= 0])
    keep  = np.array([i for i, r in enumerate(rows) if r >= 0])
    out[keep] = counts[found]
    print(f"      matched {len(found):,} / {len(sample_ids):,} samples from store")
    return out


def per_sample_median_normalise(x: np.ndarray) -> np.ndarray:
    med = np.nanmedian(x, axis=1, keepdims=True)
    med = np.where(med == 0, np.nan, med)
    return x / med


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = q
    return np.clip(out, 0, 1)


def parse_gff3(gff_path: Path) -> pd.DataFrame:
    rows = []
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip().split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append({
                "chrom":  f[0],
                "start":  int(f[3]),
                "end":    int(f[4]),
                "strand": f[6],
                "gene":   attrs.get("Name") or attrs.get("gene") or attrs.get("ID", ""),
                "biotype": attrs.get("gene_biotype", ""),
            })
    return pd.DataFrame(rows)


def lift_bin_to_gene(bin_idx: int, gff_chrom: pd.DataFrame, bin_size: int = 1000) -> str:
    start, end = bin_idx * bin_size + 1, (bin_idx + 1) * bin_size
    hits = gff_chrom[(gff_chrom["start"] <= end) & (gff_chrom["end"] >= start)]
    if hits.empty:
        return "intergenic"
    return ";".join(hits["gene"].astype(str).unique())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store-dir",  type=Path,
                    default=Path("data/inputs/KpSC-expansion-10k-mq20-1000bp-npy"))
    ap.add_argument("--results-dir", type=Path,
                    default=Path("data/results/33_kpsc_expansion_10k"))
    ap.add_argument("--meta", type=Path,
                    default=Path("assets/kpsc_expansion_metadata_runlevel.tsv"))
    ap.add_argument("--cabbage", type=Path,
                    default=Path("assets/cabbage_kpsc_phenotypes.tsv"))
    ap.add_argument("--gff", type=Path,
                    default=Path("assets/HS11286.gff3"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("data/results/cnv_scan_phase_e"))
    ap.add_argument("--crr-amp", type=float, default=1.75,
                    help="bin-level amplification threshold")
    ap.add_argument("--min-recurrence-st", type=int, default=2,
                    help="minimum distinct STs in which a bin must be amplified")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load observed + reconstructed depth ────────────────────────────
    sample_ids = np.load(args.results_dir / "sample_ids.npy", allow_pickle=True)
    sample_ids = np.array([str(s) for s in sample_ids])
    print(f"[1/6] loading observed depth for {len(sample_ids):,} samples…")
    x = load_observed_depth(args.store_dir, sample_ids)
    x = per_sample_median_normalise(x)
    print(f"      observed shape={x.shape}, nan-frac={np.isnan(x).mean():.3f}")

    xhat = np.load(args.results_dir / "reconstructions.npy").astype(np.float32)
    print(f"[2/6] reconstructions shape={xhat.shape}")

    # Per-bin copy ratio
    with np.errstate(invalid="ignore", divide="ignore"):
        r = x / np.where(xhat == 0, np.nan, xhat)
    print(f"      copy-ratio computed; valid-frac={(~np.isnan(r)).mean():.3f}")

    # ── 2. Per-bin prevalence of amplification ────────────────────────────
    print(f"[3/6] per-bin prevalence (r > {args.crr_amp})…")
    n_amp = (r > args.crr_amp).sum(axis=0)
    n_obs = (~np.isnan(r)).sum(axis=0)
    prev  = n_amp / np.maximum(n_obs, 1)

    # ── 3. Per-bin × per-antibiotic Mann-Whitney U ───────────────────────
    print("[4/6] phenotype linkage + per-bin Mann-Whitney scan…")
    meta = pd.read_csv(args.meta, sep="\t")
    if "sample_accession" not in meta.columns:
        print("ERROR: meta lacks sample_accession", file=sys.stderr); sys.exit(1)
    meta = meta[["sample_id", "sample_accession", "ST"]].rename(
        columns={"sample_accession": "biosample"})
    cabbage = pd.read_csv(args.cabbage, sep="\t",
                          usecols=["BioSample_ID", "antibiotic_name",
                                   "updated_phenotype", "resistance_phenotype"])
    cabbage["pheno"] = (cabbage["updated_phenotype"].fillna(cabbage["resistance_phenotype"])
                        .astype(str).str.lower().str[:1].str.upper())
    cabbage = cabbage[cabbage["pheno"].isin(["R", "S"])]
    pheno = (cabbage.rename(columns={"BioSample_ID": "biosample"})
             .merge(meta, on="biosample", how="inner"))
    # row order = sample_ids
    sid_idx = pd.Series(np.arange(len(sample_ids)), index=sample_ids)
    pheno["row"] = pheno["sample_id"].map(sid_idx)
    pheno = pheno.dropna(subset=["row"])
    pheno["row"] = pheno["row"].astype(int)
    print(f"      phenotype rows after join: {len(pheno):,}  "
          f"({pheno['antibiotic_name'].nunique()} antibiotics)")

    drug_results = []
    for drug, g in pheno.groupby("antibiotic_name"):
        r_idx = g.loc[g["pheno"] == "R", "row"].values
        s_idx = g.loc[g["pheno"] == "S", "row"].values
        if len(r_idx) < 10 or len(s_idx) < 10:
            continue
        rR = r[r_idx, :]
        rS = r[s_idx, :]
        med_R = np.nanmedian(rR, axis=0)
        med_S = np.nanmedian(rS, axis=0)
        # vectorised Mann-Whitney U per bin (one-sided: R > S)
        pvals = np.full(rR.shape[1], np.nan)
        for b in range(rR.shape[1]):
            a = rR[:, b]; bb = rS[:, b]
            a = a[~np.isnan(a)]; bb = bb[~np.isnan(bb)]
            if len(a) < 5 or len(bb) < 5: continue
            try:
                _, p = stats.mannwhitneyu(a, bb, alternative="greater")
                pvals[b] = p
            except ValueError:
                pass
        q = np.where(np.isnan(pvals), np.nan, bh_fdr(np.nan_to_num(pvals, nan=1.0)))
        drug_results.append(pd.DataFrame({
            "antibiotic": drug,
            "bin": np.arange(rR.shape[1]),
            "n_R": len(r_idx), "n_S": len(s_idx),
            "med_r_R": med_R, "med_r_S": med_S,
            "p": pvals, "q": q,
            "prev_amp": prev,
            "n_amp": n_amp,
        }))
    scan = pd.concat(drug_results, ignore_index=True)

    # ── 4. Cross-ST recurrence ────────────────────────────────────────────
    print("[5/6] cross-ST amplification recurrence…")
    st_map = meta.set_index("sample_id")["ST"].reindex(sample_ids).values
    df_amp = (r > args.crr_amp)
    st_series = pd.Series(st_map)
    n_sts = np.zeros(r.shape[1], dtype=np.int32)
    for b in range(r.shape[1]):
        n_sts[b] = st_series[df_amp[:, b]].nunique()
    scan["n_sts_amplified"] = scan["bin"].map(dict(zip(range(r.shape[1]), n_sts)))

    # ── 5. Gene lift ──────────────────────────────────────────────────────
    print("[6/6] gene-coordinate lift…")
    gene_map = {}
    if args.gff.exists():
        gff = parse_gff3(args.gff)
        gff_chrom = gff[gff["chrom"].str.contains("NC_016845")]
        for b in scan["bin"].unique():
            gene_map[int(b)] = lift_bin_to_gene(int(b), gff_chrom)
    else:
        print(f"  (skipping gene lift — {args.gff} not present)")
    scan["gene"] = scan["bin"].map(gene_map).fillna("?")

    # ── 6. Rank ───────────────────────────────────────────────────────────
    scan["neg_log10_q"] = -np.log10(np.clip(scan["q"].fillna(1), 1e-300, 1))
    scan["score"] = (scan["neg_log10_q"]
                     * np.maximum(scan["n_sts_amplified"], 1)
                     * np.maximum(scan["prev_amp"], 0.001))
    scan = scan.sort_values(["antibiotic", "score"], ascending=[True, False])

    # Save full + filtered tables
    full_path = args.out_dir / "scan_full.tsv"
    scan.to_csv(full_path, sep="\t", index=False)
    print(f"  full scan → {full_path}  ({len(scan):,} rows)")

    sig = scan[(scan["q"] < 0.05)
               & (scan["n_sts_amplified"] >= args.min_recurrence_st)]
    sig_path = args.out_dir / "scan_significant.tsv"
    sig.to_csv(sig_path, sep="\t", index=False)
    print(f"  significant (q<0.05, ≥{args.min_recurrence_st} STs) → {sig_path}  "
          f"({len(sig):,} rows across {sig['antibiotic'].nunique()} drugs)")

    # Top-50 candidates per drug
    top = scan.groupby("antibiotic", group_keys=False).head(50)
    top_path = args.out_dir / "scan_top50_per_drug.tsv"
    top.to_csv(top_path, sep="\t", index=False)
    print(f"  top-50 per drug → {top_path}")


if __name__ == "__main__":
    main()
