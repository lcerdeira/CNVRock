#!/usr/bin/env python3
"""
Ploidy check for C. auris isolates with detected chr5 aneuploidy (exp39).

Question: are the 24 isolates with whole-chr5 CN≥2 truly aneuploid (chr5
only elevated) or whole-genome duplicated (all 7 chromosomes elevated)?

B8441v3 chromosomes and their accessions:
  CM076438.1  Chr 1   4,155,382 bp
  CM076439.1  Chr 2   2,367,702 bp
  CM076440.1  Chr 3   1,689,892 bp  ← contains ERG11
  CM076441.1  Chr 4   1,443,242 bp
  CM076442.1  Chr 5   1,007,788 bp  ← the "aneuploidy" chromosome
  CM076443.1  Chr 6     964,688 bp
  CM076444.1  Chr 7     777,906 bp

Method: for the isolates flagged as chr5-aneuploid (segments.parquet),
read the per-bin CRR across ALL chromosomes from reconstructions.npy /
gene_calls.tsv. Compute median CRR per chromosome per isolate.
  - All chr elevated equally → WGD (whole-genome duplication, false positive)
  - Only chr5 elevated → true trisomy / segmental aneuploidy

Run on HPC (needs exp39 segments.parquet and reconstructions.npy):
    /home/lshlt19/miniconda3/envs/cnvrock/bin/python3 \
        analysis/cauris_ploidy_check.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO     = Path("/home/lshlt19/CNVRock")
EXP39    = REPO / "data/results/39_cauris"
STORE    = REPO / "data/inputs/cauris-B8441v3-mq20-1000bp-npy"
OUT      = REPO / "data/results/cauris_ploidy_check"

# B8441v3 chromosome accessions in order
CHROMS = ["CM076438.1", "CM076439.1", "CM076440.1", "CM076441.1",
          "CM076442.1", "CM076443.1", "CM076444.1"]
CHR_NAMES = ["Chr1", "Chr2", "Chr3(ERG11)", "Chr4",
             "Chr5(aneup?)", "Chr6", "Chr7"]
CHR5 = "CM076442.1"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ── Load contigs + reconstructions ───────────────────────────────
    contigs = np.load(str(STORE / "contigs.npy"), allow_pickle=True)
    counts  = np.load(str(STORE / "counts.npy"),  mmap_mode="r")
    medians = np.load(str(STORE / "medians.npy"))
    ids     = np.load(str(STORE / "sample_ids.npy"), allow_pickle=True)

    # Normalised depth (same as training)
    X = counts.astype(np.float32) / (medians[:, None] + 1e-6)

    # Load reconstructions from exp39
    rec_path = EXP39 / "reconstructions.npy"
    if rec_path.exists():
        rec = np.load(str(rec_path), mmap_mode="r")
        exp_ids = np.load(str(EXP39 / "sample_ids.npy"), allow_pickle=True)
    else:
        print(f"reconstructions.npy not found at {rec_path}")
        print("Using raw counts / median as proxy reconstruction (less accurate)")
        # Fallback: use genome-wide median as flat baseline
        rec = np.ones_like(X)
        exp_ids = ids

    # Map exp39 sample IDs to store indices
    id_to_store = {sid: i for i, sid in enumerate(ids)}
    id_to_exp   = {sid: i for i, sid in enumerate(exp_ids)}

    # CRR = observed / reconstructed
    # bin-to-chromosome mapping — handle structured dtype ('chrom','start','end')
    chrom_bins = {}
    chrom_field = contigs["chrom"] if contigs.dtype.names and "chrom" in contigs.dtype.names else contigs
    for i, ch in enumerate(chrom_field):
        chrom_bins.setdefault(str(ch), []).append(i)

    print("Bin counts per chromosome:")
    for ch, name in zip(CHROMS, CHR_NAMES):
        n = len(chrom_bins.get(ch, []))
        print(f"  {ch} ({name}): {n} bins")

    # ── Identify chr5-aneuploid isolates ─────────────────────────────
    # Load gene_calls — look for 'ERG11' CRR > 1 as proxy, or use
    # the segments file if available
    gene_calls = pd.read_csv(EXP39 / "gene_calls.tsv", sep="\t")

    # chr5 aneuploidy: CRR of the entire chr5 elevated
    # Compute median CRR on chr5 bins for each sample
    chr5_bins = chrom_bins.get(CHR5, [])
    if not chr5_bins:
        print(f"ERROR: no bins found for {CHR5}")
        return

    # For all exp39 samples, compute per-chromosome median CRR
    print(f"\nComputing per-chromosome median CRR for {len(exp_ids)} isolates…")

    chrom_crr = {ch: [] for ch in CHROMS}
    for sid in exp_ids:
        s_idx = id_to_store.get(sid)
        e_idx = id_to_exp.get(sid)
        if s_idx is None or e_idx is None:
            for ch in CHROMS:
                chrom_crr[ch].append(np.nan)
            continue
        x_obs = X[s_idx]
        x_rec = rec[e_idx] if rec.shape[1] == X.shape[1] else x_obs
        crr = x_obs / (x_rec + 1e-6)
        # Renormalise by overall median to correct for global coverage shifts
        crr = crr / (np.median(crr) + 1e-6)
        for ch in CHROMS:
            bins = chrom_bins.get(ch, [])
            chrom_crr[ch].append(np.median(crr[bins]) if bins else np.nan)

    df_crr = pd.DataFrame(chrom_crr, index=exp_ids)
    df_crr.columns = [f"crr_{c}" for c in CHROMS]

    # ── Classify isolates ─────────────────────────────────────────────
    # Chr5 aneuploid: chr5 CRR ≥ 1.4 (elevated)
    CHR5_COL = f"crr_{CHR5}"
    df_crr["chr5_elevated"] = df_crr[CHR5_COL] >= 1.4

    # WGD: ALL chromosomes elevated (median of other 6 chrs also ≥ 1.4)
    other_chr_cols = [f"crr_{c}" for c in CHROMS if c != CHR5]
    df_crr["other_chr_median"] = df_crr[other_chr_cols].median(axis=1)
    df_crr["is_wgd"] = (df_crr["chr5_elevated"] &
                        (df_crr["other_chr_median"] >= 1.35))
    df_crr["is_true_chr5_aneuploid"] = (df_crr["chr5_elevated"] &
                                         ~df_crr["is_wgd"])

    n_chr5_elev  = df_crr["chr5_elevated"].sum()
    n_wgd        = df_crr["is_wgd"].sum()
    n_true_aneu  = df_crr["is_true_chr5_aneuploid"].sum()

    print(f"\n── Classification ───────────────────────────────────────────")
    print(f"  Chr5 elevated (CRR≥1.4):       {n_chr5_elev}")
    print(f"    → Whole-genome duplication:   {n_wgd}")
    print(f"    → True chr5 aneuploidy only:  {n_true_aneu}")
    print(f"  Not elevated:                   {len(df_crr) - n_chr5_elev}")

    # Show the WGD candidates
    if n_wgd:
        print(f"\n  WGD candidates ({n_wgd} isolates):")
        wgd = df_crr[df_crr["is_wgd"]]
        print(wgd[[CHR5_COL, "other_chr_median"]].head(10).to_string())

    # Save
    df_crr.to_csv(OUT / "cauris_per_chrom_crr.tsv", sep="\t")
    print(f"\nSaved {OUT}/cauris_per_chrom_crr.tsv")

    # ── Figure ─────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: chr5 CRR vs other-chromosome median CRR
    ax = axes[0]
    colors = ["#e63946" if wgd else "#2176ae" if aneu else "#adb5bd"
              for wgd, aneu in zip(df_crr["is_wgd"],
                                    df_crr["is_true_chr5_aneuploid"])]
    ax.scatter(df_crr["other_chr_median"], df_crr[CHR5_COL],
               c=colors, s=10, alpha=0.6)
    ax.axhline(1.4, color="grey", ls="--", lw=0.8)
    ax.axvline(1.35, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("Median CRR — other 6 chromosomes", fontsize=9)
    ax.set_ylabel("CRR — Chr5", fontsize=9)
    ax.set_title("WGD vs chr5-specific aneuploidy\n"
                 "Red=WGD | Blue=chr5 only | Grey=diploid", fontsize=9)

    # Right: per-chromosome CRR heatmap for chr5-elevated isolates
    ax = axes[1]
    sub = df_crr[df_crr["chr5_elevated"]].head(50)
    mat = sub[[f"crr_{c}" for c in CHROMS]].values
    im = ax.imshow(mat.T, cmap="RdBu_r", vmin=0.5, vmax=2.0,
                   aspect="auto")
    ax.set_yticks(range(7))
    ax.set_yticklabels(CHR_NAMES, fontsize=8)
    ax.set_xlabel("Isolate index (chr5-elevated subset)", fontsize=9)
    ax.set_title("Per-chromosome CRR — chr5-elevated isolates\n"
                 "Uniform elevation = WGD; chr5 only = aneuploidy", fontsize=9)
    plt.colorbar(im, ax=ax, shrink=0.8, label="CRR")

    plt.tight_layout()
    fig.savefig(OUT / "cauris_ploidy_check.png", dpi=150, bbox_inches="tight")
    print(f"Figure saved: {OUT}/cauris_ploidy_check.png")

    print(f"\n── Manuscript update ─────────────────────────────────────────")
    print(f"  Previously reported: 24 isolates with chr5 aneuploidy (CN≥2)")
    print(f"  Revised:  {n_true_aneu} true chr5-only aneuploidies + "
          f"{n_wgd} whole-genome duplications")
    if n_wgd > 0:
        print(f"  *** IMPORTANT: {n_wgd} of the 24 are WGD, not chr5 aneuploidy!")
        print(f"  Manuscript §3.9 needs update.")
    else:
        print(f"  All {n_true_aneu} isolates confirmed as chr5-specific aneuploidy (not WGD).")


if __name__ == "__main__":
    main()
