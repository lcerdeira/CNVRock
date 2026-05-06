"""
Gene-level CNV caller — version 06: KpSC chromosomal genes.

Identical calling logic to 05_gene_cnv_caller (fractional gene coverage for
long genes, CRR fallback for short genes) adapted for K. pneumoniae Species
Complex (KpSC) targets on the HS11286 chromosome (NC_016845.1, RefSeq).

Genes of interest (all chromosomally encoded in KpSC):
  blaSHV   — ancestral chromosomal beta-lactamase; amplification increases
              resistance to penicillins and early cephalosporins.
  ramR     — repressor of acrAB efflux operon; deletion → MDR via AcrAB-TolC.

NOT included — requires sequence-level variant calling, not read-depth CNV:
  ompK35   — loss-of-function via frameshift/truncation (not copy-number change);
              CRR is identical between functional and disrupted samples.
  ompK36   — same reason as ompK35.
  Use Kleborate or a short-read variant caller (e.g. snippy) for porin status.

Config keys used (same semantics as 05_gene_cnv_caller):
    cnv_crr_gate_threshold        — veto gate for long-gene fractional path
    cnv_crr_amp_threshold         — CRR threshold for short-gene (fallback) calls
    cnv_min_cn1_proportion        — min fraction of chrom bins in CN=1 for sanity
    cnv_min_confidence            — minimum HMM posterior to count a segment
    cnv_flank_padding             — bp on each side of gene for CRR baseline
    cnv_crr_min_bins_fallback     — genes with < this many bins use CRR-only path
    cnv_min_gene_coverage_fraction — min fraction of gene bins covered by CN>=2
"""

import os
import gc
import time

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# KpSC genes of interest — HS11286 chromosome (NC_016845.1)
# ---------------------------------------------------------------------------

CHROM = "NC_016845.1"  # HS11286 chromosome (RefSeq GCF download) — confirm with: samtools view -H *.bam | grep @SQ

GENES_OF_INTEREST = [
    # blaSHV-11: KPHS_25220, 861 bp, minus strand (GFF 1-based coords)
    {"call_id": "blaSHV",  "contig": CHROM, "start": 2549403, "end": 2550263},
    # ramR removed: disruption in KpSC is almost exclusively via IS element insertion
    # or point mutation — not physical deletion.  crr_ramR never drops below 0.70
    # across 545 samples, confirming no read-depth-detectable deletions exist.
    # Evaluate ramR status using Kleborate GT (assembly + BLAST-based).
    # ompK35 and ompK36 removed for the same reason (frameshift/truncation).
]


# ---------------------------------------------------------------------------
# Calling logic (identical to 05_gene_cnv_caller)
# ---------------------------------------------------------------------------

def run_cnv_calls(store_path, out_dir, cfg):
    """Compute gene-level CNV calls for all KpSC samples and write gene_calls.tsv."""
    min_cn1_proportion     = cfg["cnv_min_cn1_proportion"]
    min_confidence         = cfg["cnv_min_confidence"]
    flank_padding          = cfg["cnv_flank_padding"]
    crr_amp_threshold      = cfg["cnv_crr_amp_threshold"]
    crr_min_bins_fallback  = cfg["cnv_crr_min_bins_fallback"]
    crr_gate_threshold     = cfg["cnv_crr_gate_threshold"]
    min_gene_coverage_frac = cfg["cnv_min_gene_coverage_fraction"]

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)
    segments   = pd.read_parquet(os.path.join(out_dir, "segments.parquet"))

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    del contigs
    gc.collect()

    n                = len(sample_ids)
    low_cov          = cfg.get("hmm_low_cov_threshold", 10)
    sample_mean      = counts.astype(float).mean(axis=1, keepdims=True)
    safe_recon       = np.where(recons.astype(float) >= low_cov,
                                recons.astype(float), sample_mean)
    copy_ratios = counts.astype(float) / (safe_recon + 1e-6)
    del recons, safe_recon, sample_mean
    gc.collect()

    print(f"Computing gene CNV calls for {n} KpSC samples…", flush=True)
    t0 = time.time()

    gene_contigs  = {g["contig"] for g in GENES_OF_INTEREST}
    segs_filtered = segments[segments["chrom"].isin(gene_contigs)]
    empty_segs    = pd.DataFrame(columns=["chrom", "x0", "x1", "cn", "confidence"])
    segs_by_sample = {sid: grp for sid, grp in segs_filtered.groupby("sample_id")}
    del segments, segs_filtered
    gc.collect()

    gene_rows = []
    for gene in GENES_OF_INTEREST:
        contig  = gene["contig"]
        g_start = gene["start"]
        g_end   = gene["end"]
        call_id = gene["call_id"]

        chrom_mask  = chroms == contig
        s           = starts[chrom_mask]
        gene_mask_l = (s >= g_start) & (s <= g_end)
        flank_mask_l = (s < g_start - flank_padding) | (s > g_end + flank_padding)
        n_chrom      = int(chrom_mask.sum())
        n_gene_bins  = int(gene_mask_l.sum())
        fallback_eligible = (n_gene_bins < crr_min_bins_fallback)

        gene_bins = s[gene_mask_l]

        cr_chrom   = copy_ratios[:, chrom_mask]
        mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)
        mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            crr_all = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)

        for idx, sid in enumerate(sample_ids):
            cn          = -1
            sanity_ok   = False
            sample_segs = segs_by_sample.get(sid, empty_segs)
            chrom_segs  = sample_segs[sample_segs["chrom"] == contig]
            if len(chrom_segs) > 0:
                hc1 = chrom_segs[
                    (chrom_segs["cn"] == 1) & (chrom_segs["confidence"] >= min_confidence)
                ]
                if len(hc1) > 0:
                    x0 = hc1["x0"].values
                    x1 = hc1["x1"].values
                    covered = ((s[:, None] >= x0) & (s[:, None] < x1)).any(axis=1)
                    if covered.sum() / n_chrom >= min_cn1_proportion:
                        sanity_ok = True
                        if fallback_eligible:
                            spanning = chrom_segs[
                                (chrom_segs["x0"] <= g_start) & (chrom_segs["x1"] >= g_end)
                            ]
                            if len(spanning) > 0:
                                cn = int(spanning.iloc[0]["cn"])

            crr_val = float(crr_all[idx]) if np.isfinite(crr_all[idx]) else None

            if sanity_ok:
                if fallback_eligible:
                    if cn in (-1, 1) and crr_val is not None and crr_val >= crr_amp_threshold:
                        cn = 2
                else:
                    gene_segs_amp = chrom_segs[
                        (chrom_segs["cn"] >= 2) &
                        (chrom_segs["confidence"] >= min_confidence) &
                        (chrom_segs["x0"] < g_end) &
                        (chrom_segs["x1"] > g_start)
                    ]
                    if n_gene_bins > 0 and len(gene_segs_amp) > 0:
                        x0_a = gene_segs_amp["x0"].values
                        x1_a = gene_segs_amp["x1"].values
                        covered_g = (
                            (gene_bins[:, None] >= x0_a) & (gene_bins[:, None] < x1_a)
                        ).any(axis=1)
                        if covered_g.sum() / n_gene_bins >= min_gene_coverage_frac:
                            cn = 2
                            if crr_val is None or crr_val < crr_gate_threshold:
                                cn = 1

            gene_rows.append({"sample_id": sid, "call_id": call_id, "cn": cn, "crr": crr_val})

        print(f"  Done {call_id} | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, segs_by_sample
    gc.collect()

    genes = [g["call_id"] for g in GENES_OF_INTEREST]
    gene_calls_long = pd.DataFrame(gene_rows)[["sample_id", "call_id", "cn", "crr"]]
    cn_wide  = gene_calls_long.pivot(index="sample_id", columns="call_id", values="cn")
    crr_wide = gene_calls_long.pivot(index="sample_id", columns="call_id", values="crr")
    cn_wide.columns.name  = None
    crr_wide.columns      = [f"crr_{g}" for g in crr_wide.columns]
    crr_wide.columns.name = None
    gene_calls_wide = pd.concat([cn_wide, crr_wide], axis=1).reset_index()

    del gene_calls_long, cn_wide, crr_wide, gene_rows
    gc.collect()

    gene_path = os.path.join(out_dir, "gene_calls.tsv")
    gene_calls_wide.to_csv(gene_path, sep="\t", index=False)
    print(f"Saved gene calls ({len(gene_calls_wide):,} samples) → {gene_path}", flush=True)
