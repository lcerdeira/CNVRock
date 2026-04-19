"""
Gene-level CNV caller — version 04: CRR gate on long genes.

Extends 03 with a CRR gate for genes that span >= cnv_crr_min_bins_fallback bins.
For these genes the HMM call is normally accepted unconditionally; in version 03 the
cnv_crr_amp_threshold only operated via the fallback path (short genes). Version 04
adds a veto: if the HMM calls CN>=2 on a long gene but the gene/flank CRR is below
cnv_crr_gate_threshold, the call is downgraded to CN=1.

Motivation: 05_conv_vae generates a broader anomaly signal than 04_conv_vae, causing
the HMM to produce many CN>=2 segments at noise level (CRT FP p50=1.12, PM2_PM3 FP
p50=1.16). The gate suppresses these HMM-driven FPs without touching short-gene logic.

For short genes (< cnv_crr_min_bins_fallback bins, e.g. GCH1) the existing fallback
logic is unchanged: CRR >= cnv_crr_amp_threshold overrides CN=1 → CN=2.

New config key:
    cnv_crr_gate_threshold  — CRR below which a long-gene HMM CN>=2 call is rejected.

Existing keys unchanged:
    cnv_min_cn1_proportion
    cnv_min_confidence
    cnv_flank_padding
    cnv_crr_amp_threshold        — still used for the short-gene fallback
    cnv_crr_min_bins_fallback    — threshold separating short/long gene treatment
"""

import os

import numpy as np
import pandas as pd


GENES_OF_INTEREST = [
    {"call_id": "MDR1",    "contig": "Pf3D7_05_v3", "start": 955955,  "end": 963095},
    {"call_id": "CRT",     "contig": "Pf3D7_07_v3", "start": 402385,  "end": 406341},
    {"call_id": "GCH1",    "contig": "Pf3D7_12_v3", "start": 974226,  "end": 976097},
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 292244,  "end": 299101},
]


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute gene-level CNV calls for all samples and write gene_calls.tsv."""
    import gc
    import time

    min_cn1_proportion    = cfg["cnv_min_cn1_proportion"]
    min_confidence        = cfg["cnv_min_confidence"]
    flank_padding         = cfg["cnv_flank_padding"]
    crr_amp_threshold     = cfg["cnv_crr_amp_threshold"]
    crr_min_bins_fallback = cfg["cnv_crr_min_bins_fallback"]
    crr_gate_threshold    = cfg["cnv_crr_gate_threshold"]

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)
    segments   = pd.read_parquet(os.path.join(out_dir, "segments.parquet"))

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    del contigs
    gc.collect()

    n           = len(sample_ids)
    copy_ratios = counts.astype(float) / (recons.astype(float) + 1e-6)
    del recons
    gc.collect()

    print(f"Computing gene CNV calls for {n} samples…", flush=True)
    t0 = time.time()

    gene_contigs   = {g["contig"] for g in GENES_OF_INTEREST}
    segs_filtered  = segments[segments["chrom"].isin(gene_contigs)]
    empty_segs     = pd.DataFrame(columns=["chrom", "x0", "x1", "cn", "confidence"])
    segs_by_sample = {sid: grp for sid, grp in segs_filtered.groupby("sample_id")}
    del segments, segs_filtered
    gc.collect()

    gene_rows = []
    for gene in GENES_OF_INTEREST:
        contig  = gene["contig"]
        g_start = gene["start"]
        g_end   = gene["end"]
        call_id = gene["call_id"]

        chrom_mask    = chroms == contig
        s             = starts[chrom_mask]
        gene_mask_l   = (s >= g_start) & (s <= g_end)
        flank_mask_l  = (s < g_start - flank_padding) | (s > g_end + flank_padding)
        n_chrom       = int(chrom_mask.sum())
        n_gene_bins   = int(gene_mask_l.sum())
        fallback_eligible = (n_gene_bins < crr_min_bins_fallback)

        cr_chrom   = copy_ratios[:, chrom_mask]
        mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)
        mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)
        crr_all    = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)

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
                        spanning = chrom_segs[
                            (chrom_segs["x0"] <= g_start) & (chrom_segs["x1"] >= g_end)
                        ]
                        if len(spanning) > 0:
                            cn = int(spanning.iloc[0]["cn"])

            crr_val = float(crr_all[idx]) if np.isfinite(crr_all[idx]) else None

            if sanity_ok:
                if fallback_eligible:
                    # Short gene: CRR can override a CN=1 or failed HMM call.
                    if cn in (-1, 1) and crr_val is not None and crr_val >= crr_amp_threshold:
                        cn = 2
                else:
                    # Long gene: veto HMM amplification calls with low CRR.
                    if cn >= 2 and (crr_val is None or crr_val < crr_gate_threshold):
                        cn = 1

            gene_rows.append({"sample_id": sid, "call_id": call_id, "cn": cn, "crr": crr_val})

        print(f"  Done {call_id} | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, segs_by_sample
    gc.collect()

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
