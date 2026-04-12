"""
Gene-level CNV caller — version 03: per-gene CRR fallback thresholds.

Identical to 02 except: cnv_crr_amp_threshold is read per-gene from config
rather than as a single global value. This allows CRT and MDR1 to use higher
thresholds (reducing FPs) while GCH1 and PM2_PM3 keep the lower threshold
needed to catch their weaker signal.

Exp 04 FP/TP CRR distributions motivated the per-gene split:
    CRT:     FP p75=1.58, TP p10=1.64 → clean gap; threshold raised to 1.60
    MDR1:    FP p90=1.48, TP p10=1.51 → clean gap; threshold raised to 1.50
    GCH1:    FP p90=1.40, TP p10=1.37 → overlap;   threshold kept at 1.30
    PM2_PM3: FP p90=1.46, TP p10=1.40 → overlap;   threshold kept at 1.30

Reads from cfg:
    cnv_min_cn1_proportion           — fraction of chromosome bins in CN=1
    cnv_min_confidence               — confidence for CN=1 sanity check
    cnv_flank_padding                — bp padding for CRR flank region
    cnv_crr_amp_threshold_crt        — per-gene CRR override threshold
    cnv_crr_amp_threshold_gch1
    cnv_crr_amp_threshold_mdr1
    cnv_crr_amp_threshold_pm2_pm3
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


def call_gene_cnv(data, segments, gene, min_cn1_proportion, min_confidence, flank_padding, crr_amp_threshold):
    """Call copy number for a single gene of interest.

    Parameters
    ----------
    data               : DataFrame with columns chrom, start, copy_ratio
    segments           : DataFrame with columns chrom, x0, x1, cn, confidence
                         (pre-filtered to this sample)
    gene               : dict with keys call_id, contig, start, end
    min_cn1_proportion : fraction of chromosome bins that must be covered by
                         high-confidence CN=1 to proceed with calling
    min_confidence     : confidence threshold for the CN=1 sanity check
    flank_padding      : bp upstream / downstream of the gene used as flank
    crr_amp_threshold  : gene/flank CRR above which non-amp HMM result is overridden to CN=2

    Returns
    -------
    dict with keys: call_id, cn (-1 = failed), crr (None if uncallable)
    """
    contig  = gene["contig"]
    g_start = gene["start"]
    g_end   = gene["end"]

    result = {"call_id": gene["call_id"], "cn": -1, "crr": None}

    chrom_data = data[data["chrom"] == contig]
    if len(chrom_data) == 0:
        return result

    s  = chrom_data["start"].values
    cr = chrom_data["copy_ratio"].values
    gene_mask  = (s >= g_start) & (s <= g_end)
    flank_mask = (s < g_start - flank_padding) | (s > g_end + flank_padding)
    mean_gene  = float(np.nanmean(cr[gene_mask]))  if gene_mask.any()  else None
    mean_flank = float(np.nanmean(cr[flank_mask])) if flank_mask.any() else None
    if mean_gene is not None and mean_flank and mean_flank > 0:
        result["crr"] = mean_gene / mean_flank

    if segments is None or len(segments) == 0:
        return result

    chrom_segs = segments[segments["chrom"] == contig]
    if len(chrom_segs) == 0:
        return result

    # Chromosome sanity: most bins should be CN=1 with high confidence
    high_conf_cn1 = chrom_segs[
        (chrom_segs["cn"] == 1) & (chrom_segs["confidence"] >= min_confidence)
    ]
    if len(high_conf_cn1) == 0:
        return result
    x0 = high_conf_cn1["x0"].values
    x1 = high_conf_cn1["x1"].values
    covered = ((s[:, None] >= x0) & (s[:, None] < x1)).any(axis=1)
    if covered.sum() / len(s) < min_cn1_proportion:
        return result

    # Single segment must span the entire gene
    spanning = chrom_segs[(chrom_segs["x0"] <= g_start) & (chrom_segs["x1"] >= g_end)]
    if len(spanning) > 0:
        result["cn"] = int(spanning.iloc[0]["cn"])

    if result["cn"] in (-1, 1) and result["crr"] is not None and result["crr"] >= crr_amp_threshold:
        result["cn"] = 2

    return result


def call_all_genes(data, segments, min_cn1_proportion, min_confidence, flank_padding, crr_amp_thresholds):
    """Call copy number for every gene in GENES_OF_INTEREST.

    Parameters
    ----------
    crr_amp_thresholds : dict mapping call_id → threshold, e.g.
                         {"CRT": 1.60, "GCH1": 1.30, "MDR1": 1.50, "PM2_PM3": 1.30}

    Returns a list of dicts (one per gene), suitable for display or batch storage.
    """
    return [
        call_gene_cnv(data, segments, gene,
                      min_cn1_proportion=min_cn1_proportion,
                      min_confidence=min_confidence,
                      flank_padding=flank_padding,
                      crr_amp_threshold=crr_amp_thresholds[gene["call_id"]])
        for gene in GENES_OF_INTEREST
    ]


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute gene-level CNV calls for all samples and write gene_calls.tsv.

    Config keys read
    ----------------
    cnv_min_cn1_proportion           — fraction of chromosome bins that must be CN=1
    cnv_min_confidence               — confidence threshold for CN=1 sanity check
    cnv_flank_padding                — bp padding defining the flank region for CRR
    cnv_crr_amp_threshold_crt        — per-gene CRR fallback thresholds
    cnv_crr_amp_threshold_gch1
    cnv_crr_amp_threshold_mdr1
    cnv_crr_amp_threshold_pm2_pm3

    Loads
    -----
    store_path/contigs.npy        — bin positions: structured array (chrom, start, end)
    store_path/counts.npy         — raw read counts: (n_samples, n_bins)
    out_dir/reconstructions.npy   — denormalised reconstructions: (n_samples, n_bins)
    out_dir/sample_ids.npy        — sample ID strings: (n_samples,)
    out_dir/segments.parquet      — HMM segments from run_hmm_all_samples

    Writes
    ------
    out_dir/gene_calls.tsv — wide format:
        columns: sample_id, CRT, GCH1, MDR1, PM2_PM3, crr_CRT, crr_GCH1, crr_MDR1, crr_PM2_PM3
    """
    import gc
    import time

    min_cn1_proportion = cfg["cnv_min_cn1_proportion"]
    min_confidence     = cfg["cnv_min_confidence"]
    flank_padding      = cfg["cnv_flank_padding"]
    # Per-gene CRR fallback thresholds — keyed by call_id
    crr_amp_thresholds = {
        "CRT":     cfg["cnv_crr_amp_threshold_crt"],
        "GCH1":    cfg["cnv_crr_amp_threshold_gch1"],
        "MDR1":    cfg["cnv_crr_amp_threshold_mdr1"],
        "PM2_PM3": cfg["cnv_crr_amp_threshold_pm2_pm3"],
    }

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
    copy_ratios = counts.astype(float) / (recons.astype(float) + 1e-6)   # (n, n_bins)
    del recons
    gc.collect()

    print(f"Computing gene CNV calls for {n} samples…", flush=True)
    t0 = time.time()

    # Pre-filter segments to the 4 relevant chromosomes, then group by sample
    gene_contigs   = {g["contig"] for g in GENES_OF_INTEREST}
    segs_filtered  = segments[segments["chrom"].isin(gene_contigs)]
    empty_segs     = pd.DataFrame(columns=["chrom", "x0", "x1", "cn", "confidence"])
    segs_by_sample = {sid: grp for sid, grp in segs_filtered.groupby("sample_id")}
    del segments, segs_filtered
    gc.collect()

    gene_rows = []
    for gene in GENES_OF_INTEREST:
        contig              = gene["contig"]
        g_start             = gene["start"]
        g_end               = gene["end"]
        call_id             = gene["call_id"]
        crr_amp_threshold   = crr_amp_thresholds[call_id]

        # Pre-compute bin masks for this gene — identical for every sample
        chrom_mask   = chroms == contig
        s            = starts[chrom_mask]
        gene_mask_l  = (s >= g_start) & (s <= g_end)
        flank_mask_l = (s < g_start - flank_padding) | (s > g_end + flank_padding)
        n_chrom      = int(chrom_mask.sum())

        # CRR: vectorised over all n samples in one shot
        cr_chrom   = copy_ratios[:, chrom_mask]                      # (n, n_chrom_bins)
        mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)   # (n,)
        mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)   # (n,)
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
            if sanity_ok and cn in (-1, 1) and crr_val is not None and crr_val >= crr_amp_threshold:
                cn = 2
            gene_rows.append({"sample_id": sid, "call_id": call_id, "cn": cn, "crr": crr_val})

        print(f"  Done {call_id} | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, counts, segs_by_sample
    gc.collect()

    gene_calls_long = pd.DataFrame(gene_rows)[["sample_id", "call_id", "cn", "crr"]]

    # Pivot to wide format: one row per sample, columns CRT/GCH1/MDR1/PM2_PM3 + crr_*
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
