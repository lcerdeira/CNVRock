"""
Gene-level CNV caller — version 03: bin-count-gated CRR fallback.

Identical to 02 except: the CRR fallback is only applied to genes that span
fewer than cnv_crr_min_bins_fallback bins in the reference panel. Genes with
enough bins are handled by the HMM alone; the fallback is reserved for ultra-
short loci (e.g. GCH1, ~2 bins) where Viterbi is physically unable to
overcome the self_transition penalty regardless of boundary placement.

Rationale: the exp 04 CRR fallback threshold of 1.30 rescued GCH1 and
PM2_PM3 recall but introduced ~190 CRT and ~128 MDR1 false positives by
overriding correct HMM CN=1 calls. These genes have 4 and 7 bins
respectively — enough for the HMM to call them accurately with the
boundary fix. Disabling the fallback for well-covered genes removes the
FPs without losing any TPs, and generalises to the whole genome without
requiring per-gene threshold tuning against labelled data.

Reads from cfg:
    cnv_min_cn1_proportion       — fraction of chromosome bins in CN=1
    cnv_min_confidence           — confidence for CN=1 sanity check
    cnv_flank_padding            — bp padding for CRR flank region
    cnv_crr_amp_threshold        — CRR above which fallback overrides CN=1
    cnv_crr_min_bins_fallback    — only apply fallback when gene spans
                                   fewer than this many bins (default: 3)
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


def call_gene_cnv(data, segments, gene, min_cn1_proportion, min_confidence,
                  flank_padding, crr_amp_threshold, crr_min_bins_fallback):
    """Call copy number for a single gene of interest.

    Parameters
    ----------
    data                  : DataFrame with columns chrom, start, copy_ratio
    segments              : DataFrame with columns chrom, x0, x1, cn, confidence
    gene                  : dict with keys call_id, contig, start, end
    min_cn1_proportion    : fraction of chromosome bins that must be CN=1
    min_confidence        : confidence threshold for the CN=1 sanity check
    flank_padding         : bp upstream / downstream used as flank for CRR
    crr_amp_threshold     : CRR above which fallback overrides CN=1
    crr_min_bins_fallback : fallback only applied when gene spans fewer bins

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
    n_gene_bins = int(gene_mask.sum())
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

    # CRR fallback: only for genes too short for the HMM to resolve reliably
    fallback_eligible = (n_gene_bins < crr_min_bins_fallback)
    if fallback_eligible and result["cn"] in (-1, 1) and result["crr"] is not None \
            and result["crr"] >= crr_amp_threshold:
        result["cn"] = 2

    return result


def call_all_genes(data, segments, min_cn1_proportion, min_confidence,
                   flank_padding, crr_amp_threshold, crr_min_bins_fallback):
    """Call copy number for every gene in GENES_OF_INTEREST.

    Returns a list of dicts (one per gene), suitable for display or batch storage.
    """
    return [
        call_gene_cnv(data, segments, gene,
                      min_cn1_proportion=min_cn1_proportion,
                      min_confidence=min_confidence,
                      flank_padding=flank_padding,
                      crr_amp_threshold=crr_amp_threshold,
                      crr_min_bins_fallback=crr_min_bins_fallback)
        for gene in GENES_OF_INTEREST
    ]


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute gene-level CNV calls for all samples and write gene_calls.tsv.

    Config keys read
    ----------------
    cnv_min_cn1_proportion    — fraction of chromosome bins that must be CN=1
    cnv_min_confidence        — confidence threshold for CN=1 sanity check
    cnv_flank_padding         — bp padding defining the flank region for CRR
    cnv_crr_amp_threshold     — CRR threshold for fallback override
    cnv_crr_min_bins_fallback — apply fallback only when gene has fewer bins

    Loads / Writes — same as 02_gene_cnv_caller.py
    """
    import gc
    import time

    min_cn1_proportion   = cfg["cnv_min_cn1_proportion"]
    min_confidence       = cfg["cnv_min_confidence"]
    flank_padding        = cfg["cnv_flank_padding"]
    crr_amp_threshold    = cfg["cnv_crr_amp_threshold"]
    crr_min_bins_fallback = cfg["cnv_crr_min_bins_fallback"]

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
            if sanity_ok and fallback_eligible and cn in (-1, 1) \
                    and crr_val is not None and crr_val >= crr_amp_threshold:
                cn = 2
            gene_rows.append({"sample_id": sid, "call_id": call_id, "cn": cn, "crr": crr_val})

        print(f"  Done {call_id} | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, counts, segs_by_sample
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
