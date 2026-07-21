"""
Plasmid gene CNV caller — version 07.

Computes Plasmid Copy Number (PCN) for accessory resistance genes carried on
plasmids:
    blaKPC-2    — carbapenemase; any copy = clinically resistant
    blaCTX-M-15 — ESBL; any copy = ESBL phenotype
    blaNDM-1    — metallo-carbapenemase; any copy = pan-resistant

Architecture
------------
Plasmid genes are NOT sent through the VAE/HMM pipeline. Instead, PCN is
computed directly from raw read depth:

    PCN = mean(plasmid_gene_bins) / chromosomal_sample_median

where chromosomal_sample_median is the per-sample median depth across the
chromosome (stored in the chromosomal NPY store as _medians).

Interpretation
--------------
PCN thresholds (tunable via config):
    < pcn_absent_threshold          → gene absent  (cn = 0)
    pcn_absent_threshold to         → 1 copy        (cn = 1)
    pcn_absent_threshold   pcn_amp  → amplified     (cn = 2)
    ≥ pcn_amp_threshold             → amplified     (cn = 2)

For plasmid genes, cn = 0 means the sample lacks the plasmid (normal for most
KpSC). cn = 1 means one plasmid copy. cn = 2 means amplified (multiple copies
of the plasmid or tandem duplication of the gene).

The "positive" clinical event differs from chromosomal genes:
    chromosomal: amplification (cn ≥ 2) or deletion (cn = 0)
    plasmid:     presence (cn ≥ 1) is the resistance-conferring event

Config keys
-----------
    plasmid_store_path       — path to KpSC-plasmid-1000bp-npy/
    plasmid_gene_coords_path — path to assets/plasmid_refs/plasmid_gene_coords.tsv
    pcn_absent_threshold     — PCN below this → gene absent  (default 0.2)
    pcn_amp_threshold        — PCN above this → amplified     (default 1.5)

Per-gene threshold override
---------------------------
    plasmid_gene_coords.tsv may include an optional `absent_threshold` column.
    When present for a gene, it overrides pcn_absent_threshold for that gene only.
    Useful when one plasmid reference has more cross-mapping noise than others
    (e.g. blaCTX-M-15 reference attracts reads from other CTX-M family members).

Outputs
-------
    out_dir/plasmid_gene_calls.tsv
        Columns: sample_id, blaKPC-2, blaCTX-M-15, blaNDM-1,
                 pcn_blaKPC-2, pcn_blaCTX-M-15, pcn_blaNDM-1
        cn values: -1 = uncallable, 0 = absent, 1 = single copy, 2 = amplified
"""

import csv
import gc
import os
import time

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Calling logic
# ---------------------------------------------------------------------------

def run_plasmid_cnv_calls(out_dir: str, cfg: dict) -> None:
    """Compute plasmid gene PCN and CN calls for all samples."""
    plasmid_store   = cfg["plasmid_store_path"]
    coords_path     = cfg["plasmid_gene_coords_path"]
    absent_thresh   = float(cfg.get("pcn_absent_threshold", 0.20))
    amp_thresh      = float(cfg.get("pcn_amp_threshold",    1.50))
    chrom_store     = cfg["store_path"]   # chromosomal store for sample medians

    # ── Load plasmid store ───────────────────────────────────────────────────
    plasmid_counts  = np.load(os.path.join(plasmid_store, "counts.npy"))       # (n, p_bins)
    plasmid_contigs = pd.DataFrame(
        np.load(os.path.join(plasmid_store, "contigs.npy"), allow_pickle=True)
    )
    plasmid_ids     = np.load(os.path.join(plasmid_store, "sample_ids.npy"), allow_pickle=True)

    # ── Load chromosomal medians for PCN normalisation ───────────────────────
    medians_path = os.path.join(chrom_store, "medians.npy")
    if os.path.exists(medians_path):
        chrom_medians = np.load(medians_path)               # (n_chrom_samples,)
        chrom_ids     = np.load(
            os.path.join(chrom_store, "sample_ids.npy"), allow_pickle=True
        )
        median_map = dict(zip(chrom_ids, chrom_medians))
    else:
        # Compute from chromosomal counts on the fly (slower)
        chrom_counts = np.load(os.path.join(chrom_store, "counts.npy"))
        chrom_ids    = np.load(
            os.path.join(chrom_store, "sample_ids.npy"), allow_pickle=True
        )
        chrom_medians_arr = np.median(chrom_counts.astype(np.float32), axis=1)
        median_map = dict(zip(chrom_ids, chrom_medians_arr))
        del chrom_counts
        gc.collect()

    # ── Load gene coordinates ────────────────────────────────────────────────
    genes = []
    with open(coords_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Per-gene absent_threshold overrides global config value when set
            gene_absent = row.get("absent_threshold", "").strip()
            genes.append({
                "gene":             row["gene"],
                "contig":           row["contig"],
                "start":            int(row["start"]),
                "end":              int(row["end"]),
                "absent_threshold": float(gene_absent) if gene_absent else absent_thresh,
            })

    if not genes:
        raise RuntimeError(f"No gene coordinates found in {coords_path}")

    print(f"Plasmid CNV calling: {len(genes)} genes × {len(plasmid_ids)} samples", flush=True)
    t0 = time.time()

    p_chroms = plasmid_contigs["chrom"].values
    p_starts = plasmid_contigs["start"].values.astype(float)
    p_ends   = plasmid_contigs["end"].values.astype(float)

    rows = []
    for gene_info in genes:
        contig  = gene_info["contig"]
        g_start = gene_info["start"]
        g_end   = gene_info["end"]
        gene    = gene_info["gene"]

        # Match this gene's exact (contig, start, end) store column, not just
        # the contig. When several genes share a plasmid contig — e.g.
        # blaZ/blaR1/blaI/ermB all on pI258 in S. aureus — a contig-only mask
        # would pool them and hand every gene the same whole-contig depth
        # (blaZ at 72 % carriage and ermB at 1 % would then look identical).
        # Each store column is already a per-gene aggregate, so the exact
        # coordinate selects the right one. For one-gene-per-contig references
        # (KpSC) the contig appears once and every start/end matches, so this is
        # identical to the previous contig-only mask.
        gene_mask = ((p_chroms == contig)
                     & (p_starts == g_start) & (p_ends == g_end))
        n_gene_bins = int(gene_mask.sum())

        if n_gene_bins == 0:
            # Contig not present in plasmid store — this gene not yet available
            print(f"  WARNING: no bins for {gene} on {contig} — filling cn=-1", flush=True)
            for sid in plasmid_ids:
                rows.append({"sample_id": sid, "gene": gene, "cn": -1, "pcn": None})
            continue

        gene_counts = plasmid_counts[:, gene_mask].astype(float)  # (n, n_gene_bins)
        mean_depth  = gene_counts.mean(axis=1)                     # (n,)

        gene_absent_thresh = gene_info["absent_threshold"]
        for idx, sid in enumerate(plasmid_ids):
            chrom_med = median_map.get(sid)
            if chrom_med is None or chrom_med < 1.0:
                cn  = -1
                pcn = None
            else:
                pcn = float(mean_depth[idx]) / float(chrom_med)
                if pcn < gene_absent_thresh:
                    cn = 0
                elif pcn >= amp_thresh:
                    cn = 2
                else:
                    cn = 1
            rows.append({"sample_id": sid, "gene": gene, "cn": cn,
                         "pcn": round(pcn, 4) if pcn is not None else None})

        n_absent   = sum(1 for r in rows if r["gene"] == gene and r["cn"] == 0)
        n_present  = sum(1 for r in rows if r["gene"] == gene and r["cn"] == 1)
        n_amp      = sum(1 for r in rows if r["gene"] == gene and r["cn"] == 2)
        n_unc      = sum(1 for r in rows if r["gene"] == gene and r["cn"] == -1)
        print(
            f"  {gene}: absent={n_absent}  1-copy={n_present}  "
            f"amplified={n_amp}  uncallable={n_unc}  "
            f"[{time.time()-t0:.0f}s]",
            flush=True,
        )

    del plasmid_counts
    gc.collect()

    # ── Wide format output ───────────────────────────────────────────────────
    gene_names  = [g["gene"] for g in genes]
    long_df     = pd.DataFrame(rows)
    cn_wide     = long_df.pivot(index="sample_id", columns="gene", values="cn")
    pcn_wide    = long_df.pivot(index="sample_id", columns="gene", values="pcn")
    cn_wide.columns.name  = None
    pcn_wide.columns      = [f"pcn_{g}" for g in pcn_wide.columns]
    pcn_wide.columns.name = None
    result = pd.concat([cn_wide, pcn_wide], axis=1).reset_index()

    out_path = os.path.join(out_dir, "plasmid_gene_calls.tsv")
    result.to_csv(out_path, sep="\t", index=False)
    print(f"Saved plasmid gene calls ({len(result):,} samples) → {out_path}", flush=True)
