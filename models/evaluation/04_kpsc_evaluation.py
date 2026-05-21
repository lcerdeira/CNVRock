"""
Evaluation — version 04: KpSC gene calls vs AMRFinder+ / Kleborate ground truth.

Ground truth source:
  AMRFinder+ (amrfinder --plus) — blaSHV, blaKPC, blaCTX-M, blaNDM, and Phase C
  plasmid genes (qnrB1, blaOXA-48, aac6-Ib-cr).

  When kpsc_kleborate_gt_path is set, the blaSHV column is overridden with
  Kleborate BLAST copy counts (more sensitive for tandem amplification).

NOT evaluated (no read-depth signal):
  ompK35 / ompK36 — loss-of-function is via frameshift/truncation, not physical
    deletion; CRR stays at 1.0. Use Kleborate sequence-level calls.
  ramR — disruption in KpSC is almost exclusively via IS element insertion or
    point mutation; crr_ramR never drops below 0.70 across 545 samples.
    Use Kleborate GT for ramR status.

AMRFinder+ TSV columns (kpsc_gt_path):

    sample_id   — matches sample_ids in gene_calls.tsv / plasmid_gene_calls.tsv
    blaSHV      — integer copy count (0 = absent, 1 = single copy, 2+ = amplified)
    blaKPC      — integer copy count (0 = absent, ≥1 = present on plasmid)
    blaCTX-M    — integer copy count (0 = absent, ≥1 = present on plasmid)
    blaNDM      — integer copy count (0 = absent, ≥1 = present on plasmid)

Callers produce:
    Chromosomal (gene_calls.tsv):
        cn = -1  → uncallable (HMM sanity check failed or low coverage)
        cn =  0  → deletion called
        cn =  1  → normal (single copy, baseline)
        cn =  2  → amplification called
    Plasmid (plasmid_gene_calls.tsv):
        cn = -1  → uncallable (no chromosomal median available)
        cn =  0  → gene absent (PCN < absent_threshold)
        cn =  1  → gene present, single copy
        cn =  2  → gene amplified (PCN ≥ amp_threshold)

Gene modes
----------
    "amp"      — chromosomal; positive = cn ≥ 2   (extra copies increase resistance)
    "del"      — chromosomal; positive = cn = 0   (gene loss = resistance)
    "presence" — plasmid;    positive = cn ≥ 1   (any copy = clinically resistant)

Reads from cfg:
    kpsc_gt_path             — path to AMRFinder+ ground truth TSV (required)
    kpsc_kleborate_gt_path   — path to Kleborate GT TSV (optional); overrides
                               blaSHV column from AMRFinder with BLAST copy counts
    kpsc_meta_path           — path to sample metadata TSV with ST and Species
                               columns (optional)
    eval_min_group_n         — minimum evaluable samples per group (default 10)
    eval_holdout_ids_path    — path to plain-text file of held-out sample IDs
                               (one per line); when set, metrics are computed only
                               on those samples (held-out validation mode)
"""

import datetime
import gc
import os

import pandas as pd


# Chromosomal genes (from gene_calls.tsv, produced by VAE + HMM)
# ompK35, ompK36, ramR excluded: disruption mechanism is sequence-level (frameshift /
# IS insertion / point mutation), not copy-number change — invisible to read-depth CNV.
CHROM_GENES = ["blaSHV"]

# Plasmid accessory genes (from plasmid_gene_calls.tsv, produced by PCN caller)
# GT column names match GENE_PATTERNS keys in get_amrfinder_gt.py
PLASMID_GENES = ["blaKPC", "blaCTX-M", "blaNDM", "qnrB1", "blaOXA-48", "aac6-Ib-cr",
                 "dfrA12", "dfrA14", "sul1", "sul2", "aac3-IIa"]

GENES     = CHROM_GENES + PLASMID_GENES
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
Q_LABELS  = ["p10", "p25", "p50", "p75", "p90"]

# Gene modes:
#   "amp"      — positive = cn ≥ 2  (chromosomal amplification)
#   "del"      — positive = cn = 0  (chromosomal deletion)
#   "presence" — positive = cn ≥ 1  (plasmid gene present, any copy count)
GENE_MODE = {
    "blaSHV":     "amp",       # chromosomal, extra copies increase resistance
    "blaKPC":     "presence",  # plasmid carbapenemase; presence = resistance
    "blaCTX-M":   "presence",  # plasmid ESBL; presence = ESBL phenotype
    "blaNDM":     "presence",  # plasmid metallo-carbapenemase; presence = resistance
    "qnrB1":      "presence",  # plasmid PMQR; presence = quinolone resistance
    "blaOXA-48":  "presence",  # plasmid OXA-type carbapenemase
    "aac6-Ib-cr": "presence",  # plasmid aminoglycoside/fluoroquinolone resistance
    "dfrA12":     "presence",  # plasmid trimethoprim DHFR
    "dfrA14":     "presence",  # plasmid trimethoprim DHFR
    "sul1":       "presence",  # plasmid sulfonamide DHPS
    "sul2":       "presence",  # plasmid sulfonamide DHPS
    "aac3-IIa":   "presence",  # plasmid aminoglycoside (gentamicin); aac(3)-II family
}

# Map GT column name → call column name (plasmid calls use gene name from
# plasmid_gene_coords.tsv, e.g. "blaKPC-2", while GT uses "blaKPC")
PLASMID_GT_TO_CALL = {
    "blaKPC":     "blaKPC-2",
    "blaCTX-M":   "blaCTX-M-15",
    "blaNDM":     "blaNDM-1",
    "qnrB1":      "qnrB1",      # same name in GT and calls
    "blaOXA-48":  "blaOXA-48",
    "aac6-Ib-cr": "aac6-Ib-cr",
    "dfrA12":     "dfrA12",
    "dfrA14":     "dfrA14",
    "sul1":       "sul1",
    "sul2":       "sul2",
    "aac3-IIa":   "aac3-IIa",
}


def _amrfinder_to_gt(copy_counts: pd.Series, gene: str) -> pd.Series:
    """Convert AMRFinder copy count to ground-truth label (-1/0/1).

    Returns:
        -1  missing / unknown
         0  not the event of interest (normal)
         1  positive (the resistance-associated copy number state)
    """
    mode = GENE_MODE[gene]
    if mode == "amp":
        return copy_counts.map(
            lambda v: -1 if pd.isna(v) else (1 if v >= 2 else 0)
        )
    elif mode == "del":
        return copy_counts.map(
            lambda v: -1 if pd.isna(v) else (1 if v == 0 else 0)
        )
    else:  # "presence"
        return copy_counts.map(
            lambda v: -1 if pd.isna(v) else (1 if v >= 1 else 0)
        )


def _cn_to_pred(cn: pd.Series, gene: str) -> pd.Series:
    """Convert CN call to predicted label (-1/0/1)."""
    mode = GENE_MODE[gene]
    if mode == "amp":
        return cn.map(lambda v: -1 if v == -1 else (1 if v > 1 else 0))
    elif mode == "del":
        return cn.map(lambda v: -1 if v == -1 else (1 if v == 0 else 0))
    else:  # "presence"
        return cn.map(lambda v: -1 if v == -1 else (1 if v >= 1 else 0))


def _metrics(gt: pd.Series, pred_gt: pd.Series) -> dict:
    """Compute classification metrics for one (gene, group) slice."""
    from sklearn.metrics import matthews_corrcoef  # noqa: PLC0415
    n         = len(gt)
    eval_mask = (gt != -1) & (pred_gt != -1)
    n_eval    = int(eval_mask.sum())
    m = {
        "n":                  n,
        "gt_missing_rate":    round((gt == -1).sum() / n, 2) if n else None,
        "pred_missing_rate":  round((pred_gt == -1).sum() / n, 2) if n else None,
        "delta":              round(((gt != -1) & (pred_gt == -1)).sum() / n, 2) if n else None,
        "n_eval": n_eval,
        "mcc": None, "fnr": None, "ppv": None,
    }
    if n_eval > 0:
        gt_e, pr_e = gt[eval_mask], pred_gt[eval_mask]
        y_true = (gt_e > 0).astype(int)
        y_pred = (pr_e > 0).astype(int)
        if y_true.nunique() == 2:
            m["mcc"] = round(float(matthews_corrcoef(y_true, y_pred)), 2)
        pos_mask = gt_e == 1
        tp = ((gt_e == 1) & (pr_e == 1)).sum()
        fp = ((gt_e == 0) & (pr_e == 1)).sum()
        if pos_mask.sum() > 0:
            m["fnr"] = round((pr_e[pos_mask] == 0).sum() / pos_mask.sum(), 2)
        if (tp + fp) > 0:
            m["ppv"] = round(float(tp) / float(tp + fp), 2)
    return m


def _crr_quantiles(series: pd.Series) -> dict:
    sub = series.dropna()
    if len(sub) == 0:
        return {"n": 0, **{q: None for q in Q_LABELS}}
    qs = sub.quantile(QUANTILES).round(2).tolist()
    return {"n": len(sub), **dict(zip(Q_LABELS, qs))}


def run_evaluation(out_dir, cfg):
    """Evaluate KpSC gene calls against AMRFinder+ ground truth.

    Config keys read
    ----------------
    kpsc_gt_path     — path to AMRFinder+ TSV (required)
    kpsc_meta_path   — path to sample metadata TSV with Species, ST columns (optional)
    eval_min_group_n — min n_eval for per-ST/species groups (default 10)

    Loads
    -----
    out_dir/gene_calls.tsv   — wide-format calls from run_cnv_calls
    out_dir/segments.parquet — HMM segments for segment diagnostics

    Writes
    ------
    out_dir/evaluation.txt
    """
    kpsc_gt_path          = cfg["kpsc_gt_path"]
    kleborate_gt_path     = cfg.get("kpsc_kleborate_gt_path")
    kpsc_meta_path        = cfg.get("kpsc_meta_path")
    min_group_n           = int(cfg.get("eval_min_group_n", 10))
    holdout_ids_path      = cfg.get("eval_holdout_ids_path")
    holdout_ids: set      = set()
    if holdout_ids_path:
        with open(holdout_ids_path) as _f:
            holdout_ids = {ln.strip() for ln in _f if ln.strip()}
        print(f"Hold-out evaluation mode: restricting metrics to "
              f"{len(holdout_ids)} held-out samples", flush=True)

    # ── Chromosomal calls (VAE + HMM) ───────────────────────────────────────
    chrom_calls = pd.read_csv(os.path.join(out_dir, "gene_calls.tsv"), sep="\t")
    # Only cast gene columns that are actually present (caller version may vary)
    cast_cols = {g: pd.Int64Dtype() for g in CHROM_GENES if g in chrom_calls.columns}
    chrom_calls = chrom_calls.astype(cast_cols)

    # ── Plasmid calls (PCN caller) — optional ────────────────────────────────
    plasmid_path = os.path.join(out_dir, "plasmid_gene_calls.tsv")
    has_plasmid_calls = os.path.exists(plasmid_path)
    if has_plasmid_calls:
        plasmid_calls = pd.read_csv(plasmid_path, sep="\t")
        # Rename blaKPC-2 → blaKPC etc. so they match GT column names
        rename_map = {v: k for k, v in PLASMID_GT_TO_CALL.items()}
        plasmid_calls = plasmid_calls.rename(columns=rename_map)
        pcn_rename = {f"pcn_{v}": f"pcn_{k}" for k, v in PLASMID_GT_TO_CALL.items()}
        plasmid_calls = plasmid_calls.rename(columns=pcn_rename)
        chrom_calls = chrom_calls.merge(plasmid_calls, on="sample_id", how="left")
    else:
        print("plasmid_gene_calls.tsv not found — plasmid gene metrics will be N/A.")

    # ── Ground truth: AMRFinder+ copy counts ────────────────────────────────
    available_genes = CHROM_GENES + (PLASMID_GENES if has_plasmid_calls else [])
    all_gt_cols = pd.read_csv(kpsc_gt_path, sep="\t", nrows=0).columns.tolist()
    gt_cols = ["sample_id"] + [g for g in available_genes if g in all_gt_cols]
    gt = pd.read_csv(kpsc_gt_path, sep="\t", usecols=gt_cols)

    # ── Kleborate GT: override blaSHV if path is set ────────────────────────
    # blaSHV: BLAST copy count from get_kleborate_gt.py overrides AMRFinder
    #         (AMRFinder misses tandem amplifications caught by read depth)
    KLEBORATE_GENES = ["blaSHV"]
    kleborate_gt_path_resolved = (
        kleborate_gt_path
        if kleborate_gt_path and os.path.isabs(kleborate_gt_path)
        else os.path.join(os.path.dirname(kpsc_gt_path), os.path.basename(kleborate_gt_path))
        if kleborate_gt_path else None
    )
    if kleborate_gt_path_resolved and os.path.exists(kleborate_gt_path_resolved):
        kleb_cols_avail = pd.read_csv(kleborate_gt_path_resolved, sep="\t", nrows=0).columns.tolist()
        kleb_use = ["sample_id"] + [g for g in KLEBORATE_GENES if g in kleb_cols_avail]
        kleb_gt = pd.read_csv(kleborate_gt_path_resolved, sep="\t", usecols=kleb_use)
        # Drop Kleborate-covered columns from AMRFinder GT, merge in Kleborate values
        drop_cols = [g for g in KLEBORATE_GENES if g in gt.columns and g in kleb_gt.columns]
        gt = gt.drop(columns=drop_cols).merge(kleb_gt, on="sample_id", how="left")
        print(f"Loaded Kleborate GT → overriding {drop_cols} from AMRFinder.", flush=True)
    elif kleborate_gt_path:
        print(f"kpsc_kleborate_gt_path set but file not found: {kleborate_gt_path_resolved}. "
              "Using AMRFinder GT for blaSHV.", flush=True)

    # Suffix gt columns with _gt to avoid collision with prediction columns
    df = gt.merge(chrom_calls, on="sample_id", suffixes=("_gt", ""))
    del chrom_calls, gt
    gc.collect()

    # Restrict to held-out samples when running hold-out validation
    if holdout_ids:
        n_before = len(df)
        df = df[df["sample_id"].isin(holdout_ids)].reset_index(drop=True)
        print(f"Filtered to {len(df)} / {n_before} held-out samples", flush=True)

    has_meta = False
    if kpsc_meta_path:
        meta = pd.read_csv(kpsc_meta_path, sep="\t", usecols=["sample_id", "Species", "ST"])
        df = meta.merge(df, on="sample_id")
        has_meta = True
        del meta
        gc.collect()

    # Only evaluate genes present in both GT and calls
    eval_genes = [g for g in GENES if f"{g}_gt" in df.columns and g in df.columns]
    missing_genes = [g for g in GENES if g not in eval_genes]
    if missing_genes:
        print(f"Skipping genes not in merged data: {missing_genes}")

    # Pre-compute labels — ground truth uses _gt suffix, predictions use bare gene name
    gt_label   = {g: _amrfinder_to_gt(df[f"{g}_gt"], g) for g in eval_genes}
    pred_label = {g: _cn_to_pred(df[g], g)               for g in eval_genes}

    gene_results = {}
    crr_results  = {}

    for gene in eval_genes:
        gt_s   = gt_label[gene]
        pred_s = pred_label[gene]
        # Chromosomal genes have crr_{gene}; plasmid genes have pcn_{gene}
        if gene in PLASMID_GENES:
            crr = df.get(f"pcn_{gene}", pd.Series(dtype=float, index=df.index))
        else:
            crr = df[f"crr_{gene}"]

        gene_r = {"overall": _metrics(gt_s, pred_s)}

        if has_meta:
            gene_r["by_species"] = {
                sp: _metrics(grp_gt, grp_pred)
                for sp, (grp_gt, grp_pred) in {
                    sp: (gt_s[df["Species"] == sp], pred_s[df["Species"] == sp])
                    for sp in df["Species"].dropna().unique()
                }.items()
            }
            st_results = {}
            for st, grp in df.groupby("ST"):
                m = _metrics(gt_s[grp.index], pred_s[grp.index])
                if m["n_eval"] >= min_group_n:
                    st_results[str(st)] = m
            gene_r["by_st"] = st_results

        gene_results[gene] = gene_r

        crr_label = "pcn" if gene in PLASMID_GENES else "crr"
        crr_results[gene] = {
            "by_pred": {
                label: _crr_quantiles(crr[pred_s == val])
                for val, label in [(-1, "failed"), (0, "pred_normal"), (1, "pred_event")]
            },
            "_label": crr_label,
        }
        eval_mask = (gt_s != -1) & (pred_s != -1)
        crr_results[gene]["by_outcome"] = {
            label: _crr_quantiles(crr[eval_mask & (gt_s == tv) & (pred_s == pv)])
            for (tv, pv), label in [
                ((0, 0), "TN"), ((0, 1), "FP"), ((1, 0), "FN"), ((1, 1), "TP")
            ]
        }

    del df
    gc.collect()

    seg_diag = _segment_diagnostics(out_dir)

    # ── Format text report ───────────────────────────────────────────────────
    W = 64

    def _fmt(v):
        return f"{v:.2f}" if isinstance(v, float) else ("N/A" if v is None else str(v))

    def _crr_row(label, d):
        qs = " ".join(f"{_fmt(d[q]):>5}" for q in Q_LABELS)
        return f"  {label:<12} {d['n']:>6}  {qs}"

    lines = [
        "=" * W,
        "GUIDANCE",
        "-" * W,
        "Ground truth: AMRFinder+ (blaSHV, blaKPC, blaCTX-M, blaNDM, qnrB1,",
        "  blaOXA-48, aac6-Ib-cr); blaSHV overridden with Kleborate BLAST copy",
        "  counts if kpsc_kleborate_gt_path is set.",
        "  blaSHV: positive = extra chromosomal copy (count >= 2).",
        "  blaKPC/blaCTX-M/blaNDM/qnrB1/blaOXA-48/aac6-Ib-cr: positive = present.",
        "FNR: fraction of true events called as normal. Primary optimisation target.",
        "  FN CRR/PCN p50 >> 1.0 (amp/presence) or << 1.0 (del) = signal present but",
        "  being discarded — try tuning thresholds or HMM state initialisation.",
        "PPV: precision. Low PPV acceptable if FNR is the priority.",
        "  Assembly-derived ground truth is imperfect — apparent FPs may be real.",
        "delta: model-added missingness on GT-callable samples (HMM sanity failures).",
        "Note: bacteria are haploid. CN=1 is the normal chromosomal baseline.",
        "",
        "=" * W,
        "KpSC experiment evaluation",
        f"Generated : {datetime.datetime.utcnow().isoformat()}",
        f"Out dir   : {out_dir}",
        f"Ground truth: AMRFinder+ (blaSHV, plasmid genes)"
        + (f" + Kleborate (blaSHV copy counts)" if kleborate_gt_path_resolved and os.path.exists(kleborate_gt_path_resolved or "") else ""),
        "=" * W,
        "",
        "OVERALL",
        "-" * W,
        f"{'Gene':<12} {'Mode':<9} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'call_rate':>10} {'n_eval':>8}",
        "-" * W,
    ]
    for gene in eval_genes:
        m = gene_results[gene]["overall"]
        call_rate = round(1.0 - (m["pred_missing_rate"] or 0.0), 2)
        lines.append(
            f"{gene:<12} {GENE_MODE[gene]:<9} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
            f"{_fmt(m['ppv']):>5} {_fmt(call_rate):>10} {m['n_eval']:>8}"
        )

    lines += [
        "", "MISSINGNESS", "-" * W,
        f"{'Gene':<12} {'gt_miss':>10} {'pred_miss':>10} {'delta':>8}", "-" * W,
    ]
    for gene in eval_genes:
        m = gene_results[gene]["overall"]
        lines.append(
            f"{gene:<12} {_fmt(m['gt_missing_rate']):>10} "
            f"{_fmt(m['pred_missing_rate']):>10} {_fmt(m['delta']):>8}"
        )

    q_header = "  " + " ".join(f"{q:>5}" for q in Q_LABELS)
    lines += [
        "", "CRR/PCN BY PREDICTED LABEL",
        "  (chromosomal genes: CRR = gene/flank copy ratio; plasmid genes: PCN = depth/chrom_median)",
        "-" * W,
        f"  {'label':<12} {'n':>6}  {q_header.strip()}",
    ]
    for gene in eval_genes:
        metric_lbl = crr_results[gene]["_label"].upper()
        lines.append(f"  — {gene}  ({GENE_MODE[gene]})  [{metric_lbl}]")
        for label, d in crr_results[gene]["by_pred"].items():
            lines.append(_crr_row(label, d))

    lines += [
        "", "CRR/PCN BY CALL OUTCOME  (evaluable samples only)",
        "-" * W,
        f"  {'outcome':<12} {'n':>6}  {q_header.strip()}",
    ]
    for gene in eval_genes:
        metric_lbl = crr_results[gene]["_label"].upper()
        lines.append(f"  — {gene}  ({GENE_MODE[gene]})  [{metric_lbl}]")
        for label, d in crr_results[gene]["by_outcome"].items():
            lines.append(_crr_row(label, d))

    if any("by_species" in gene_results[g] for g in eval_genes):
        for gene in eval_genes:
            if "by_species" not in gene_results[gene]:
                continue
            lines += [
                "", f"BY SPECIES — {gene}  ({GENE_MODE[gene]})", "-" * W,
                f"  {'Species':<22} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'n_eval':>8}",
                f"  {'-' * (W - 2)}",
            ]
            for sp, m in sorted(gene_results[gene]["by_species"].items()):
                lines.append(
                    f"  {sp:<22} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
                    f"{_fmt(m['ppv']):>5} {m['n_eval']:>8}"
                )

    if any(gene_results[g].get("by_st") for g in eval_genes):
        for gene in eval_genes:
            st_r = gene_results[gene].get("by_st", {})
            if not st_r:
                continue
            lines += [
                "", f"BY SEQUENCE TYPE — {gene}  (n_eval ≥ {min_group_n})", "-" * W,
                f"  {'ST':<12} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'n_eval':>8}",
                f"  {'-' * (W - 2)}",
            ]
            for st, m in sorted(st_r.items()):
                lines.append(
                    f"  {st:<12} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
                    f"{_fmt(m['ppv']):>5} {m['n_eval']:>8}"
                )

    if seg_diag is not None:
        trans = seg_diag["transitions_percentiles"]
        trans_str = "  ".join(
            f"{lbl}={trans[lbl]:.1f}" for lbl in ["p10", "p25", "p50", "p75", "p90"]
        )
        lines += [
            "", "SEGMENT DIAGNOSTICS", "-" * W,
            "  (K. pneumoniae is haploid — CN=1 is the normal baseline.)",
            f"  callability (samples with any non-CN1 segment): {seg_diag['callability']:.3f}",
            f"  within-chrom CN transitions per sample (percentiles):",
            f"    {trans_str}",
            f"  n_samples analysed: {seg_diag['n_samples']}",
        ]
    else:
        lines += ["", "SEGMENT DIAGNOSTICS", "-" * W,
                  "  segments.parquet not found — skipped."]

    out_path = os.path.join(out_dir, "evaluation.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved evaluation → {out_path}", flush=True)
    return {"genes": gene_results, "crr": crr_results}


def _segment_diagnostics(out_dir):
    """Return callability and CN transition percentiles from segments.parquet."""
    seg_path = os.path.join(out_dir, "segments.parquet")
    if not os.path.exists(seg_path):
        return None

    segs = pd.read_parquet(seg_path, columns=["sample_id", "chrom", "x0", "cn"])

    sample_has_cnv = segs.groupby("sample_id")["cn"].apply(lambda x: (x != 1).any())
    callability    = float(sample_has_cnv.mean())
    n_samples      = len(sample_has_cnv)

    segs_s = segs.sort_values(["sample_id", "chrom", "x0"])
    segs_s["prev_cn"]  = segs_s.groupby(["sample_id", "chrom"])["cn"].shift(1)
    segs_s["is_trans"] = (
        segs_s["prev_cn"].notna() & (segs_s["cn"] != segs_s["prev_cn"])
    )
    trans_per_sample = segs_s.groupby("sample_id")["is_trans"].sum()
    q_labels = ["p10", "p25", "p50", "p75", "p90"]
    trans_pcts = dict(zip(
        q_labels,
        trans_per_sample.quantile([0.10, 0.25, 0.50, 0.75, 0.90]).round(1).tolist(),
    ))

    del segs, segs_s
    gc.collect()

    return {
        "callability":             round(callability, 3),
        "transitions_percentiles": trans_pcts,
        "n_samples":               n_samples,
    }
