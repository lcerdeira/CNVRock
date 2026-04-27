"""
Evaluation — version 04: KpSC gene calls vs AMRFinder+ ground truth.

Ground truth source: NCBI AMRFinder+ (amrfinder --plus) run on AllTheBacteria
KpSC assemblies. The expected TSV has columns:

    sample_id   — matches sample_ids in gene_calls.tsv
    blaSHV      — integer copy count (0 = absent, 1 = single copy, 2+ = amplified)
    ompK35      — integer copy count (0 = absent means gene deleted)
    ompK36      — integer copy count (0 = absent means gene deleted)
    ramR        — integer copy count (0 = absent means gene disrupted/deleted)

Callers produce:
    cn = -1  → uncallable (HMM sanity check failed or low coverage)
    cn =  0  → deletion called (copy number < 1 by HMM)
    cn =  1  → normal (single copy, baseline)
    cn =  2  → amplification called

Ground truth interpretation per gene:
    blaSHV:  AMRFinder count >= 2  → amplified (copy > 1)
             AMRFinder count == 1  → normal (single ancestral copy)
             AMRFinder count == 0  → absent (rare; treat as deletion)
    ompK35:  AMRFinder count == 0  → deleted (loss of function = resistance event)
             AMRFinder count >= 1  → present (normal)
    ompK36:  same interpretation as ompK35
    ramR:    AMRFinder count == 0  → deleted/disrupted
             AMRFinder count >= 1  → present (normal)

Reads from cfg:
    kpsc_gt_path       — path to AMRFinder+ ground truth TSV (required)
    kpsc_meta_path     — path to sample metadata TSV with ST and Species columns
                         (optional; stratified tables omitted if absent)
    eval_min_group_n   — minimum evaluable samples per group to print (default 10)
"""

import datetime
import gc
import os

import pandas as pd


GENES     = ["blaSHV", "ompK35", "ompK36", "ramR"]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
Q_LABELS  = ["p10", "p25", "p50", "p75", "p90"]

# For each gene: what AMRFinder copy count constitutes "positive" (amplification
# or deletion depending on biology)?
# "amp" genes: positive = count >= 2 (extra copies → elevated expression/resistance)
# "del" genes: positive = count == 0 (loss of gene = resistance)
GENE_MODE = {
    "blaSHV": "amp",   # chromosomal, extra copies increase resistance
    "ompK35": "del",   # porin loss → impermeability resistance
    "ompK36": "del",
    "ramR":   "del",   # repressor deletion → efflux upregulation
}


def _amrfinder_to_gt(copy_counts: pd.Series, gene: str) -> pd.Series:
    """Convert AMRFinder copy count to ground-truth label (-1/0/1)."""
    # -1 → missing/unknown in AMRFinder output
    # 0  → normal (not the event of interest)
    # 1  → positive (the resistance-associated copy number state)
    if GENE_MODE[gene] == "amp":
        return copy_counts.map(
            lambda v: -1 if pd.isna(v) else (1 if v >= 2 else 0)
        )
    else:  # "del"
        return copy_counts.map(
            lambda v: -1 if pd.isna(v) else (1 if v == 0 else 0)
        )


def _cn_to_pred(cn: pd.Series, gene: str) -> pd.Series:
    """Convert VAE+HMM CN call to predicted label (-1/0/1)."""
    if GENE_MODE[gene] == "amp":
        return cn.map(lambda v: -1 if v == -1 else (1 if v > 1 else 0))
    else:  # "del"
        return cn.map(lambda v: -1 if v == -1 else (1 if v == 0 else 0))


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
    kpsc_gt_path   = cfg["kpsc_gt_path"]
    kpsc_meta_path = cfg.get("kpsc_meta_path")
    min_group_n    = int(cfg.get("eval_min_group_n", 10))

    wide = (
        pd.read_csv(os.path.join(out_dir, "gene_calls.tsv"), sep="\t")
        .astype({g: pd.Int64Dtype() for g in GENES})
    )

    # Ground truth: AMRFinder+ copy counts per sample and gene
    gt_cols = ["sample_id"] + GENES
    gt = pd.read_csv(kpsc_gt_path, sep="\t", usecols=gt_cols)
    df = gt.merge(wide, on="sample_id")
    del wide, gt
    gc.collect()

    has_meta = False
    if kpsc_meta_path:
        meta = pd.read_csv(kpsc_meta_path, sep="\t", usecols=["sample_id", "Species", "ST"])
        df = meta.merge(df, on="sample_id")
        has_meta = True
        del meta
        gc.collect()

    # Pre-compute labels
    gt_label   = {g: _amrfinder_to_gt(df[g], g)     for g in GENES}
    pred_label = {g: _cn_to_pred(df[g + "_call"] if g + "_call" in df.columns else df[g], g)
                  for g in GENES}
    # gene_calls.tsv uses gene name as column (blaSHV, ompK35, etc.)
    pred_label = {g: _cn_to_pred(df[g], g) for g in GENES}

    gene_results = {}
    crr_results  = {}

    for gene in GENES:
        gt_s   = gt_label[gene]
        pred_s = pred_label[gene]
        crr    = df[f"crr_{gene}"]

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

        crr_results[gene] = {
            "by_pred": {
                label: _crr_quantiles(crr[pred_s == val])
                for val, label in [(-1, "failed"), (0, "pred_normal"), (1, "pred_event")]
            },
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
        "Ground truth: AMRFinder+ copy counts from AllTheBacteria KpSC assemblies.",
        "  blaSHV: positive = extra chromosomal copy (AMRFinder count >= 2).",
        "  ompK35/ompK36: positive = gene absent (AMRFinder count = 0) — deletion.",
        "  ramR: positive = gene absent/disrupted (count = 0) — efflux upregulated.",
        "FNR: fraction of true events called as normal. Primary optimisation target.",
        "  FN p50 >> 1.0 (for amp genes) or << 1.0 (for del genes) = HMM signal present",
        "  but being discarded — try tuning self_transition or HMM state initialisation.",
        "PPV: precision. Low PPV acceptable if FNR is the priority.",
        "  Assembly-derived ground truth is imperfect — apparent FPs may be real.",
        "delta: model-added missingness on GT-callable samples (HMM sanity failures).",
        "  High delta with failed CRR p90 well off 1.0 → rescuable -1 calls.",
        "Note: bacteria are haploid. CN=1 is the normal baseline (same as Pf).",
        "",
        "=" * W,
        "KpSC experiment evaluation",
        f"Generated : {datetime.datetime.utcnow().isoformat()}",
        f"Out dir   : {out_dir}",
        "Ground truth: AMRFinder+ on AllTheBacteria KpSC assemblies",
        "=" * W,
        "",
        "OVERALL",
        "-" * W,
        f"{'Gene':<10} {'Mode':<5} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'call_rate':>10} {'n_eval':>8}",
        "-" * W,
    ]
    for gene in GENES:
        m = gene_results[gene]["overall"]
        call_rate = round(1.0 - (m["pred_missing_rate"] or 0.0), 2)
        lines.append(
            f"{gene:<10} {GENE_MODE[gene]:<5} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
            f"{_fmt(m['ppv']):>5} {_fmt(call_rate):>10} {m['n_eval']:>8}"
        )

    lines += [
        "", "MISSINGNESS", "-" * W,
        f"{'Gene':<10} {'gt_miss':>10} {'pred_miss':>10} {'delta':>8}", "-" * W,
    ]
    for gene in GENES:
        m = gene_results[gene]["overall"]
        lines.append(
            f"{gene:<10} {_fmt(m['gt_missing_rate']):>10} "
            f"{_fmt(m['pred_missing_rate']):>10} {_fmt(m['delta']):>8}"
        )

    q_header = "  " + " ".join(f"{q:>5}" for q in Q_LABELS)
    lines += ["", "CRR BY PREDICTED LABEL  (CRR = gene/flank copy ratio)", "-" * W,
              f"  {'label':<12} {'n':>6}  {q_header.strip()}"]
    for gene in GENES:
        lines.append(f"  — {gene}  ({GENE_MODE[gene]})")
        for label, d in crr_results[gene]["by_pred"].items():
            lines.append(_crr_row(label, d))

    lines += [
        "", "CRR BY CALL OUTCOME  (evaluable samples only)",
        "-" * W,
        f"  {'outcome':<12} {'n':>6}  {q_header.strip()}",
    ]
    for gene in GENES:
        lines.append(f"  — {gene}  ({GENE_MODE[gene]})")
        for label, d in crr_results[gene]["by_outcome"].items():
            lines.append(_crr_row(label, d))

    if any("by_species" in gene_results[g] for g in GENES):
        for gene in GENES:
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

    if any(gene_results[g].get("by_st") for g in GENES):
        for gene in GENES:
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
