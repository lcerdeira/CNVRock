"""
Evaluation — version 05: A. baumannii gene calls vs AMRFinder ground truth.

Ground truth source:
  assets/abaumannii_amrfinder_gt.tsv  (built by data/setup/build_abaumannii_amrfinder_gt.py)

Evaluable genes (absent in a subset of strains → MCC is meaningful):
  blaOXA-23       — acquired OXA-23 carbapenemase (~32% prevalence)
  blaOXA-24-like  — acquired OXA-24/40/72 family   (~3.5%)
  blaOXA-58-like  — acquired OXA-58 family          (~6.5%)

  All three genes are tracked as CHROMOSOMAL loci in gene_coords.tsv.
  Positive call = cn ≥ 1  ("presence mode" — any copy is clinically relevant).

NOT evaluated with binary GT (intrinsic, universal):
  adeB / adeJ — AdeABC/AdeIJK RND efflux pumps; ALL A. baumannii carry these.
                Amplification is the resistance-relevant event; CRR distribution
                is reported but no MCC (no binary absent/present GT).
  adeA/C/R/S/I/K — operon members / regulators, same reason.
  blaOXA-69     — intrinsic OXA-51-like, universal.

Reads from cfg:
  kpsc_gt_path     — path to abaumannii_amrfinder_gt.tsv (required)
  out_dir          — experiment output directory with gene_calls.tsv
  eval_min_group_n — minimum evaluable samples per metric (default 10)
"""

import datetime
import os

import pandas as pd

# Genes evaluated against AMRFinder GT (presence mode: cn ≥ 1 = positive)
PRESENCE_GENES = ["blaOXA-23", "blaOXA-24-like", "blaOXA-58-like"]

# Amplification-report-only genes (CRR distribution reported, no binary GT)
AMP_REPORT_GENES = ["adeB", "adeJ"]

QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
Q_LABELS  = ["p10", "p25", "p50", "p75", "p90"]


def _mcc(tp, tn, fp, fn):
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return (tp * tn - fp * fn) / denom if denom else 0.0


def _metrics(gt_series: pd.Series, pred_series: pd.Series, label: str) -> dict:
    """Compute binary metrics given two integer series (0/1/-1=missing)."""
    df = pd.DataFrame({"gt": gt_series, "pred": pred_series}).dropna()
    df = df[(df["gt"] != -1) & (df["pred"] != -1)]
    n = len(df)
    if n == 0:
        return {"gene": label, "n": 0, "tp": 0, "tn": 0, "fp": 0, "fn": 0,
                "sensitivity": float("nan"), "specificity": float("nan"),
                "ppv": float("nan"), "mcc": float("nan")}
    tp = int(((df["gt"] == 1) & (df["pred"] == 1)).sum())
    tn = int(((df["gt"] == 0) & (df["pred"] == 0)).sum())
    fp = int(((df["gt"] == 0) & (df["pred"] == 1)).sum())
    fn = int(((df["gt"] == 1) & (df["pred"] == 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv  = tp / (tp + fp) if (tp + fp) else float("nan")
    mcc  = _mcc(tp, tn, fp, fn)
    return {"gene": label, "n": n, "tp": tp, "tn": tn, "fp": fp, "fn": fn,
            "sensitivity": round(sens, 4), "specificity": round(spec, 4),
            "ppv": round(ppv, 4), "mcc": round(mcc, 4)}


def run(cfg: dict) -> None:
    """
    Entry point called by the training harness.

    cfg keys used:
      kpsc_gt_path     — AMRFinder GT TSV path (required)
      out_dir          — directory containing gene_calls.tsv
      eval_min_group_n — int (default 10)
    """
    gt_path    = cfg["kpsc_gt_path"]
    out_dir    = cfg["out_dir"]
    min_n      = int(cfg.get("eval_min_group_n", 10))

    calls_path = os.path.join(out_dir, "gene_calls.tsv")
    if not os.path.exists(calls_path):
        print(f"[eval05] gene_calls.tsv not found at {calls_path} — skipping")
        return

    calls = pd.read_csv(calls_path, sep="\t")
    gt    = pd.read_csv(gt_path,    sep="\t", dtype=str)

    # Resolve accession column: some GT files use 'sample_id', others 'accession'
    if "sample_id" not in gt.columns and "accession" in gt.columns:
        gt = gt.rename(columns={"accession": "sample_id"})

    merged = calls.merge(gt, on="sample_id", how="inner", suffixes=("", "_gt"))
    n_total = len(calls)
    n_matched = len(merged)
    print(f"[eval05] gene_calls: {n_total}   matched to GT: {n_matched}")

    lines = []
    lines.append("=" * 68)
    lines.append("CNVRock Evaluation — Acinetobacter baumannii (v05)")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"Samples with GT: {n_matched} / {n_total}")
    lines.append("=" * 68)

    # ── 1. Presence/absence evaluation for acquired OXA genes ─────────────
    lines.append("\n── Acquired OXA carbapenemase presence/absence ──────────────")
    lines.append(f"{'Gene':<20} {'n':>6} {'TP':>5} {'TN':>5} {'FP':>5} {'FN':>5}"
                 f"  {'Sens':>6} {'Spec':>6} {'PPV':>6} {'MCC':>6}")
    lines.append("-" * 68)

    for gene in PRESENCE_GENES:
        if gene not in merged.columns:
            lines.append(f"  {gene}: not in gene_calls.tsv — skipped")
            continue
        if gene + "_gt" not in merged.columns and gene not in gt.columns:
            lines.append(f"  {gene}: GT column missing — skipped")
            continue

        gt_col   = gene + "_gt" if gene + "_gt" in merged.columns else gene
        call_col = gene

        # Predicted: cn >= 1 = present (presence mode)
        pred = merged[call_col].apply(
            lambda x: 1 if x >= 1 else (0 if x == 0 else -1))
        gt_vals = merged[gt_col].astype(int)

        m = _metrics(gt_vals, pred, gene)
        if m["n"] < min_n:
            lines.append(f"  {gene}: n={m['n']} < min_group_n={min_n} — skipped")
            continue
        lines.append(
            f"  {gene:<18} {m['n']:>6} {m['tp']:>5} {m['tn']:>5} "
            f"{m['fp']:>5} {m['fn']:>5}  "
            f"{m['sensitivity']:>6.3f} {m['specificity']:>6.3f} "
            f"{m['ppv']:>6.3f} {m['mcc']:>6.3f}"
        )

    # ── 2. CRR distribution for intrinsic efflux genes ────────────────────
    lines.append("\n── AdeABC / AdeIJK efflux amplification (CRR distribution) ──")
    crr_cols = {g: f"crr_{g}" for g in AMP_REPORT_GENES}
    for gene, crr_col in crr_cols.items():
        if crr_col not in calls.columns:
            lines.append(f"  {crr_col}: not found in gene_calls.tsv — skipped")
            continue
        vals = calls[crr_col].dropna()
        if len(vals) < min_n:
            continue
        qs = vals.quantile(QUANTILES)
        amp_rate = (vals >= 1.75).mean()
        q_str = "  ".join(f"{q_labels}={v:.3f}"
                          for q_labels, v in zip(Q_LABELS, qs))
        lines.append(f"  {gene:<6}  n={len(vals):,}  amp_rate(CRR≥1.75)={amp_rate:.3f}")
        lines.append(f"         {q_str}")

    # ── 3. Uncallable rate ────────────────────────────────────────────────
    lines.append("\n── Uncallable (cn = -1) rates ───────────────────────────────")
    all_tracked = PRESENCE_GENES + AMP_REPORT_GENES
    for gene in all_tracked:
        if gene not in calls.columns:
            continue
        tot = len(calls)
        unc = int((calls[gene] == -1).sum())
        lines.append(f"  {gene:<20}  uncallable: {unc:,} / {tot:,}"
                     f"  ({100*unc/tot:.1f}%)")

    lines.append("\n" + "=" * 68)
    report = "\n".join(lines)
    print(report)

    out_path = os.path.join(out_dir, "evaluation.txt")
    with open(out_path, "w") as fh:
        fh.write(report + "\n")
    print(f"\n[eval05] Wrote {out_path}")
