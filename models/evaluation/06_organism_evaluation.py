"""
Evaluation — version 06: organism-agnostic presence evaluation.

04_kpsc_evaluation hardcodes the KpSC gene lists, a family->representative map,
and a Kleborate override. A. baumannii needed a bespoke fork (05_). This
version derives everything from the config so E. coli, S. aureus and any later
organism reuse one evaluator:

  - the panel to score is whatever presence columns the AMRFinder GT file
    carries (amrfinder_gt_<organism>.tsv);
  - variant calls are aggregated to the GT's family level using the family
    map (gene_families.tsv), falling back to prefix matching for families
    whose members share a stem (blaKPC -> blaKPC-2; blaCTX-M -> blaCTX-M-*);
  - a family is called present if ANY member variant is called present, which
    matches the read-level MQ=0 aggregation the plasmid caller already does
    and the "blaCTX-M present" convention of clinical reporting.

Ground truth is AMRFinder presence (0/1). Chromosomal intrinsic genes (norA,
ampC, acrB) are absent from the GT by construction and are simply not scored
here — their signal is the CRR distribution, reported elsewhere.

Config keys:
  out_dir                    — results dir with gene_calls.tsv (+ plasmid_gene_calls.tsv)
  kpsc_gt_path               — AMRFinder presence GT (family-level columns)
  plasmid_family_map_path    — optional gene_families.tsv (family <tab> members csv)

Output: out_dir/evaluation.txt and out_dir/evaluation_metrics.tsv
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd


def _resolve(out_dir: str, path: str) -> str:
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(out_dir, path))


def _load_family_map(path: str | None) -> dict[str, list[str]]:
    """family -> [member variant call_ids]; {} if no file."""
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", comment="#")
    col_fam = "family" if "family" in df.columns else df.columns[0]
    col_mem = "members" if "members" in df.columns else df.columns[1]
    out: dict[str, list[str]] = {}
    for _, r in df.iterrows():
        fam = str(r[col_fam]).strip()
        mem = [m.strip() for m in str(r[col_mem]).split(",") if m.strip()]
        if fam:
            out[fam] = mem or [fam]
    return out


def _members_for(family: str, fam_map: dict[str, list[str]],
                 call_cols: set[str]) -> list[str]:
    """Resolve a GT family to the call columns that belong to it."""
    members = set(fam_map.get(family, []))
    # always include an exact-named call column
    if family in call_cols:
        members.add(family)
    # prefix fallback: blaKPC -> blaKPC-2, blaNDM -> blaNDM-1/-5, qnrB -> qnrB1
    stem = family[:-5] if family.endswith("-like") else family
    for c in call_cols:
        if c == stem or c.startswith(stem + "-") or c.startswith(stem):
            # guard against blaSHV matching blaSHV-12 only, not blaSHVxyz
            if c == stem or c[len(stem):len(stem) + 1] in ("", "-") or c.startswith(stem):
                members.add(c)
    # keep only members that are real call columns
    return sorted(m for m in members if m in call_cols)


def _metrics(gt: pd.Series, pred: pd.Series) -> dict:
    from sklearn.metrics import matthews_corrcoef  # noqa: PLC0415
    mask = gt.notna() & pred.notna()
    y_true = (gt[mask] > 0).astype(int)
    y_pred = (pred[mask] > 0).astype(int)
    n = int(mask.sum())
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    npos = int((y_true == 1).sum())
    mcc = (round(float(matthews_corrcoef(y_true, y_pred)), 3)
           if y_true.nunique() > 1 and y_pred.nunique() > 1 else None)
    return {
        "mcc": mcc,
        "fnr": round(fn / npos, 3) if npos else None,
        "ppv": round(tp / (tp + fp), 3) if (tp + fp) else None,
        "n_pos": npos,
        "n_eval": n,
    }


def run_evaluation(out_dir: str, cfg: dict) -> pd.DataFrame:
    gt_path = _resolve(out_dir, cfg["kpsc_gt_path"])
    fam_map = _load_family_map(
        _resolve(out_dir, cfg["plasmid_family_map_path"])
        if cfg.get("plasmid_family_map_path") else None)

    gt = pd.read_csv(gt_path, sep="\t")
    gt_genes = [c for c in gt.columns if c not in ("sample_id", "biosample")]
    # Prefix GT columns so a gene present under the same name in both the calls
    # and the GT (chromosomal blaSHV; mecA/blaZ in S. aureus) does not collide
    # on merge and get silently renamed with _x/_y suffixes.
    gt = gt.rename(columns={g: f"__gt__{g}" for g in gt_genes})

    calls = pd.read_csv(os.path.join(out_dir, "gene_calls.tsv"), sep="\t")
    plasmid_path = os.path.join(out_dir, "plasmid_gene_calls.tsv")
    if os.path.exists(plasmid_path):
        plasmid = pd.read_csv(plasmid_path, sep="\t")
        calls = calls.merge(plasmid, on="sample_id", how="outer")
    else:
        print("plasmid_gene_calls.tsv not found — plasmid genes will be N/A.")

    # binary presence call columns (exclude the crr_/pcn_ magnitude columns)
    call_cols = {c for c in calls.columns
                 if c != "sample_id" and not c.startswith(("crr_", "pcn_"))}

    df = calls.merge(gt, on="sample_id", how="inner")
    rows = []
    for fam in gt_genes:
        truth = df[f"__gt__{fam}"]
        members = _members_for(fam, fam_map, call_cols)
        if not members:
            rows.append({"gene": fam, "mcc": None, "fnr": None, "ppv": None,
                         "n_pos": int((truth > 0).sum()), "n_eval": 0,
                         "note": "no matching call column"})
            continue
        # family present if ANY member variant is called present
        family_pred = (df[members] > 0).any(axis=1).astype(int)
        m = _metrics(truth, family_pred)
        m["gene"] = fam
        m["note"] = "+".join(members) if len(members) > 1 else ""
        rows.append(m)

    metrics = pd.DataFrame(rows)[
        ["gene", "mcc", "fnr", "ppv", "n_pos", "n_eval", "note"]]
    metrics.to_csv(os.path.join(out_dir, "evaluation_metrics.tsv"),
                   sep="\t", index=False)

    lines = ["=" * 64, "ORGANISM EVALUATION (AMRFinder presence)", "=" * 64,
             f"{'gene':<16}{'MCC':>7}{'FNR':>7}{'PPV':>7}{'n_pos':>8}{'n_eval':>8}",
             "-" * 64]
    for _, r in metrics.iterrows():
        def fmt(x): return f"{x:.3f}" if isinstance(x, float) else ("—" if x is None else str(x))
        lines.append(f"{r['gene']:<16}{fmt(r['mcc']):>7}{fmt(r['fnr']):>7}"
                     f"{fmt(r['ppv']):>7}{r['n_pos']:>8}{r['n_eval']:>8}")
    report = "\n".join(lines)
    with open(os.path.join(out_dir, "evaluation.txt"), "w") as fh:
        fh.write(report + "\n")
    print(report)
    return metrics


if __name__ == "__main__":
    import argparse
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    exp_dir = os.path.dirname(os.path.abspath(args.config))
    out_dir = os.path.normpath(os.path.join(exp_dir, cfg["out_dir"]))
    run_evaluation(out_dir, cfg)
