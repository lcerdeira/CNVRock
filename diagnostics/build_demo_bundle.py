#!/usr/bin/env python3
"""
Build a small demo bundle from a full experiment results directory, so the
Streamlit diagnostics app (cnvrock.streamlit.app) can render real outputs
without the gitignored multi-hundred-MB reconstructions.npy.

Subsamples the first N sample IDs (alphabetical) and slices every per-sample
artefact to those IDs:
  latents.npy, reconstructions.npy, sample_ids.npy   (row-sliced by index)
  gene_calls.tsv, plasmid_gene_calls.tsv             (row-filtered by sample_id)
  segments.parquet                                   (filtered by sample_id col)
  evaluation.txt, training_log.json                  (copied verbatim)

Run on HPC (where the full results live):
  python3 diagnostics/build_demo_bundle.py \
      --src  data/results/33_kpsc_expansion_10k \
      --name 33_kpsc_10k_demo \
      --n    200 \
      --label "exp 33 (KpSC 10K)"
"""
from __future__ import annotations
import argparse, os, shutil, json
from pathlib import Path
import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src",  required=True, help="full experiment results dir")
    ap.add_argument("--name", required=True, help="demo bundle subdir name")
    ap.add_argument("--n",    type=int, default=200)
    ap.add_argument("--label", default="")
    ap.add_argument("--out",  default="diagnostics/demo")
    ap.add_argument("--store", default="",
                    help="inputs npy store (counts/contigs) for the CNV viewer")
    ap.add_argument("--showcase", default="",
                    help="semicolon list of 'ID:label' showcase isolates to "
                         "force-include and tag (e.g. 'SRR1:ERG11 4x;SRR2:chr5')")
    a = ap.parse_args()

    src = Path(a.src)
    out = Path(a.out) / a.name
    out.mkdir(parents=True, exist_ok=True)

    # ── Parse showcase isolates (force-included + labelled) ───────────────
    showcase: dict[str, str] = {}
    for tok in a.showcase.split(";"):
        tok = tok.strip()
        if ":" in tok:
            sid, lab = tok.split(":", 1)
            showcase[sid.strip()] = lab.strip()

    # ── Select sample IDs: showcase first, then alphabetical fill to N ────
    ids = np.load(src / "sample_ids.npy", allow_pickle=True)
    ids_str = ids.astype(str)
    sc_present = [s for s in showcase if s in set(ids_str)]
    alpha = sorted(set(ids_str) - set(sc_present))
    fill = alpha[: max(0, a.n - len(sc_present))]
    keep_ids_list = sc_present + fill            # showcase guaranteed in
    pos = {s: i for i, s in enumerate(ids_str)}
    keep_idx = np.array([pos[s] for s in keep_ids_list])
    keep_ids = ids[keep_idx]
    keep_set = set(map(str, keep_ids))
    print(f"{a.name}: {len(ids)} → {len(keep_ids)} samples "
          f"({len(sc_present)} showcase + {len(fill)} fill)")

    # ── Row-sliced numpy arrays ───────────────────────────────────────────
    np.save(out / "sample_ids.npy", keep_ids)
    for fn in ("latents.npy", "reconstructions.npy"):
        if (src / fn).exists():
            arr = np.load(src / fn)
            np.save(out / fn, arr[keep_idx])
            print(f"  {fn}: {arr.shape} → {arr[keep_idx].shape}")

    # ── Row-filtered TSVs (index = sample_id) ─────────────────────────────
    for fn in ("gene_calls.tsv", "plasmid_gene_calls.tsv"):
        if (src / fn).exists():
            df = pd.read_csv(src / fn, sep="\t")
            idcol = "sample_id" if "sample_id" in df.columns else df.columns[0]
            sub = df[df[idcol].astype(str).isin(keep_set)]
            sub.to_csv(out / fn, sep="\t", index=False)
            print(f"  {fn}: {len(df)} → {len(sub)} rows")

    # ── segments.parquet (filter by sample_id column) — optional ──────────
    seg = src / "segments.parquet"
    if seg.exists():
        try:
            s = pd.read_parquet(seg)
            idcol = next((c for c in s.columns if "sample" in c.lower()), None)
            if idcol:
                s = s[s[idcol].astype(str).isin(keep_set)]
            s.to_parquet(out / "segments.parquet")
            print(f"  segments.parquet: {len(s)} rows")
        except Exception as e:
            print(f"  segments.parquet SKIPPED (no parquet engine): {e}")

    # ── Copy small metadata files verbatim ────────────────────────────────
    for fn in ("evaluation.txt", "training_log.json"):
        if (src / fn).exists():
            shutil.copy2(src / fn, out / fn)

    # ── Inputs store slice (counts + contigs) for the CNV-profile viewer ──
    # Sliced in the SAME order as keep_ids so a single sample_ids.npy serves
    # both results (latents/recon) and inputs (counts) via DataFrame .loc.
    if a.store:
        store = Path(a.store)
        st_ids = np.load(store / "sample_ids.npy", allow_pickle=True).astype(str)
        st_pos = {s: i for i, s in enumerate(st_ids)}
        counts = np.load(store / "counts.npy", mmap_mode="r")
        rows, ok_ids = [], []
        for sid in keep_ids.astype(str):
            if sid in st_pos:
                rows.append(np.asarray(counts[st_pos[sid]]))
                ok_ids.append(sid)
        if rows:
            np.save(out / "counts.npy", np.vstack(rows))
            shutil.copy2(store / "contigs.npy", out / "contigs.npy")
            # store-aligned sample_ids for load_inputs (subset that had counts)
            np.save(out / "counts_sample_ids.npy", np.array(ok_ids, dtype=object))
            print(f"  inputs store: {len(rows)}/{len(keep_ids)} samples "
                  f"× {np.vstack(rows).shape[1]} bins")
        else:
            print("  inputs store: no overlapping samples — skipped")

    # ── showcase.json (curated isolates with descriptive labels) ──────────
    sc_kept = {s: showcase[s] for s in keep_ids.astype(str) if s in showcase}
    if sc_kept:
        (out / "showcase.json").write_text(json.dumps(sc_kept, indent=2))
        print(f"  showcase: {len(sc_kept)} curated isolates → "
              f"{list(sc_kept.items())[:3]}")

    # ── README ────────────────────────────────────────────────────────────
    (out / "README.md").write_text(
        f"# Demo: {a.label or a.name} — {len(keep_ids)}-sample subsample\n\n"
        f"First {len(keep_ids)} sample IDs (alphabetical) from "
        f"`{a.src}/`, sliced for the Streamlit diagnostics app. The full "
        f"results directory is gitignored because of reconstructions.npy size.\n\n"
        f"Generated by `diagnostics/build_demo_bundle.py`.\n")

    # report bundle size
    total = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    print(f"  bundle size: {total/1e6:.1f} MB → {out}")


if __name__ == "__main__":
    main()
