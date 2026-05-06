"""
Merge Phase B plasmid read counts (blaCTX-M-15, blaNDM-1) into the existing
plasmid NPY store (which currently has only blaKPC-2 / NC_016846.1).

The remap_unmapped_to_plasmids.sh SLURM array job produces per-sample TSV
files in data/inputs/plasmid_remap_counts/, one line each:
    <sample_id>  <ctxm_count>  <ndm_count>

This script reads those files, aligns counts to the existing sample order,
appends two new bins to counts.npy and contigs.npy, and saves in place.
Old files are backed up as *.bak before being replaced.

Usage:
    python data/setup/merge_plasmid_counts.py \\
        --counts-dir data/inputs/plasmid_remap_counts/ \\
        --store-path data/inputs/KpSC-plasmid-1000bp-npy/
"""

import argparse
import os

import numpy as np
import pandas as pd


# New bin definitions — ±500 bp padding around each gene (matches blaKPC-2 convention)
# blaCTX-M-15: MK552109.1:119392-120264
# blaNDM-1:    MZ606384.2:90937-91746
NEW_BINS = [
    {"chrom": "MK552109.1", "start": 118892, "end": 120764},
    {"chrom": "MZ606384.2", "start":  90437, "end":  92246},
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts-dir", required=True,
                        help="Directory of <accession>.plasmid_counts.tsv files "
                             "from remap_unmapped_to_plasmids.sh.")
    parser.add_argument("--store-path", required=True,
                        help="Existing plasmid NPY store directory to update in place.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print diagnostics without writing files.")
    args = parser.parse_args()

    # ── Load existing store ───────────────────────────────────────────────────
    store = args.store_path
    sample_ids = np.load(os.path.join(store, "sample_ids.npy"), allow_pickle=True)
    counts_old = np.load(os.path.join(store, "counts.npy"))
    contigs_old = np.load(os.path.join(store, "contigs.npy"), allow_pickle=True)

    n_samples  = len(sample_ids)
    n_old_bins = counts_old.shape[1]
    print(f"Existing store: {n_samples} samples × {n_old_bins} bins")
    print(f"  Contigs: {list(pd.DataFrame(contigs_old)['chrom'])}")

    # ── Read per-sample remap counts ─────────────────────────────────────────
    sid_to_idx = {str(sid): i for i, sid in enumerate(sample_ids)}
    ctxm_counts = np.zeros(n_samples, dtype=np.int32)
    ndm_counts  = np.zeros(n_samples, dtype=np.int32)

    n_found = 0
    n_missing = 0
    counts_dir = args.counts_dir

    for fname in sorted(os.listdir(counts_dir)):
        if not fname.endswith(".plasmid_counts.tsv"):
            continue
        acc = fname.replace(".plasmid_counts.tsv", "")
        idx = sid_to_idx.get(acc)
        if idx is None:
            continue   # sample not in NPY store (e.g. failed QC)
        fpath = os.path.join(counts_dir, fname)
        with open(fpath) as f:
            line = f.read().strip()
        parts = line.split("\t")
        if len(parts) < 3:
            print(f"  WARNING: malformed line in {fname}: {line!r}")
            n_missing += 1
            continue
        ctxm_counts[idx] = int(float(parts[1]))
        ndm_counts[idx]  = int(float(parts[2]))
        n_found += 1

    # Report samples in store without a remap output file
    for sid, idx in sid_to_idx.items():
        tsv = os.path.join(counts_dir, f"{sid}.plasmid_counts.tsv")
        if not os.path.exists(tsv):
            n_missing += 1

    print(f"\nLoaded remap counts for {n_found}/{n_samples} samples "
          f"({n_missing} missing — will be 0).")

    n_ctxm_pos = int((ctxm_counts > 0).sum())
    n_ndm_pos  = int((ndm_counts > 0).sum())
    ctxm_p50   = int(np.median(ctxm_counts[ctxm_counts > 0])) if n_ctxm_pos else 0
    ndm_p50    = int(np.median(ndm_counts[ndm_counts > 0]))  if n_ndm_pos  else 0
    print(f"  CTX-M: {n_ctxm_pos} samples with counts > 0  (median when present: {ctxm_p50})")
    print(f"  NDM:   {n_ndm_pos} samples with counts > 0   (median when present: {ndm_p50})")

    if args.dry_run:
        print("\nDry run — no files written.")
        return

    # ── Build updated arrays ──────────────────────────────────────────────────
    counts_new = np.concatenate(
        [counts_old, ctxm_counts[:, None], ndm_counts[:, None]], axis=1
    )

    contigs_df = pd.DataFrame(contigs_old)
    new_rows   = pd.DataFrame(NEW_BINS)
    contigs_new = pd.concat([contigs_df, new_rows], ignore_index=True)
    # Convert back to structured numpy array matching original dtype
    contigs_arr = np.array(
        [tuple(row) for _, row in contigs_new.iterrows()],
        dtype=contigs_old.dtype,
    )

    # ── Backup and save ───────────────────────────────────────────────────────
    for fname in ("counts.npy", "contigs.npy"):
        src = os.path.join(store, fname)
        bak = os.path.join(store, fname + ".bak")
        if not os.path.exists(bak):
            import shutil
            shutil.copy2(src, bak)
            print(f"Backed up {fname} → {fname}.bak")

    np.save(os.path.join(store, "counts.npy"),  counts_new)
    np.save(os.path.join(store, "contigs.npy"), contigs_arr)

    print(f"\nUpdated store: {n_samples} samples × {counts_new.shape[1]} bins")
    print(f"  Bins: {list(contigs_new['chrom'])}")
    print(f"  Saved → {store}")


if __name__ == "__main__":
    main()
