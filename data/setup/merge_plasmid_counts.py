"""
Merge plasmid remap read counts into the plasmid NPY store.

Reads per-sample TSV files produced by remap_unmapped_to_plasmids.sh and
adds one bin per gene (from plasmid_gene_coords.tsv, excluding blaKPC-2)
to the existing plasmid NPY store.

The store is updated in-place.  Old files are backed up as *.bak on first
run.  To rebuild from scratch, restore the .bak files before re-running:
    cp data/inputs/KpSC-plasmid-1000bp-npy/counts.npy.bak \\
       data/inputs/KpSC-plasmid-1000bp-npy/counts.npy
    cp data/inputs/KpSC-plasmid-1000bp-npy/contigs.npy.bak \\
       data/inputs/KpSC-plasmid-1000bp-npy/contigs.npy

Usage:
    python data/setup/merge_plasmid_counts.py \\
        --counts-dir  data/inputs/plasmid_remap_counts/ \\
        --store-path  data/inputs/KpSC-plasmid-1000bp-npy/ \\
        --coords-tsv  assets/plasmid_refs/plasmid_gene_coords.tsv
"""

import argparse
import csv
import os
import shutil

import numpy as np
import pandas as pd


def load_gene_bins(coords_tsv: str) -> list[dict]:
    """Read non-KPC genes from plasmid_gene_coords.tsv."""
    bins = []
    with open(coords_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["gene"] == "blaKPC-2":
                continue   # already in original store
            bins.append({
                "gene":   row["gene"],
                "chrom":  row["contig"],
                "start":  int(row["start"]),
                "end":    int(row["end"]),
            })
    return bins


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts-dir", required=True,
                        help="Directory of <accession>.plasmid_counts.tsv files.")
    parser.add_argument("--store-path", required=True,
                        help="Plasmid NPY store to update in place.")
    parser.add_argument("--coords-tsv", default="assets/plasmid_refs/plasmid_gene_coords.tsv",
                        help="plasmid_gene_coords.tsv (default: assets/plasmid_refs/plasmid_gene_coords.tsv).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print diagnostics without writing files.")
    parser.add_argument("--restore-bak", action="store_true",
                        help="Restore store from *.bak files before merging "
                             "(use when re-running after adding new genes).")
    args = parser.parse_args()

    store = args.store_path

    # ── Optionally restore from backup ───────────────────────────────────────
    if args.restore_bak:
        for fname in ("counts.npy", "contigs.npy"):
            bak = os.path.join(store, fname + ".bak")
            if os.path.exists(bak):
                shutil.copy2(bak, os.path.join(store, fname))
                print(f"Restored {fname} from {fname}.bak")
            else:
                print(f"WARNING: {fname}.bak not found — cannot restore.")

    # ── Load existing store ───────────────────────────────────────────────────
    sample_ids  = np.load(os.path.join(store, "sample_ids.npy"), allow_pickle=True)
    counts_old  = np.load(os.path.join(store, "counts.npy"))
    contigs_old = np.load(os.path.join(store, "contigs.npy"), allow_pickle=True)

    n_samples  = len(sample_ids)
    n_old_bins = counts_old.shape[1]
    print(f"Existing store: {n_samples} samples × {n_old_bins} bins")
    print(f"  Contigs: {list(pd.DataFrame(contigs_old)['chrom'])}")

    # ── Read gene bins to add from coords TSV ────────────────────────────────
    gene_bins = load_gene_bins(args.coords_tsv)
    if not gene_bins:
        print("No non-KPC genes found in coords TSV — nothing to merge.")
        return
    gene_names = [g["gene"] for g in gene_bins]
    print(f"\nGenes to merge: {gene_names}")

    # ── Check which genes are already in the store ────────────────────────────
    existing_chroms = set(pd.DataFrame(contigs_old)["chrom"].tolist())
    new_bins = [g for g in gene_bins if g["chrom"] not in existing_chroms]
    already  = [g for g in gene_bins if g["chrom"] in existing_chroms]
    if already:
        print(f"  Already in store: {[g['gene'] for g in already]}")
    if not new_bins:
        print("All genes already in store — nothing to add.")
        return
    print(f"  Adding: {[g['gene'] for g in new_bins]}")

    # ── Load per-sample counts from TSV files ─────────────────────────────────
    sid_to_idx = {str(sid): i for i, sid in enumerate(sample_ids)}
    # counts_arr[gene_idx][sample_idx]
    new_counts = np.zeros((n_samples, len(new_bins)), dtype=np.int32)

    n_found = 0
    for fname in sorted(os.listdir(args.counts_dir)):
        if not fname.endswith(".plasmid_counts.tsv"):
            continue
        acc = fname.replace(".plasmid_counts.tsv", "")
        idx = sid_to_idx.get(acc)
        if idx is None:
            continue
        fpath = os.path.join(args.counts_dir, fname)
        with open(fpath) as f:
            lines = f.read().strip().splitlines()
        if len(lines) < 2:
            continue
        header_parts = lines[0].split("\t")   # sample_id gene1 gene2 ...
        value_parts  = lines[1].split("\t")   # acc       cnt1  cnt2  ...
        col_map = {col: i for i, col in enumerate(header_parts)}

        for bi, gbin in enumerate(new_bins):
            gene = gbin["gene"]
            ci   = col_map.get(gene)
            if ci is not None and ci < len(value_parts):
                new_counts[idx, bi] = int(float(value_parts[ci]))
        n_found += 1

    n_missing = n_samples - n_found
    print(f"\nLoaded counts for {n_found}/{n_samples} samples "
          f"({n_missing} missing — will be 0).")
    for bi, gbin in enumerate(new_bins):
        col = new_counts[:, bi]
        n_pos  = int((col > 0).sum())
        median = int(np.median(col[col > 0])) if n_pos else 0
        print(f"  {gbin['gene']}: {n_pos} samples > 0  (median when present: {median})")

    if args.dry_run:
        print("\nDry run — no files written.")
        return

    # ── Build updated arrays ──────────────────────────────────────────────────
    counts_new  = np.concatenate([counts_old, new_counts], axis=1)

    contigs_df  = pd.DataFrame(contigs_old)
    new_rows    = pd.DataFrame([
        {"chrom": g["chrom"], "start": g["start"], "end": g["end"]}
        for g in new_bins
    ])
    contigs_new = pd.concat([contigs_df, new_rows], ignore_index=True)
    contigs_arr = np.array(
        [tuple(row) for _, row in contigs_new.iterrows()],
        dtype=contigs_old.dtype,
    )

    # ── Backup originals (once) and save ─────────────────────────────────────
    for fname in ("counts.npy", "contigs.npy"):
        src = os.path.join(store, fname)
        bak = os.path.join(store, fname + ".bak")
        if not os.path.exists(bak):
            shutil.copy2(src, bak)
            print(f"Backed up {fname} → {fname}.bak")

    np.save(os.path.join(store, "counts.npy"),  counts_new)
    np.save(os.path.join(store, "contigs.npy"), contigs_arr)

    print(f"\nUpdated store: {n_samples} samples × {counts_new.shape[1]} bins")
    print(f"  Bins: {list(contigs_new['chrom'])}")
    print(f"  Saved → {store}")


if __name__ == "__main__":
    main()
