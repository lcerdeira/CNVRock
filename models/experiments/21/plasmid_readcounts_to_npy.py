"""
Build a plasmid NPY store from GATK CollectReadCounts TSVs.

Differences from the chromosomal readcounts_to_npy.py:
  - Includes plasmid contigs (NC_016846.1 for blaKPC-2, plus downloaded
    plasmids for blaCTX-M-15 and blaNDM-1)
  - No zero-fraction filter: plasmid bins are expected to be 0 in samples
    that lack the plasmid — this is signal, not noise.
  - Reads plasmid_readcounts/ not readcounts/

Outputs: ../../../data/inputs/KpSC-plasmid-1000bp-npy/
    counts.npy    — (n_samples, n_plasmid_bins)  raw read counts
    contigs.npy   — (n_plasmid_bins,)            structured array: chrom/start/end
    sample_ids.npy — (n_samples,)

Note: contig names in the TSV must match the CONTIGS list below.
      Verify with: head -30 plasmid_readcounts/SRR*.plasmid_counts.tsv
"""

import os
import numpy as np
import pandas as pd

from pathlib import Path
from tqdm.auto import tqdm
from concurrent.futures import ThreadPoolExecutor

PATH_TO_READ_COUNTS = "plasmid_readcounts"
OUT_DIR             = Path("../../../data/inputs/KpSC-plasmid-1000bp-npy")

# Plasmid contigs to include.
# - NC_016846.1 is the blaKPC-2 plasmid already present in HS11286.fasta.
# - blaCTX-M-15 and blaNDM-1 contig IDs are determined by get_plasmid_references.py
#   and written to assets/plasmid_refs/plasmid_gene_coords.tsv.
#   Update this list after running that script.
#
# To find the actual contig IDs from your TSVs:
#   head -40 plasmid_readcounts/SRR*.plasmid_counts.tsv | grep -v "^@" | awk -F'\t' '{print $1}' | sort -u
CONTIGS_CHROMOSOMAL = ["NC_016845.1"]   # used only to exclude from plasmid store

# Will be populated from the TSV files at runtime — any contig NOT in
# CONTIGS_CHROMOSOMAL is treated as a plasmid contig and included.
# If you want to whitelist specific plasmid contigs, set PLASMID_WHITELIST:
PLASMID_WHITELIST = None  # None = include all non-chromosomal contigs

REQUIRED_COLUMNS = {"CONTIG", "START", "END", "COUNT"}


def read_counts_tsv(path: str, plasmid_contigs: set | None = None) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Not found: {path}")
    if p.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")

    with open(path) as f:
        skip = sum(1 for line in f if line.startswith("@"))

    df = pd.read_csv(path, sep="\t", skiprows=skip)
    if df.empty:
        raise ValueError(f"Zero rows after parsing: {path}")

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")

    # Keep only plasmid contigs (exclude chromosome)
    df = df[~df["CONTIG"].isin(CONTIGS_CHROMOSOMAL)].copy()
    if plasmid_contigs is not None:
        df = df[df["CONTIG"].isin(plasmid_contigs)].copy()

    df = df.rename(columns={"CONTIG": "CHROM"})
    result = df[["CHROM", "START", "END", "COUNT"]].copy()
    result["START"] = result["START"].astype(np.uint32)
    result["END"]   = result["END"].astype(np.uint32)
    result["COUNT"] = result["COUNT"].astype(np.uint32)
    return result


def _discover_plasmid_contigs(tsv_dir: str, sample_ids: list[str]) -> set:
    """Read first 5 files to discover which plasmid contigs are present."""
    contigs: set = set()
    for sid in sample_ids[:5]:
        path = os.path.join(tsv_dir, f"{sid}.plasmid_counts.tsv")
        if not os.path.exists(path):
            continue
        try:
            df = read_counts_tsv(path)
            contigs.update(df["CHROM"].unique())
        except Exception:
            pass
    return contigs


def build_npy(out_dir: Path, tsv_dir: str, sample_ids: list[str], n_workers: int = 16):
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover plasmid contigs from the first few files
    plasmid_contigs = PLASMID_WHITELIST or _discover_plasmid_contigs(tsv_dir, sample_ids)
    if not plasmid_contigs:
        raise RuntimeError(
            "No plasmid contigs found in TSV files. "
            "Check that collect_plasmid_readcounts.sh ran successfully and "
            "that CONTIGS_CHROMOSOMAL is correctly set."
        )
    print(f"Plasmid contigs found: {sorted(plasmid_contigs)}")

    # Read first valid file to get bin layout
    first_df = None
    for sid in sample_ids:
        path = os.path.join(tsv_dir, f"{sid}.plasmid_counts.tsv")
        try:
            df = read_counts_tsv(path, plasmid_contigs)
            if not df.empty:
                first_df = df
                n_bins = len(df)
                break
        except Exception:
            pass
    if first_df is None:
        raise RuntimeError("Could not read any valid plasmid TSV file.")

    print(f"Plasmid bins per sample: {n_bins}")

    contigs_dtype = np.dtype([("chrom", object), ("start", np.uint32), ("end", np.uint32)])
    contigs_arr   = np.empty(n_bins, dtype=contigs_dtype)
    contigs_arr["chrom"] = first_df["CHROM"].values
    contigs_arr["start"] = first_df["START"].values.astype(np.uint32)
    contigs_arr["end"]   = first_df["END"].values.astype(np.uint32)

    np.save(out_dir / "contigs.npy",    contigs_arr)
    np.save(out_dir / "sample_ids.npy", np.array(sample_ids, dtype=object))

    n_samples = len(sample_ids)
    counts    = np.zeros((n_samples, n_bins), dtype=np.uint32)
    missing   = []

    def read_one(args):
        idx, sid = args
        path = os.path.join(tsv_dir, f"{sid}.plasmid_counts.tsv")
        if not os.path.exists(path):
            return idx, None
        try:
            df = read_counts_tsv(path, plasmid_contigs)
            if len(df) != n_bins:
                raise RuntimeError(f"Bin count mismatch: expected {n_bins}, got {len(df)}")
            return idx, df["COUNT"].values.astype(np.uint32)
        except Exception as exc:
            raise RuntimeError(f"Failed reading {sid}: {exc}") from exc

    jobs = list(enumerate(sample_ids))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for idx, data in tqdm(pool.map(read_one, jobs), total=n_samples):
            if data is None:
                missing.append(sample_ids[idx])
            else:
                counts[idx, :] = data

    if missing:
        print(f"WARNING: {len(missing)} samples had no plasmid TSV (will be all-zero rows):")
        for s in missing[:10]:
            print(f"  {s}")

    # No zero-fraction filter — plasmid absence is signal, not noise.
    zero_frac_per_bin = (counts == 0).mean(axis=0)
    n_always_zero = int((zero_frac_per_bin == 1.0).sum())
    print(f"Always-zero bins (no sample has reads): {n_always_zero}/{n_bins}")

    np.save(out_dir / "counts.npy", counts)
    print(f"Done. counts shape: {counts.shape}  →  {out_dir}/counts.npy")
    print(f"  Samples with any plasmid coverage: "
          f"{int((counts.sum(axis=1) > 0).sum())}/{n_samples}")
    return counts


# ── Run ───────────────────────────────────────────────────────────────────────
tsv_files  = sorted(f for f in os.listdir(PATH_TO_READ_COUNTS) if f.endswith(".plasmid_counts.tsv"))
sample_ids = [f.replace(".plasmid_counts.tsv", "") for f in tsv_files]
print(f"Found {len(sample_ids)} plasmid TSV files.")

build_npy(
    out_dir    = OUT_DIR,
    tsv_dir    = PATH_TO_READ_COUNTS,
    sample_ids = sample_ids,
)
