import zarr
import os
import numpy as np
import pandas as pd

from pathlib import Path
from tqdm.auto import tqdm
from concurrent.futures import ThreadPoolExecutor

PATH_TO_READ_COUNTS = "/lustre/scratch127/gsu/malariagen/production/runs/cnv/readcounts"
CONTIGS = [
    "Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3",
    "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3",
    "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3",
    "Pf3D7_13_v3", "Pf3D7_14_v3", "Pf3D7_API_v3", "Pf3D7_MIT_v3"
]
CONTIG_TO_IDX = {c: i for i, c in enumerate(CONTIGS)}

def read_counts_tsv(path: str) -> pd.DataFrame:
    with open(path) as f:
        skip = sum(1 for line in f if line.startswith("@"))
    df = pd.read_csv(path, sep="\t", skiprows=skip)
    df = df[df["CONTIG"].isin(CONTIG_TO_IDX)]
    df["POS"]        = ((df["START"] + df["END"]) / 2).astype(np.uint32)
    df["CONTIG_IDX"] = df["CONTIG"].map(CONTIG_TO_IDX).astype(np.uint8)
    return df[["CONTIG_IDX", "POS", "COUNT"]]

def build_store(store_path: str, tsv_dir: str, sample_ids: list[str], n_workers: int = 64):
    # Read first sample to define coordinates
    first_df = read_counts_tsv(os.path.join(tsv_dir, f"{sample_ids[0]}.counts.tsv"))
    n_bins   = len(first_df)
    n_samples = len(sample_ids)

    store = zarr.open(store_path, mode="w")
    store.array("contig_idx", first_df["CONTIG_IDX"].values, dtype="u1", chunks=False)
    store.array("pos",        first_df["POS"].values,        dtype="u4", chunks=False)
    counts = store.zeros(
        "counts",
        shape=(n_samples, n_bins),
        chunks=(100, n_bins),
        dtype="u4",
        compressor=zarr.Blosc(cname="lz4", clevel=5, shuffle=zarr.Blosc.SHUFFLE),
    )
    store.attrs.update({
        "contig_names": CONTIGS,
        "bin_size":     100,
        "genome":       "PlasmoDB-54_Pfalciparum3D7",
        "n_bins":       n_bins,
        "sample_ids":   sample_ids,
    })

    def read_one(args):
        idx, sample_id = args
        df = read_counts_tsv(os.path.join(tsv_dir, f"{sample_id}.counts.tsv"))
        assert len(df) == n_bins, f"Bin count mismatch for {sample_id}: expected {n_bins}, got {len(df)}"
        return idx, df["COUNT"].values.astype(np.uint32)

    jobs = list(enumerate(sample_ids))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for idx, data in tqdm(pool.map(read_one, jobs), total=n_samples):
            counts[idx, :] = data

    print(f"Done. Store shape: {counts.shape}")
    return store


# Run
tsv_files  = sorted(f for f in os.listdir(PATH_TO_READ_COUNTS) if f.endswith(".counts.tsv"))
sample_ids = [f.replace(".counts.tsv", "") for f in tsv_files]

build_store(
    store_path = "../../../data/raw/Pf9-53973-samples-100bp.zarr",
    tsv_dir    = PATH_TO_READ_COUNTS,
    sample_ids = sample_ids,
)