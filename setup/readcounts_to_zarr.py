import zarr
import os
import numpy as np
import pandas as pd

from pathlib import Path
from tqdm.auto import tqdm
from concurrent.futures import ThreadPoolExecutor

PATH_TO_READ_COUNTS = "readcounts"
CONTIGS = [
    "Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3",
    "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3",
    "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3",
    "Pf3D7_13_v3", "Pf3D7_14_v3", "Pf3D7_API_v3", "Pf3D7_MIT_v3"
]
CONTIG_TO_IDX = {c: i for i, c in enumerate(CONTIGS)}

REQUIRED_COLUMNS = {"CONTIG", "START", "END", "COUNT"}

def read_counts_tsv(path: str) -> pd.DataFrame:
    p = Path(path)

    # ── existence & size ──────────────────────────────────────────────────────
    if not p.exists():
        raise FileNotFoundError(f"Readcounts file not found: {path}")
    if p.stat().st_size == 0:
        raise ValueError(f"Readcounts file is empty (0 bytes): {path}")

    # ── parse ─────────────────────────────────────────────────────────────────
    with open(path) as f:
        skip = sum(1 for line in f if line.startswith("@"))

    try:
        df = pd.read_csv(path, sep="\t", skiprows=skip)
    except Exception as exc:
        raise ValueError(f"Failed to parse {path}: {exc}") from exc

    if df.empty:
        raise ValueError(f"Readcounts file parsed to zero rows: {path}")

    # ── schema ────────────────────────────────────────────────────────────────
    missing_cols = REQUIRED_COLUMNS - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Missing columns {missing_cols} in {path} "
            f"(found: {list(df.columns)})"
        )

    df = df[df["CONTIG"].isin(CONTIG_TO_IDX)].copy()
    if df.empty:
        raise ValueError(
            f"No rows matched known contigs in {path}. "
            f"Sample contigs present: {df['CONTIG'].unique().tolist()}"
        )

    df["POS"]        = ((df["START"] + df["END"]) / 2).astype(np.uint32)
    df["CONTIG_IDX"] = df["CONTIG"].map(CONTIG_TO_IDX).astype(np.uint8)
    result = df[["CONTIG_IDX", "POS", "COUNT"]]

    # ── integrity: no NaNs (catches truncated rows) ───────────────────────────
    null_counts = result.isnull().sum()
    if null_counts.any():
        bad = null_counts[null_counts > 0].to_dict()
        raise ValueError(
            f"NaN values found in {path} after parsing — "
            f"file may be truncated. Affected columns: {bad}"
        )

    return result


def validate_all_files(tsv_dir: str, sample_ids: list[str], n_workers: int = 64) -> int:
    """
    Read every file once, check shape consistency, and return the expected bin count.
    Raises a RuntimeError listing ALL bad files rather than stopping at the first.
    """
    print(f"Pre-flight: validating {len(sample_ids)} files …")
    errors: list[str] = []
    bin_counts: dict[str, int] = {}

    def validate_one(sample_id: str):
        path = os.path.join(tsv_dir, f"{sample_id}.counts.tsv")
        try:
            df = read_counts_tsv(path)
            return sample_id, len(df), None
        except Exception as exc:
            return sample_id, None, str(exc)

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for sample_id, n_bins, err in tqdm(
            pool.map(validate_one, sample_ids), total=len(sample_ids)
        ):
            if err:
                errors.append(f"  {sample_id}: {err}")
            else:
                bin_counts[sample_id] = n_bins

    if errors:
        raise RuntimeError(
            f"{len(errors)} file(s) failed validation:\n" + "\n".join(errors)
        )

    # ── all counts must agree ─────────────────────────────────────────────────
    unique_shapes = set(bin_counts.values())
    if len(unique_shapes) > 1:
        from collections import Counter
        shape_freq = Counter(bin_counts.values())
        offenders = [
            f"  {sid}: {n}" for sid, n in bin_counts.items()
            if n != shape_freq.most_common(1)[0][0]
        ]
        raise RuntimeError(
            f"Inconsistent bin counts across files {dict(shape_freq)}.\n"
            f"Outliers:\n" + "\n".join(offenders)
        )

    n_bins = unique_shapes.pop()
    print(f"Pre-flight passed: {len(sample_ids)} files, {n_bins} bins each.")
    return n_bins


def build_store(store_path: str, tsv_dir: str, sample_ids: list[str], n_workers: int = 64):
    # ── validate everything before creating the store ─────────────────────────
    n_bins    = validate_all_files(tsv_dir, sample_ids, n_workers)
    n_samples = len(sample_ids)

    # Re-read first sample for coordinates (already validated above)
    first_df = read_counts_tsv(os.path.join(tsv_dir, f"{sample_ids[0]}.counts.tsv"))

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
        "bin_size":     1000,
        "genome":       "PlasmoDB-54_Pfalciparum3D7",
        "n_bins":       n_bins,
        "sample_ids":   sample_ids,
    })

    def read_one(args):
        idx, sample_id = args
        path = os.path.join(tsv_dir, f"{sample_id}.counts.tsv")
        try:
            df = read_counts_tsv(path)
        except Exception as exc:
            # Re-raise with sample context so the traceback is actionable
            raise RuntimeError(f"Failed reading sample '{sample_id}' ({path}): {exc}") from exc

        if len(df) != n_bins:
            raise RuntimeError(
                f"Bin count changed between validation and write for '{sample_id}': "
                f"expected {n_bins}, got {len(df)}. File may have been modified mid-run."
            )
        return idx, df["COUNT"].values.astype(np.uint32)

    jobs = list(enumerate(sample_ids))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for idx, data in tqdm(pool.map(read_one, jobs), total=n_samples):
            counts[idx, :] = data

    print(f"Done. Store shape: {counts.shape}")
    return store


# ── Run ───────────────────────────────────────────────────────────────────────
tsv_files  = sorted(f for f in os.listdir(PATH_TO_READ_COUNTS) if f.endswith(".counts.tsv"))
sample_ids = [f.replace(".counts.tsv", "") for f in tsv_files]

build_store(
    store_path = "../data/Pf9-53973-samples-1000bp.zarr",
    tsv_dir    = PATH_TO_READ_COUNTS,
    sample_ids = sample_ids,
)