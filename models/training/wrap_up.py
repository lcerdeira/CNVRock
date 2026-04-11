"""
Wrap-up script: runs after a model has finished training.

Steps
-----
1. Inference — encode every sample with mu (deterministic), decode, denormalise.
   Writes: latents.npy, reconstructions.npy, sample_ids.npy
2. HMM segmentation — fit a per-chromosome Gaussian HMM on the copy-ratio
   (raw counts / reconstruction) for every sample.
   Writes: segments.parquet  [sample_id, chrom, x0, x1, cn, confidence]

Usage (standalone — re-run wrap-up on an existing checkpoint):
    python -m training.wrap_up path/to/config.yaml [path/to/checkpoint.pth]

Imported by train.py:
    from training.wrap_up import run_inference, run_hmm_all_samples
"""

import argparse
import json
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader


def _write_json(path, data):
    """Atomically write a JSON file (safe for concurrent readers)."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(data, f)
        os.replace(tmp, path)
    except Exception:
        os.unlink(tmp)
        raise


# ---------------------------------------------------------------------------
# 1. Inference
# ---------------------------------------------------------------------------

def run_inference(model, dataset, device, out_dir, batch_size=128):
    """Encode every sample with mu (deterministic) and save outputs.

    Outputs written to out_dir:
        latents.npy          — (n_samples, latent_dim)  mu vectors
        reconstructions.npy  — (n_samples, n_bins)      raw count space (denormalised)
        sample_ids.npy       — (n_samples,)             sample ID strings
    """
    os.makedirs(out_dir, exist_ok=True)
    model.eval()

    dl = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    all_mu    = []
    all_recon = []

    with torch.no_grad():
        for batch in dl:
            x     = batch.to(device)
            mu, _ = model.enc(x)
            recon = model.dec(mu)                       # deterministic: use mu
            recon_denorm = torch.pow(2, recon) - 1      # inverse of log2(count+1)
            all_mu.append(mu.cpu().numpy())
            all_recon.append(recon_denorm.cpu().numpy())

    latents = np.concatenate(all_mu,    axis=0)         # (n_samples, latent_dim)
    recons  = np.concatenate(all_recon, axis=0)         # (n_samples, n_bins)

    np.save(os.path.join(out_dir, "latents.npy"),         latents)
    np.save(os.path.join(out_dir, "reconstructions.npy"), recons)
    np.save(os.path.join(out_dir, "sample_ids.npy"),      np.array(dataset.sample_ids))

    print(f"Saved latents         {latents.shape} → {out_dir}/latents.npy",       flush=True)
    print(f"Saved reconstructions {recons.shape}  → {out_dir}/reconstructions.npy", flush=True)
    n = len(dataset.sample_ids)
    print(f"Saved sample_ids      ({n},) → {out_dir}/sample_ids.npy",             flush=True)


# ---------------------------------------------------------------------------
# 2. HMM segmentation (pure numpy — no streamlit dependency)
# ---------------------------------------------------------------------------

def _merge_short_runs(states, min_len):
    """Absorb runs shorter than min_len into their preceding (or following) neighbour."""
    states  = states.copy()
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(states):
            j = i
            while j < len(states) and states[j] == states[i]:
                j += 1
            if (j - i) < min_len:
                if i > 0:
                    fill = states[i - 1]
                elif j < len(states):
                    fill = states[j]
                else:
                    fill = states[i]
                states[i:j] = fill
                changed = True
            i = j
    return states


def fit_hmm_segments(positions, copy_ratios, n_states=6, self_transition=0.99, min_seg_bins=3):
    """Fit a Gaussian HMM to copy ratios for a single chromosome.

    Parameters
    ----------
    positions    : 1-D array of genomic start positions (bp)
    copy_ratios  : 1-D array of copy ratios (input / reconstruction)

    Returns
    -------
    DataFrame with columns: x0, x1, cn, confidence
    Returns an empty DataFrame if there are too few observations.
    """
    from hmmlearn.hmm import GaussianHMM

    valid = np.isfinite(copy_ratios) & (copy_ratios > 0)
    pos   = positions[valid]
    obs   = np.clip(copy_ratios[valid], 0, 8).reshape(-1, 1)

    if len(obs) < n_states * 3:
        return pd.DataFrame(columns=["x0", "x1", "cn", "confidence"])

    off      = (1 - self_transition) / (n_states - 1)
    transmat = np.full((n_states, n_states), off)
    np.fill_diagonal(transmat, self_transition)

    model = GaussianHMM(n_components=n_states, covariance_type="diag",
                        n_iter=200, random_state=42,
                        params="smc", init_params="smc")
    model.transmat_ = transmat
    model.fit(obs)

    states     = model.predict(obs)
    posteriors = model.predict_proba(obs)
    means      = model.means_.flatten()

    states      = _merge_short_runs(states, min_len=min_seg_bins)
    state_to_cn = np.clip(np.round(means).astype(int), 0, n_states - 1)

    bin_size = float(np.median(np.diff(pos))) if len(pos) > 1 else 300.0
    cn_seq   = state_to_cn[states]
    bin_conf = posteriors[np.arange(len(states)), states]
    obs_flat = obs.flatten()

    def _segment_confidence(sl):
        hmm_post  = float(bin_conf[sl].mean())
        raw_std   = float(np.std(obs_flat[sl]))
        stability = 1.0 / (1.0 + (raw_std / 0.3) ** 2)
        return hmm_post * stability

    rows, seg_start = [], 0
    for i in range(1, len(cn_seq)):
        if cn_seq[i] != cn_seq[seg_start]:
            rows.append({"x0": pos[seg_start], "x1": pos[i],
                         "cn": int(cn_seq[seg_start]),
                         "confidence": _segment_confidence(slice(seg_start, i))})
            seg_start = i
    rows.append({"x0": pos[seg_start], "x1": pos[-1] + bin_size,
                 "cn": int(cn_seq[seg_start]),
                 "confidence": _segment_confidence(slice(seg_start, None))})

    return pd.DataFrame(rows)


def _hmm_one_sample(args):
    """Fit one HMM across all valid contigs of a single sample.

    Bins are excluded if input OR reconstruction < low_cov_threshold, or if the
    copy ratio is invalid.  Remaining bins are split into contigs wherever the
    gap between consecutive starts exceeds 1.5× the median bin spacing, matching
    the hypervariable-region gaps visible in the plot.  Each contig is treated as
    an independent HMM sequence via hmmlearn's `lengths` parameter so transitions
    never cross contig or chromosome boundaries.
    """
    sid, copy_ratio, raw_counts, recon_counts, chroms, starts, \
        n_states, self_transition, low_cov_threshold = args
    from hmmlearn.hmm import GaussianHMM

    # Build list of (chrom, pos_arr, cr_arr) per valid contig
    contigs = []
    for chrom in np.unique(chroms):
        mask = chroms == chrom
        pos  = starts[mask]
        cr   = copy_ratio[mask]
        inp  = raw_counts[mask]
        rec  = recon_counts[mask]

        valid = (
            np.isfinite(cr) & (cr > 0) &
            (inp >= low_cov_threshold) &
            (rec >= low_cov_threshold)
        )
        pos_v = pos[valid]
        cr_v  = np.clip(cr[valid], 0, 8)

        if len(pos_v) < 2:
            continue

        bin_size   = float(np.median(np.diff(pos_v)))
        breaks     = np.where(np.diff(pos_v) > 1.5 * bin_size)[0] + 1
        idx_groups = np.split(np.arange(len(pos_v)), breaks)

        for idx in idx_groups:
            if len(idx) >= n_states * 3:
                contigs.append((chrom, pos_v[idx], cr_v[idx]))

    if not contigs:
        return pd.DataFrame(columns=["sample_id", "chrom", "x0", "x1", "cn", "confidence"])

    all_obs = np.concatenate([c[2] for c in contigs]).reshape(-1, 1)
    lengths = np.array([len(c[2]) for c in contigs])

    off      = (1 - self_transition) / (n_states - 1)
    transmat = np.full((n_states, n_states), off)
    np.fill_diagonal(transmat, self_transition)

    model = GaussianHMM(n_components=n_states, covariance_type="diag",
                        n_iter=200, random_state=42,
                        params="s", init_params="")
    model.transmat_  = transmat
    model.means_     = np.arange(n_states, dtype=float).reshape(-1, 1)
    model.covars_    = np.full((n_states, 1), 0.25)
    model.startprob_ = np.full(n_states, 1.0 / n_states)
    model.fit(all_obs, lengths=lengths)

    all_states     = model.predict(all_obs, lengths=lengths)
    all_posteriors = model.predict_proba(all_obs, lengths=lengths)

    segs, ptr = [], 0
    for chrom, pos, cr_arr in contigs:
        length   = len(cr_arr)
        states_c = all_states[ptr:ptr + length]
        post_c   = all_posteriors[ptr:ptr + length]
        ptr     += length

        bin_conf = post_c[np.arange(length), states_c]
        bin_size = float(np.median(np.diff(pos))) if length > 1 else 1000.0

        def _seg_conf(sl, _bc=bin_conf, _oc=cr_arr):
            hmm_post  = float(_bc[sl].mean())
            stability = 1.0 / (1.0 + (float(np.std(_oc[sl])) / 0.3) ** 2)
            return hmm_post * stability

        rows, seg_start = [], 0
        for i in range(1, length):
            if states_c[i] != states_c[seg_start]:
                rows.append({"x0": pos[seg_start], "x1": pos[i],
                             "cn": int(states_c[seg_start]),
                             "confidence": _seg_conf(slice(seg_start, i))})
                seg_start = i
        rows.append({"x0": pos[seg_start], "x1": pos[-1] + bin_size,
                     "cn": int(states_c[seg_start]),
                     "confidence": _seg_conf(slice(seg_start, None))})

        df = pd.DataFrame(rows)
        df.insert(0, "chrom",     chrom)
        df.insert(0, "sample_id", sid)
        segs.append(df)

    return pd.concat(segs, ignore_index=True) if segs else pd.DataFrame(
        columns=["sample_id", "chrom", "x0", "x1", "cn", "confidence"]
    )


# ---------------------------------------------------------------------------
# 3. Gene-level CNV calls
# ---------------------------------------------------------------------------

GENES_OF_INTEREST = [
    {"call_id": "MDR1",    "contig": "Pf3D7_05_v3", "start": 955955,  "end": 963095},
    {"call_id": "CRT",     "contig": "Pf3D7_07_v3", "start": 402385,  "end": 406341},
    {"call_id": "GCH1",    "contig": "Pf3D7_12_v3", "start": 974226,  "end": 976097},
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 292244,  "end": 299101},
]


def call_gene_cnv(data, segments, gene, min_cn1_proportion=0.8, min_confidence=0.7,
                  flank_padding=100_000):
    """Call copy number for a single gene of interest.  See utils.py for full docstring."""
    contig  = gene["contig"]
    g_start = gene["start"]
    g_end   = gene["end"]

    result = {
        "call_id":          gene["call_id"],
        "cn":               -1,
        "crr": None,
    }

    chrom_data = data[data["chrom"] == contig]
    if len(chrom_data) == 0:
        return result

    s  = chrom_data["start"].values
    cr = chrom_data["copy_ratio"].values
    gene_mask  = (s >= g_start) & (s <= g_end)
    flank_mask = (s < g_start - flank_padding) | (s > g_end + flank_padding)
    mean_gene  = float(np.nanmean(cr[gene_mask]))  if gene_mask.any()  else None
    mean_flank = float(np.nanmean(cr[flank_mask])) if flank_mask.any() else None
    if mean_gene is not None and mean_flank and mean_flank > 0:
        result["crr"] = mean_gene / mean_flank

    if segments is None or len(segments) == 0:
        return result

    chrom_segs = segments[segments["chrom"] == contig]
    if len(chrom_segs) == 0:
        return result

    # Chromosome sanity: most bins should be CN=1 with high confidence
    high_conf_cn1 = chrom_segs[(chrom_segs["cn"] == 1) & (chrom_segs["confidence"] >= min_confidence)]
    covered = np.zeros(len(s), dtype=bool)
    for _, seg in high_conf_cn1.iterrows():
        covered |= (s >= seg["x0"]) & (s < seg["x1"])
    if covered.sum() / len(s) < min_cn1_proportion:
        return result

    # Single segment must span the entire gene
    spanning = chrom_segs[(chrom_segs["x0"] <= g_start) & (chrom_segs["x1"] >= g_end)]
    if len(spanning) > 0:
        result["cn"] = int(spanning.iloc[0]["cn"])

    return result


def run_hmm_all_samples(store_path, out_dir, n_states=6, self_transition=0.99,
                        low_cov_threshold=10, n_jobs=-1):
    """Fit HMM on every sample/chromosome in parallel and write segments.parquet.

    Loads
    -----
    store_path/contigs.npy        — bin positions: structured array (chrom, start, end)
    store_path/counts.npy         — raw read counts: (n_samples, n_bins)
    out_dir/reconstructions.npy   — denormalised reconstructions: (n_samples, n_bins)
    out_dir/sample_ids.npy        — sample ID strings: (n_samples,)

    Writes
    ------
    out_dir/segments.parquet   — columns: sample_id, chrom, x0, x1, cn, confidence
    out_dir/hmm_progress.json  — live progress for monitoring
    """
    from joblib import Parallel, delayed
    from tqdm import tqdm

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    n      = len(sample_ids)

    copy_ratios  = counts.astype(float) / (recons.astype(float) + 1e-6)   # (n, n_bins)
    progress_path = os.path.join(out_dir, "hmm_progress.json")

    _write_json(progress_path, {"status": "running", "current": 0, "total": n,
                                "elapsed_s": 0.0, "eta_s": None})

    args = [
        (sid, copy_ratios[idx], counts[idx].astype(float), recons[idx].astype(float),
         chroms, starts, n_states, self_transition, low_cov_threshold)
        for idx, sid in enumerate(sample_ids)
    ]

    print(f"Fitting HMM for {n} samples (n_jobs={n_jobs})...", flush=True)
    t0      = time.time()
    results = []

    for i, result in enumerate(
        Parallel(n_jobs=n_jobs, return_as="generator")(
            delayed(_hmm_one_sample)(a) for a in args
        ),
        start=1,
    ):
        results.append(result)
        elapsed = time.time() - t0
        eta     = (elapsed / i) * (n - i)
        _write_json(progress_path, {
            "status":    "running",
            "current":   i,
            "total":     n,
            "elapsed_s": round(elapsed, 1),
            "eta_s":     round(eta, 1),
        })
        if i % 500 == 0 or i == n:
            print(f"  HMM {i}/{n} | elapsed {elapsed:.0f}s | eta {eta:.0f}s", flush=True)

    segments = pd.concat(results, ignore_index=True)
    out_path = os.path.join(out_dir, "segments.parquet")
    segments.to_parquet(out_path, index=False)

    _write_json(progress_path, {
        "status":    "done",
        "current":   n,
        "total":     n,
        "elapsed_s": round(time.time() - t0, 1),
        "eta_s":     0.0,
    })
    print(f"Saved segments ({len(segments):,} rows) → {out_path}", flush=True)

    # Gene-level CNV calls for every sample
    print("Computing gene CNV calls…", flush=True)
    gene_rows = []
    for idx, sid in enumerate(sample_ids):
        # Reconstruct per-sample data DataFrame matching the process_sample schema
        cr = copy_ratios[idx]
        sample_data = pd.DataFrame({
            "chrom":          chroms,
            "start":          starts,
            "copy_ratio":     cr,
        })
        sample_segs = segments[segments["sample_id"] == sid]
        for call in [call_gene_cnv(sample_data, sample_segs, g) for g in GENES_OF_INTEREST]:
            call["sample_id"] = sid
            gene_rows.append(call)

    gene_calls = pd.DataFrame(gene_rows)[
        ["sample_id", "call_id", "cn", "crr"]
    ]
    gene_path = os.path.join(out_dir, "gene_calls.parquet")
    gene_calls.to_parquet(gene_path, index=False)
    print(f"Saved gene calls ({len(gene_calls):,} rows) → {gene_path}", flush=True)


# ---------------------------------------------------------------------------
# Entry point (standalone re-run of wrap-up on an existing checkpoint)
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Re-run inference + HMM wrap-up on an existing checkpoint."
    )
    parser.add_argument("config",     help="Path to experiment config.yaml")
    parser.add_argument("checkpoint", nargs="?",
                        help="Path to checkpoint.pth (default: out_dir/checkpoint.pth)")
    args = parser.parse_args()

    # ── resolve paths ───────────────────────────────────────────────────────
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import yaml  # noqa: PLC0415

    from architectures import ConvVAE  # noqa: PLC0415
    from training.dataset import ReadCountDataset  # noqa: PLC0415

    config_dir = os.path.dirname(os.path.abspath(args.config))
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    def resolve(path):
        return path if os.path.isabs(path) else os.path.join(config_dir, path)

    store_path      = resolve(cfg["store_path"])
    out_dir         = resolve(cfg["out_dir"])
    checkpoint_path = args.checkpoint or os.path.join(out_dir, "checkpoint.pth")

    # ── device ──────────────────────────────────────────────────────────────
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    print(f"Device: {device}", flush=True)

    # ── load model ──────────────────────────────────────────────────────────
    model = ConvVAE(latent_dim=cfg["latent_dim"]).to(device)
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    print(f"Loaded checkpoint: {checkpoint_path}", flush=True)

    # ── inference ───────────────────────────────────────────────────────────
    ds = ReadCountDataset(store_path, normalise=cfg.get("normalise", True))
    run_inference(model, ds, device, out_dir, batch_size=cfg["batch_size"])

    # ── HMM segmentation ────────────────────────────────────────────────────
    print("Fitting HMM segments...", flush=True)
    run_hmm_all_samples(store_path, out_dir)


if __name__ == "__main__":
    main()
