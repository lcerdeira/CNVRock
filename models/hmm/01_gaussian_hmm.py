"""
HMM segmentation — version 01: per-sample Gaussian HMM on copy ratios.

Fits a single multi-state Gaussian HMM across all valid contigs of each
sample.  Contigs are treated as independent sequences via hmmlearn's
`lengths` parameter so transitions never cross contig or chromosome
boundaries.

Reads from cfg:
    hmm_n_states          — number of HMM copy-number states
    hmm_self_transition   — diagonal of the transition matrix (stickiness)
    hmm_low_cov_threshold — minimum raw or reconstructed count to include a bin
    hmm_n_jobs            — joblib parallelism (-1 = all cores)
"""

import os
import time

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
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


def fit_hmm_segments(positions, copy_ratios, n_states, self_transition, min_seg_bins=2):
    """Fit a Gaussian HMM to copy ratios for a single chromosome.

    Parameters
    ----------
    positions      : 1-D array of genomic start positions (bp)
    copy_ratios    : 1-D array of copy ratios (input / reconstruction)
    n_states       : int — number of HMM states
    self_transition: float — diagonal of the transition matrix
    min_seg_bins   : int — minimum run length before merging into neighbour

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
# Interactive single-sample fit (used by the diagnostics app)
# ---------------------------------------------------------------------------

def fit_hmm_sample(data, n_states, self_transition, low_cov_threshold):
    """Fit one HMM across all valid contigs of a single sample.

    This is the interactive counterpart to run_hmm_all_samples — it takes a
    pre-processed per-sample DataFrame rather than raw NPY arrays, so it can
    be called directly from the diagnostics app.

    Parameters
    ----------
    data              : DataFrame with columns chrom, start, copy_ratio,
                        input, reconstruction  (output of process_sample)
    n_states          : int — number of HMM states
    self_transition   : float — diagonal of transition matrix
    low_cov_threshold : float — exclude bins with input or reconstruction below this

    Returns
    -------
    DataFrame with columns: chrom, x0, x1, cn, confidence
    """
    from hmmlearn.hmm import GaussianHMM

    contigs = []
    for chrom in data["chrom"].unique():
        sub = data[data["chrom"] == chrom]
        pos = sub["start"].values.astype(float)
        cr  = sub["copy_ratio"].values.astype(float)
        inp = sub["input"].values.astype(float)
        rec = sub["reconstruction"].values.astype(float)

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
        return pd.DataFrame(columns=["chrom", "x0", "x1", "cn", "confidence"])

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
        df.insert(0, "chrom", chrom)
        segs.append(df)

    return pd.concat(segs, ignore_index=True) if segs else pd.DataFrame(
        columns=["chrom", "x0", "x1", "cn", "confidence"]
    )


# ---------------------------------------------------------------------------
# Public entry point (batch)
# ---------------------------------------------------------------------------

def run_hmm_all_samples(store_path, out_dir, cfg):
    """Fit HMM on every sample in parallel and write segments.parquet.

    Config keys read
    ----------------
    hmm_n_states          — number of Gaussian HMM states
    hmm_self_transition   — diagonal of transition matrix (higher = stickier segments)
    hmm_low_cov_threshold — bins with raw or reconstructed count below this are excluded
    hmm_n_jobs            — joblib parallelism (-1 = all cores)

    Loads
    -----
    store_path/contigs.npy        — bin positions: structured array (chrom, start, end)
    store_path/counts.npy         — raw read counts: (n_samples, n_bins)
    out_dir/reconstructions.npy   — denormalised reconstructions: (n_samples, n_bins)
    out_dir/sample_ids.npy        — sample ID strings: (n_samples,)

    Writes
    ------
    out_dir/segments.parquet   — columns: sample_id, chrom, x0, x1, cn, confidence
    """
    import gc
    from joblib import Parallel, delayed

    n_states          = cfg["hmm_n_states"]
    self_transition   = cfg["hmm_self_transition"]
    low_cov_threshold = cfg["hmm_low_cov_threshold"]
    n_jobs            = cfg["hmm_n_jobs"]

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    del contigs
    gc.collect()

    n = len(sample_ids)
    copy_ratios = counts.astype(float) / (recons.astype(float) + 1e-6)   # (n, n_bins)

    args = [
        (sid, copy_ratios[idx], counts[idx].astype(float), recons[idx].astype(float),
         chroms, starts, n_states, self_transition, low_cov_threshold)
        for idx, sid in enumerate(sample_ids)
    ]
    del counts, recons
    gc.collect()

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
        if i % 500 == 0 or i == n:
            elapsed = time.time() - t0
            eta     = (elapsed / i) * (n - i)
            print(f"  HMM {i}/{n} | elapsed {elapsed:.0f}s | eta {eta:.0f}s", flush=True)

    del args, copy_ratios
    gc.collect()

    # Filter out empty DataFrames before concat to avoid FutureWarning on all-NA columns
    results  = [r for r in results if len(r) > 0]
    segments = pd.concat(results, ignore_index=True)
    del results
    gc.collect()

    out_path = os.path.join(out_dir, "segments.parquet")
    segments.to_parquet(out_path, index=False)
    print(f"Saved segments ({len(segments):,} rows) → {out_path}", flush=True)
