#!/usr/bin/env python3
"""
Select nested stratified subsets of the KpSC expansion cohort for scale studies.

Emits FOUR strictly-nested manifests so 5K ⊂ 10K ⊂ 20K ⊂ 80K — any sample in
the smaller set is also in every larger set. This lets exp 32/33/34/35 isolate
"more training data" as the single variable.

Inputs
------
    assets/kpsc_expansion_kleborate_gt.tsv  — Kleborate v3 ground truth (88,128
                                              rows, identified by `strain` =
                                              BioSample accession)
    assets/kpsc_expansion_metadata.tsv      — Bridges BioSample → run accession
                                              (multi-run BioSamples have a
                                              comma-separated `sample_id` field;
                                              we keep the first run)
    assets/ena_url_manifest.tsv             — FASTQ URLs (one row per run)

Outputs
-------
    assets/kpsc_expansion_subset_5k.tsv      —  5,000 samples
    assets/kpsc_expansion_subset_10k.tsv     — 10,000 samples (superset of 5k)
    assets/kpsc_expansion_subset_20k.tsv     — 20,000 samples (superset of 10k)
    assets/kpsc_expansion_subset_80k.tsv     — all KpSC core samples with FASTQ
                                                URLs (≈77,906 — superset of 20k)
    assets/kpsc_expansion_subset_5k_meta.tsv — per-sample Kleborate annotations
                                                (only the 5k file; bigger ones
                                                share the same metadata bridge)

Stratification (identical to original 5k selector)
--------------------------------------------------
    1. Restrict to KpSC core species: K. pneumoniae, K. quasipneumoniae (both
       subspecies), K. variicola subsp. variicola/tropica, K. africana.
    2. Bridge BioSample → run accession via metadata.sample_id (first run only).
    3. Inner-join with ENA manifest (must have FASTQ URLs).
    4. Stratify by species × Bla_Carb_acquired presence (carbapenemase carrier
       yes/no). Within each stratum, cap any single ST at MAX_PER_ST.
    5. Sample within strata weighted 1.5× toward carbapenemase carriers.

Nesting trick
-------------
We do ONE big weighted shuffle (without replacement) of the full eligible pool
once. The 5K is the first 5,000 of that ordering; the 10K is the first 10,000;
the 20K is the first 20,000; the 80K is everything. By construction, every
smaller set is a prefix (and therefore subset) of every larger set.

Reproducibility: deterministic via fixed random seed (--seed, default 42).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KPSC_CORE_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae subsp. quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. similipneumoniae",
    "Klebsiella variicola subsp. variicola",
    "Klebsiella variicola subsp. tropica",
    "Klebsiella africana",
}

MAX_PER_ST = 150        # cap per ST inside each stratum (avoid overweighting common clones)
CARB_OVERSAMPLE = 1.5   # weight applied to carbapenemase-carrying samples

# Subset sizes in increasing order. Each is a strict superset of the previous.
# None = "everything left after the previous cuts" (full eligible pool).
# The manuscript reports 545 / 5k / 10k / 20k / 40k for the scaling figure;
# 80k is generated for future work / supplementary extension.
SUBSET_TARGETS = [
    ("5k",   5_000),
    ("10k", 10_000),
    ("20k", 20_000),
    ("40k", 40_000),
    ("80k", None),
]

DL_COLS = ["accession", "layout", "r1_url", "r2_url"]
META_COLS = [
    "accession", "species", "ST", "Bla_acquired",
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "AGly_acquired",
    "Flq_acquired", "Tet_acquired", "Sul_acquired",
    "virulence_score", "resistance_score",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kleborate", default="assets/kpsc_expansion_kleborate_gt.tsv")
    ap.add_argument("--metadata",  default="assets/kpsc_expansion_metadata.tsv")
    ap.add_argument("--ena-manifest", default="assets/ena_url_manifest.tsv")
    ap.add_argument("--out-dir",   default="assets")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"Loading Kleborate GT: {args.kleborate}")
    kleb = pd.read_csv(args.kleborate, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(kleb):,} rows")

    print(f"Loading metadata bridge: {args.metadata}")
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    meta["run_acc"] = meta["sample_id"].str.split(",").str[0]
    print(f"  {len(meta):,} rows")

    print(f"Loading ENA URL manifest: {args.ena_manifest}")
    ena = pd.read_csv(args.ena_manifest, sep="\t", dtype=str)
    print(f"  {len(ena):,} rows")

    # ── 1. KpSC core species filter ───────────────────────────────────────────
    in_kpsc = kleb["species"].isin(KPSC_CORE_SPECIES)
    print(f"\nKpSC core species filter:  {in_kpsc.sum():,} / {len(kleb):,} retained")
    kleb = kleb[in_kpsc].copy()

    # ── 2. Bridge BioSample → run accession ───────────────────────────────────
    kleb = kleb.merge(
        meta[["sample_accession", "run_acc"]],
        left_on="strain", right_on="sample_accession", how="inner",
    )
    print(f"After BioSample→run join:  {len(kleb):,}")

    # ── 3. Inner-join with ENA URL manifest ───────────────────────────────────
    merged = kleb.merge(ena, left_on="run_acc", right_on="accession", how="inner")
    print(f"After ENA manifest join:   {len(merged):,} samples have FASTQ URLs")

    # ── 4. Carb flag + ST_filled ──────────────────────────────────────────────
    merged["has_carb"] = (
        merged["Bla_Carb_acquired"].fillna("-").ne("-")
        & merged["Bla_Carb_acquired"].fillna("-").ne("")
    )
    merged["ST_filled"] = merged["ST"].fillna("UNKNOWN").replace("", "UNKNOWN")
    print(f"\nCarbapenemase carriers: "
          f"{merged['has_carb'].sum():,} / {len(merged):,}")

    # ── 5. Per-stratum ST cap → weighted shuffle ──────────────────────────────
    # Build the eligible pool with weights. ST cap prevents any one ST from
    # dominating the prefix; the carb oversample biases earlier slots toward
    # carbapenemase carriers (the positive class CNVRock cares about).
    eligible = []
    eligible_w = []
    for (species, has_carb), group in merged.groupby(["species", "has_carb"]):
        for _, st_group in group.groupby("ST_filled"):
            if len(st_group) > MAX_PER_ST:
                # Deterministically pick MAX_PER_ST from the larger ST groups.
                # These are eligible for ALL subset sizes.
                pick = rng.choice(st_group.index, size=MAX_PER_ST, replace=False)
                eligible.extend(pick)
            else:
                eligible.extend(st_group.index)
            w = CARB_OVERSAMPLE if has_carb else 1.0
            n = min(len(st_group), MAX_PER_ST)
            eligible_w.extend([w] * n)

    eligible = np.array(eligible)
    eligible_w = np.array(eligible_w, dtype=float)
    eligible_w = eligible_w / eligible_w.sum()

    print(f"\nST-capped eligible pool:   {len(eligible):,}")

    # ── Nesting seed: lock in the already-running 5k if present ───────────────
    # We're already running the pipeline against assets/kpsc_expansion_subset_5k.tsv.
    # To avoid wasting any of those count files, every larger subset MUST
    # contain those exact 5,000 samples. We "anchor" by reading them first and
    # then sampling the remainder of the order without replacement.
    out_dir = Path(args.out_dir)
    seed_5k_path = out_dir / "kpsc_expansion_subset_5k.tsv"
    repo_5k_path = Path("assets/kpsc_expansion_subset_5k.tsv")
    anchor_5k = None
    for cand in (seed_5k_path, repo_5k_path):
        if cand.exists():
            anchor_5k = pd.read_csv(cand, sep="\t", dtype=str)
            print(f"\nAnchoring 5k to existing manifest: {cand}  "
                  f"({len(anchor_5k):,} samples)")
            break

    if anchor_5k is not None:
        # Map accession → merged index. Only keep anchors that still survive
        # current filters (ST cap can prune; if the cap shrunk, warn).
        acc_to_idx = dict(zip(merged["accession"], merged.index))
        eligible_set = set(eligible.tolist())
        anchored_indices = []
        for acc in anchor_5k["accession"]:
            idx = acc_to_idx.get(acc)
            if idx is None:
                print(f"  WARNING: anchored sample {acc} no longer in pool", file=sys.stderr)
                continue
            if idx not in eligible_set:
                # Force-include even if ST cap would drop it — anchored
                # samples already exist as count files, we never want to lose them.
                pass
            anchored_indices.append(idx)
        anchored_indices = np.array(anchored_indices)

        # Remainder = eligible \ anchored, weighted-shuffled
        anchored_set = set(anchored_indices.tolist())
        rem_mask = np.array([i not in anchored_set for i in eligible])
        rem_indices = eligible[rem_mask]
        rem_weights = eligible_w[rem_mask]
        rem_weights = rem_weights / rem_weights.sum()
        rem_ordered = rng.choice(rem_indices, size=len(rem_indices),
                                 replace=False, p=rem_weights)
        ordered = np.concatenate([anchored_indices, rem_ordered])
        print(f"  Order = [anchored 5k] + [{len(rem_ordered):,} remainder]")
    else:
        # Fresh run, no anchor.
        ordered = rng.choice(eligible, size=len(eligible),
                             replace=False, p=eligible_w)

    full_pool = ordered  # length == len(eligible) when anchor is empty,
                          # otherwise len(anchor) + (len(eligible) - len(anchor))

    # ── 6. Write nested subsets ───────────────────────────────────────────────
    # The ST-capped `ordered` array has ~38k samples max. Any target >38k
    # (e.g. 40k, 80k) needs to draw extra samples from beyond the ST cap.
    # We build ONE extended ordering: first the ST-capped pool (size 38k),
    # then everything else in random order. Taking prefixes of this extended
    # ordering preserves nesting AND lets larger tiers exceed the cap.
    full_unrestricted_idx = merged.index.to_numpy()
    seen = set(ordered.tolist())
    extras = np.array([i for i in full_unrestricted_idx if i not in seen])
    rng.shuffle(extras)
    full_pool = np.concatenate([ordered, extras])     # full 77,906 samples

    for name, target in SUBSET_TARGETS:
        if target is None:
            indices = full_pool                       # all eligible samples
        else:
            indices = full_pool[:target]              # prefix → strict superset
        subset = merged.loc[indices].copy()

        # Sanity check: every smaller subset must be a subset of every larger.
        if name == "5k":
            prev_set = set()
        cur_set = set(subset["accession"])
        if not prev_set.issubset(cur_set):
            missing = prev_set - cur_set
            raise RuntimeError(
                f"Nesting broken at {name}: {len(missing)} samples in previous "
                f"subset are missing from this one."
            )
        prev_set = cur_set

        print(f"\n── Subset {name}: {len(subset):,} samples ──")
        print("  Species × has_carb:")
        for (sp, hc), n in subset.groupby(["species", "has_carb"]).size().items():
            print(f"    {hc:1}  {n:>6,}  {sp}")
        print(f"  Unique STs: {subset['ST_filled'].nunique():,}")
        print(f"  Carb carriers: {subset['has_carb'].sum():,} "
              f"({100*subset['has_carb'].mean():.1f}%)")

        out_tsv = out_dir / f"kpsc_expansion_subset_{name}.tsv"
        subset[DL_COLS].sort_values("accession").to_csv(
            out_tsv, sep="\t", index=False
        )
        print(f"  → {out_tsv}")

        if name == "5k":
            avail = [c for c in META_COLS if c in subset.columns]
            out_meta = out_dir / "kpsc_expansion_subset_5k_meta.tsv"
            subset[avail].sort_values("accession").to_csv(
                out_meta, sep="\t", index=False
            )
            print(f"  → {out_meta}  (Kleborate cols for the 5k seed set)")

    print("\nAll subsets written. Nesting verified.")


if __name__ == "__main__":
    main()
