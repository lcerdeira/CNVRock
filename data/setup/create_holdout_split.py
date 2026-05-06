"""
Create a stratified 80/20 hold-out split for KpSC validation.

Stratification is done by a binary presence vector across the three highest-
prevalence genes (blaKPC, blaCTX-M, blaNDM), giving up to 8 strata. Within
each stratum 20% of samples are randomly assigned to the hold-out set.

Outputs:
    assets/holdout_sample_ids.txt  — one sample ID per line (held-out 20%)
    assets/train_sample_ids.txt    — one sample ID per line (training 80%)

Usage:
    python data/setup/create_holdout_split.py \
        [--gt assets/kpsc_amrfinder_ground_truth.tsv] \
        [--store data/inputs/KpSC-HS11286-1000bp-core-npy] \
        [--seed 42] \
        [--holdout-frac 0.20]
"""

import argparse
import os
import numpy as np
import pandas as pd


STRAT_GENES = ["blaKPC", "blaCTX-M", "blaNDM"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gt",           default="assets/kpsc_amrfinder_ground_truth.tsv")
    parser.add_argument("--store",        default="data/inputs/KpSC-HS11286-1000bp-core-npy")
    parser.add_argument("--seed",         type=int,   default=42)
    parser.add_argument("--holdout-frac", type=float, default=0.20)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    # Load all sample IDs present in the numpy store (authoritative order)
    sample_ids = np.load(os.path.join(args.store, "sample_ids.npy"), allow_pickle=True)
    print(f"Total samples in store: {len(sample_ids)}")

    # Load GT for stratification genes (presence/absence only)
    cols = ["sample_id"] + STRAT_GENES
    available_cols = pd.read_csv(args.gt, sep="\t", nrows=0).columns.tolist()
    use_cols = ["sample_id"] + [g for g in STRAT_GENES if g in available_cols]
    gt = pd.read_csv(args.gt, sep="\t", usecols=use_cols).set_index("sample_id")

    # Build stratification key: binary string, e.g. "010" = no KPC, has CTX-M, no NDM
    strat = pd.DataFrame(index=sample_ids)
    for gene in STRAT_GENES:
        if gene in gt.columns:
            strat[gene] = (gt.reindex(strat.index)[gene].fillna(0) >= 1).astype(int)
        else:
            strat[gene] = 0
    strat["stratum"] = strat[STRAT_GENES].astype(str).agg("".join, axis=1)

    holdout_ids = []
    train_ids   = []
    for stratum, grp in strat.groupby("stratum"):
        ids   = grp.index.tolist()
        n_out = max(1, round(len(ids) * args.holdout_frac))
        chosen = rng.choice(ids, size=n_out, replace=False).tolist()
        rest   = [s for s in ids if s not in set(chosen)]
        holdout_ids.extend(chosen)
        train_ids.extend(rest)
        print(f"  stratum={stratum}  n={len(ids)}  held_out={n_out}")

    print(f"\nHold-out: {len(holdout_ids)} samples  |  Train: {len(train_ids)} samples")

    os.makedirs("assets", exist_ok=True)
    for fname, ids in [("holdout_sample_ids.txt", holdout_ids),
                       ("train_sample_ids.txt",   train_ids)]:
        path = os.path.join("assets", fname)
        with open(path, "w") as f:
            f.write("\n".join(ids) + "\n")
        print(f"Saved {path}")


if __name__ == "__main__":
    main()
