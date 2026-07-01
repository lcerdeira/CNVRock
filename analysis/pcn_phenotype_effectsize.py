#!/usr/bin/env python3
"""
Effect sizes with bootstrap CIs for the PCN(R) vs PCN(S) contrasts (Table 5).

The genotype-phenotype summary (`pcn_phenotype_summary.tsv`) reports a
Mann-Whitney U p-value per gene-antibiotic pair. With per-pair n in the
hundreds, p-values are uninformative about magnitude. This script adds the
standard non-parametric effect size for Mann-Whitney — the rank-biserial
correlation, equivalent to the common-language effect size (probability of
superiority) rescaled to [-1, 1]:

    A       = P(PCN_R > PCN_S) + 0.5 * P(PCN_R == PCN_S)     (= AUC)
    r_rb    = 2*A - 1

with a 2 000-resample non-parametric bootstrap 95% CI (resampling the R and
S groups independently). r_rb = 0 is no effect; r_rb -> 1 means resistant
isolates almost always carry higher PCN. Writes
`pcn_phenotype_summary_effectsize.tsv`.

    python3 analysis/pcn_phenotype_effectsize.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
JOINED = REPO / "data/results/pcn_phenotype_joined.tsv"
SUMM   = REPO / "data/results/pcn_phenotype_summary.tsv"
OUT    = REPO / "data/results/pcn_phenotype_summary_effectsize.tsv"

N_BOOT = 2000
MIN_PER_GROUP = 5


def prob_superiority(r: np.ndarray, s: np.ndarray) -> float:
    """AUC = P(R>S) + 0.5 P(R==S) via the Mann-Whitney U statistic."""
    n1, n2 = len(r), len(s)
    if n1 == 0 or n2 == 0:
        return np.nan
    u, _ = stats.mannwhitneyu(r, s, alternative="two-sided")
    return u / (n1 * n2)


def main() -> None:
    j = pd.read_csv(JOINED, sep="\t",
                    usecols=["gene", "antibiotic_name", "pcn", "updated_phenotype"])
    j = j.dropna(subset=["pcn", "updated_phenotype"])
    rng = np.random.default_rng(42)

    rows = []
    for (gene, ab), g in j.groupby(["gene", "antibiotic_name"]):
        ph = g["updated_phenotype"].astype(str).str.lower()
        r = g.loc[ph.str.startswith("r"), "pcn"].to_numpy(float)
        s = g.loc[ph.str.startswith("s"), "pcn"].to_numpy(float)
        if len(r) < MIN_PER_GROUP or len(s) < MIN_PER_GROUP:
            continue
        A = prob_superiority(r, s)
        r_rb = 2 * A - 1
        boot = np.empty(N_BOOT)
        for b in range(N_BOOT):
            rb = rng.choice(r, len(r), replace=True)
            sb = rng.choice(s, len(s), replace=True)
            boot[b] = 2 * prob_superiority(rb, sb) - 1
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append(dict(gene=gene, antibiotic=ab, n_R=len(r), n_S=len(s),
                         prob_superiority=round(A, 3),
                         rank_biserial_r=round(r_rb, 3),
                         rb_ci_lo=round(lo, 3), rb_ci_hi=round(hi, 3),
                         ci_excludes_zero=bool(lo > 0 or hi < 0)))

    res = pd.DataFrame(rows).sort_values("rank_biserial_r", ascending=False)
    res.to_csv(OUT, sep="\t", index=False)
    with pd.option_context("display.width", 160, "display.max_rows", None):
        print(res.to_string(index=False))
    n_sig = int(res["ci_excludes_zero"].sum())
    print(f"\n{n_sig}/{len(res)} gene-antibiotic pairs have a 95% CI that "
          f"excludes zero effect.")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
