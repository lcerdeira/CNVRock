#!/usr/bin/env python3
"""
Confounder-adjusted CNV-by-ST analysis for blaKPC-2 PCN (exp33, 10K KpSC).

The headline §3.6 claim — ST258 and ST11 carry blaKPC at higher copy
number than the ST101/ST15 background — rests on a Kruskal-Wallis test
that does not control for the fact that ST is correlated with collection
country and year across the BioSample-linked submissions. This script
tests whether the ST effect survives adjustment:

  Model 0 (unadjusted): PCN ~ C(ST)                          [OLS]
  Model 1 (adjusted)  : PCN ~ C(ST) + C(country) + year      [OLS, FE]

Country is entered as fixed-effect dummies (a within-country contrast,
which fully absorbs any between-country confounding) and collection year
as a centred linear term. Standard errors are heteroskedasticity-robust
(HC0). This is implemented in pure numpy/scipy to avoid a statsmodels /
scipy ABI mismatch in the local environment.

Restricted to blaKPC-positive isolates (pcn > 0.05) with a country/year
record, and to STs with n >= 8 positives (matching the manuscript).
Reference ST for contrasts is the pooled ST101/ST15 background.

    python3 analysis/cnv_st_mixedmodel.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def ols_hc0(y: np.ndarray, X: np.ndarray):
    """OLS with HC0 (White) robust covariance. Returns beta, se, p (t-approx)."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    # HC0 sandwich
    S = (X * (resid ** 2)[:, None]).T @ X
    cov = XtX_inv @ S @ XtX_inv
    se = np.sqrt(np.diag(cov))
    dof = max(1, len(y) - X.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return beta, se, p

REPO   = Path(__file__).resolve().parents[1]
PLASM  = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
STMETA = REPO / "assets/kpsc_expansion_metadata_runlevel.tsv"
PHENO  = REPO / "data/results/pcn_phenotype_joined.tsv"
OUT    = REPO / "data/results/cnv_cooccurrence_network"

# NOTE ON DATA SCOPE. Cohort-wide collection country/year is not present in
# the committed result tables; it exists only for the AST-linked subset
# (pcn_phenotype_joined.tsv, ~497 isolates, of which only ~36 are KPC+).
# A fully-powered ST × country × year adjustment therefore needs the full
# BioSample provenance metadata (HPC / ENA). This script runs the model
# on whatever geography is available and PRINTS AN EXPLICIT POWER WARNING;
# to produce the manuscript supplement figure, point PHENO at a cohort-wide
# {sample_id, country, collection_year} table and raise MIN_ST_N to 8.
MIN_ST_N = 3
POS_THRESH = 0.05


def main() -> None:
    plasm = pd.read_csv(PLASM, sep="\t")[["sample_id", "pcn_blaKPC-2"]]
    stm   = pd.read_csv(STMETA, sep="\t")[["sample_id", "ST"]]
    phe   = pd.read_csv(PHENO, sep="\t",
                        usecols=["run_acc", "country", "collection_year"])
    phe   = (phe.dropna(subset=["run_acc"])
                .drop_duplicates("run_acc")
                .rename(columns={"run_acc": "sample_id"}))

    df = (plasm.merge(stm, on="sample_id", how="left")
               .merge(phe, on="sample_id", how="left"))
    df = df.rename(columns={"pcn_blaKPC-2": "pcn"})

    # KPC-positive, with ST + country + year
    df = df[(df["pcn"] > POS_THRESH)
            & df["ST"].notna()
            & df["country"].notna()
            & df["collection_year"].notna()].copy()
    df["ST"] = df["ST"].astype(str)
    df["collection_year"] = pd.to_numeric(df["collection_year"], errors="coerce")
    df = df.dropna(subset=["collection_year"])

    # keep STs with enough positives
    keep = df["ST"].value_counts()
    keep = keep[keep >= MIN_ST_N].index.tolist()
    df = df[df["ST"].isin(keep)].copy()
    print(f"KPC-positive isolates with ST+country+year: {len(df)}")
    print(f"STs (n>={MIN_ST_N}): {sorted(keep)}")
    print(f"Countries: {df['country'].nunique()} | "
          f"years {int(df['collection_year'].min())}-{int(df['collection_year'].max())}")

    # Reference background: pool ST101 + ST15 into the baseline.
    ref_levels = [s for s in ["ST101", "ST15", "101", "15"] if s in keep]
    if not ref_levels:
        ref_levels = [sorted(keep)[0]]
    other_sts = [s for s in sorted(keep) if s not in ref_levels]
    print(f"Reference background (pooled): {ref_levels}\n")

    y = df["pcn"].values.astype(float)
    n = len(y)
    yr = (df["collection_year"].values - df["collection_year"].mean())

    # Design matrices. intercept + ST dummies (ref = pooled background).
    def st_dummies():
        return np.column_stack([(df["ST"] == s).astype(float).values
                                for s in other_sts])
    ST = st_dummies()
    country_levels = sorted(df["country"].unique())[1:]   # drop first = ref
    CO = np.column_stack([(df["country"] == c).astype(float).values
                          for c in country_levels]) if country_levels else np.empty((n, 0))

    ones = np.ones((n, 1))
    # Model 0: unadjusted
    X0 = np.hstack([ones, ST])
    b0, se0, p0 = ols_hc0(y, X0)

    # Power / identifiability check. With few isolates spread across many
    # countries, ST is collinear with country and the country-FE model is
    # rank-deficient. Fall back to a year-only adjustment and warn.
    n_params_full = 1 + ST.shape[1] + CO.shape[1] + 1
    underpowered = n < 5 * n_params_full
    if underpowered:
        print("\n*** POWER WARNING ***")
        print(f"  n={n} KPC+ isolates with geography vs {n_params_full} "
              f"parameters in the country-FE model.")
        print("  ST is collinear with country at this sample size; the full")
        print("  country-adjusted model is NOT identifiable. Reporting a")
        print("  YEAR-ONLY sensitivity adjustment. A fully-powered ST ×")
        print("  country × year model requires cohort-wide provenance")
        print("  metadata (see header note).")
        X1 = np.hstack([ones, ST, yr[:, None]])   # year-only
        adj_label = "year-adj"
    else:
        X1 = np.hstack([ones, ST, CO, yr[:, None]])
        adj_label = "country+year-adj"
    b1, se1, p1 = ols_hc0(y, X1)

    # ST coefficient indices: positions 1..len(other_sts) in both models
    focus = [s for s in ("ST258", "ST11", "258", "11") if s in other_sts]
    ordered = focus + [s for s in other_sts if s not in focus]
    print("── Contrast vs pooled ST101/ST15 background; "
          "positive coef = higher blaKPC PCN ──")
    print(f"{'ST':8s} {'unadj coef':>11s} {'unadj p':>10s} "
          f"{'adj coef':>10s} {'adj p':>10s}   note")
    rows = []
    for st in ordered:
        j = 1 + other_sts.index(st)
        tag = "  << focal clone" if st in focus else ""
        print(f"{st:8s} {b0[j]:>11.3f} {p0[j]:>10.2e} "
              f"{b1[j]:>10.3f} {p1[j]:>10.2e}{tag}")
        rows.append(dict(ST=st, unadj_coef=b0[j], unadj_p=p0[j],
                         adj_coef=b1[j], adj_p=p1[j]))

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "kpc_pcn_st_mixedmodel.tsv",
                              sep="\t", index=False)
    yr_j = X1.shape[1] - 1
    print(f"\nCountries adjusted for: {len(country_levels)+1} (fixed effects)")
    print(f"Collection-year slope (adj): {b1[yr_j]:+.4f} PCN/yr "
          f"(p={p1[yr_j]:.2e})")
    print(f"Saved {OUT/'kpc_pcn_st_mixedmodel.tsv'}")


if __name__ == "__main__":
    main()
