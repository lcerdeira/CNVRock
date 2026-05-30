#!/usr/bin/env python3
"""
Stochastic birth-death model for AMR gene copy-number amplification.

Theory (Andersson & Hughes, Science 2009):
  Under a linear birth-death process where each gene copy can independently:
    - Amplify at rate λ (tandem duplication per replication)
    - De-amplify at rate μ (deletion per replication)

  The steady-state copy-number distribution is geometric:
    P(CN = k) = (1 - p) * p^(k-1),  k ≥ 1,  p = λ/(λ+μ)

  Under NEUTRAL evolution: λ ≈ μ → p ≈ 0.5
  Under POSITIVE SELECTION: λ > μ → p > 0.5 (more amplified isolates)
  Under PURIFYING SELECTION: λ < μ → p < 0.5 (low copy maintained)

  Selection coefficient estimate: s ≈ (λ - μ) / (λ + μ) = 2p - 1

For plasmid genes (continuous PCN), we discretise:
  CN=0: PCN < 0.2  (absent)
  CN=1: 0.2 ≤ PCN < 1.5  (single copy / moderate)
  CN=2: 1.5 ≤ PCN < 2.5
  CN=3: 2.5 ≤ PCN < 3.5
  CN≥4: PCN ≥ 3.5

For chromosomal blaSHV (direct CN calls):
  Use the cn column from gene_calls.tsv.

Run locally:
    python3 analysis/birth_death_cnv.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.stats import geom, kstest
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
PLASM  = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
CHROM  = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT    = REPO / "data/results/birth_death_cnv"


def pcn_to_cn(pcn: pd.Series) -> pd.Series:
    """Discretise continuous PCN to integer copy number."""
    bins = [-np.inf, 0.2, 1.5, 2.5, 3.5, np.inf]
    labels = [0, 1, 2, 3, 4]
    return pd.cut(pcn, bins=bins, labels=labels).astype(float)


def fit_birth_death(cn_counts: dict[int, int]) -> dict:
    """
    Fit geometric distribution to CN counts (restricted to CN ≥ 1, i.e.
    among isolates that carry the gene).

    Returns: p (geometric parameter), λ/μ ratio, selection coeff s,
             log-likelihood, KS p-value.
    """
    # Exclude CN=0 (absent) — birth-death model is for carriers only
    counts = {k: v for k, v in cn_counts.items() if k >= 1}
    if not counts or sum(counts.values()) < 10:
        return {}

    total = sum(counts.values())
    obs = np.array([counts.get(k, 0) for k in range(1, 6)])
    obs_norm = obs / obs.sum()

    # MLE for geometric: p̂ = 1 / mean_CN
    cn_vals = np.repeat(list(range(1, 6)), obs)
    p_hat = 1.0 / cn_vals.mean()

    # Bootstrap CI for p
    rng = np.random.default_rng(42)
    p_boot = []
    for _ in range(500):
        sample = rng.choice(cn_vals, size=len(cn_vals), replace=True)
        p_boot.append(1.0 / sample.mean())
    p_ci = np.percentile(p_boot, [2.5, 97.5])

    # Selection coefficient
    s = 2 * p_hat - 1   # = (λ-μ)/(λ+μ)
    lam_over_mu = (1 - p_hat) / p_hat   # λ/μ

    # KS test: does empirical CN dist match geometric(p_hat)?
    cdf_obs = np.cumsum(obs_norm)
    cdf_exp = geom.cdf(range(1, 6), p=p_hat)
    ks_stat = np.max(np.abs(cdf_obs - cdf_exp))

    # Log-likelihood
    ll = sum(c * np.log(geom.pmf(k, p=p_hat) + 1e-12)
             for k, c in counts.items())

    return {
        "n_carriers":   total,
        "mean_cn":      cn_vals.mean(),
        "p_hat":        p_hat,
        "p_ci_lo":      p_ci[0],
        "p_ci_hi":      p_ci[1],
        "lambda_over_mu": lam_over_mu,
        "selection_s":  s,
        "selection_interp": (
            "strong positive" if s > 0.3 else
            "moderate positive" if s > 0.1 else
            "weak positive" if s > 0 else
            "neutral" if abs(s) < 0.05 else
            "purifying"
        ),
        "ks_stat":      ks_stat,
        "log_lik":      ll,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    plasm = pd.read_csv(PLASM, sep="\t")
    chrom = pd.read_csv(CHROM, sep="\t")

    results = []

    # ── Plasmid genes ─────────────────────────────────────────────────
    pcn_cols = [c for c in plasm.columns
                if c.startswith("pcn_") and not c.endswith("_fam")]
    print("── Plasmid gene birth-death models ──────────────────────────")
    for col in sorted(pcn_cols):
        gene = col.replace("pcn_", "")
        cn = pcn_to_cn(plasm[col].fillna(0))
        cn_counts = cn.value_counts().to_dict()
        res = fit_birth_death({int(k): int(v) for k, v in cn_counts.items()})
        if res:
            res["gene"] = gene
            res["track"] = "plasmid"
            results.append(res)
            print(f"  {gene:18s}: λ/μ={res['lambda_over_mu']:.3f}  "
                  f"s={res['selection_s']:+.3f}  [{res['selection_interp']}]  "
                  f"mean_CN={res['mean_cn']:.2f}  n={res['n_carriers']}")

    # ── Chromosomal blaSHV ────────────────────────────────────────────
    print("\n── Chromosomal blaSHV birth-death model ──────────────────────")
    if "cn" in chrom.columns:
        cn_blashv = chrom["cn"].dropna()
        cn_counts = cn_blashv.value_counts().to_dict()
    elif "crr_blaSHV" in chrom.columns:
        cn_blashv = pcn_to_cn(chrom["crr_blaSHV"].fillna(1.0))
        cn_counts = cn_blashv.value_counts().to_dict()
    else:
        cn_counts = {}

    if cn_counts:
        # Get CN counts from actual gene calls column
        if "blaSHV" in chrom.columns:
            cn_col = chrom["blaSHV"].fillna(-1)
            cn_counts_direct = {int(k): int(v)
                                 for k, v in cn_col.value_counts().items()
                                 if k >= 1}
            res = fit_birth_death(cn_counts_direct)
        else:
            res = fit_birth_death({int(k): int(v)
                                   for k, v in cn_counts.items()})
        if res:
            res["gene"] = "blaSHV (chromosomal)"
            res["track"] = "chromosomal"
            results.append(res)
            print(f"  {'blaSHV (chromosomal)':25s}: "
                  f"λ/μ={res['lambda_over_mu']:.3f}  "
                  f"s={res['selection_s']:+.3f}  [{res['selection_interp']}]  "
                  f"mean_CN={res['mean_cn']:.2f}")

    # ── Null comparison: what does neutral look like? ─────────────────
    print("\n── Reference: neutral evolution (λ=μ, p=0.5) ───────────────")
    print("  Neutral: s=0.00, λ/μ=1.00, mean_CN=2.00")
    print("  Genes with s >> 0 are under positive selection for amplification")

    # ── Save & figure ─────────────────────────────────────────────────
    df = pd.DataFrame(results)
    df.sort_values("selection_s", ascending=False, inplace=True)
    df.to_csv(OUT / "birth_death_results.tsv", sep="\t", index=False)
    print(f"\nSaved {OUT}/birth_death_results.tsv")

    if not df.empty:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # λ/μ bar chart
        ax = axes[0]
        colors = ["#e63946" if s > 0.1 else "#2176ae" if s < -0.05 else "#adb5bd"
                  for s in df["selection_s"]]
        bars = ax.barh(range(len(df)), df["lambda_over_mu"], color=colors,
                       alpha=0.85)
        ax.axvline(1.0, color="black", ls="--", lw=1, label="Neutral (λ/μ=1)")
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df["gene"].str.replace("pcn_",""), fontsize=8)
        ax.set_xlabel("λ/μ  (amplification / de-amplification rate ratio)", fontsize=9)
        ax.set_title("Selective pressure on gene amplification\n"
                     "Red = positive selection | Blue = purifying", fontsize=9)
        ax.legend(fontsize=8)

        # Selection coefficient s with CI
        ax = axes[1]
        x = range(len(df))
        ax.barh(x, df["selection_s"], xerr=[df["selection_s"]-df["p_ci_lo"]*2+1,
                                             df["p_ci_hi"]*2-1-df["selection_s"]],
                color=colors, alpha=0.8, capsize=3, error_kw={"lw":1})
        ax.axvline(0, color="black", ls="-", lw=0.8)
        ax.axvline(0.1, color="grey", ls="--", lw=0.7, label="s=0.1 threshold")
        ax.set_yticks(x)
        ax.set_yticklabels(df["gene"].str.replace("pcn_",""), fontsize=8)
        ax.set_xlabel("Selection coefficient s = (λ-μ)/(λ+μ)", fontsize=9)
        ax.set_title("s > 0: amplification favoured\ns < 0: amplification purged",
                     fontsize=9)
        ax.legend(fontsize=8)

        plt.tight_layout()
        fig.savefig(OUT / "birth_death_selection.png", dpi=150,
                    bbox_inches="tight")
        print(f"Figure saved: {OUT}/birth_death_selection.png")

    print("\n── Interpretation ────────────────────────────────────────────")
    if not df.empty:
        top = df.iloc[0]
        print(f"  Strongest positive selection: {top['gene']} "
              f"(s={top['selection_s']:+.3f}, λ/μ={top['lambda_over_mu']:.2f})")
        print(f"  → Amplification of {top['gene']} is ~{top['lambda_over_mu']:.1f}× "
              f"more likely than de-amplification under current selective pressure.")
        neutral = df[df["selection_s"].abs() < 0.05]
        if len(neutral):
            print(f"  Neutral genes (s≈0): {', '.join(neutral['gene'].tolist())}")


if __name__ == "__main__":
    main()
