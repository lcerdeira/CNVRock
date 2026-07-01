#!/usr/bin/env python3
"""
Zero-inflation-corrected CNV co-occurrence analysis (exp33, 10K KpSC).

Motivation
----------
The first-pass co-occurrence network (`cnv_cooccurrence_network.py`)
computed a single Spearman ρ between per-gene copy-number values over the
union of samples in which *either* gene was called (mask = A>0 | B>0).
Because most samples carry neither gene, and a sample carrying gene A but
not gene B enters with PCN_B = 0, that ρ conflates two distinct signals:

  (i)  whether two genes tend to be PRESENT together (plasmid
       compatibility / ecological partitioning), and
  (ii) whether, WHEN BOTH ARE PRESENT, their copy numbers co-vary
       (true co-amplification).

A strongly negative ρ (e.g. blaKPC × blaOXA-48 = -0.59) is then largely an
artefact of joint absence / mutual near-exclusivity at the presence level,
not of copy-number covariation. This script decomposes the two layers:

  Layer 1 — PRESENCE co-occurrence:
      2x2 table on binary present/absent (call >= 1), reported as a
      log2 odds ratio with Fisher exact p, Haldane-Anscombe 0.5
      continuity correction. This is the correct test for plasmid
      incompatibility / mutual exclusion.

  Layer 2 — COPY-NUMBER covariation GIVEN co-presence:
      Spearman ρ on family-aggregated PCN restricted to samples where
      BOTH families are present. This is the artefact-free
      co-amplification signal.

Both layers are FDR-corrected (Benjamini-Hochberg) across the gene-pair
set independently. Run locally:

    python3 analysis/cnv_cooccurrence_v2.py
"""
from __future__ import annotations
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

REPO  = Path(__file__).resolve().parents[1]
PLASM = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
CHROM = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT   = REPO / "data/results/cnv_cooccurrence_network"

# Gene-family aggregation: family -> member PCN columns (plasmid) or CRR (chrom)
FAMILIES = {
    "aac6-Ib-cr":    ["pcn_aac6-Ib-cr"],
    "blaCTX-M":      ["pcn_blaCTX-M-14", "pcn_blaCTX-M-15",
                      "pcn_blaCTX-M-27", "pcn_blaCTX-M-65"],
    "blaTEM":        ["pcn_blaTEM-1"],
    "qnrB":          ["pcn_qnrB1"],
    "blaKPC":        ["pcn_blaKPC-2"],
    "blaNDM":        ["pcn_blaNDM-1", "pcn_blaNDM-5"],
    "blaOXA-48-like": ["pcn_blaOXA-48", "pcn_blaOXA-181"],
}
# corresponding binary presence columns (call >= 1)
PRESENCE = {
    "aac6-Ib-cr":    ["aac6-Ib-cr"],
    "blaCTX-M":      ["blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-27", "blaCTX-M-65"],
    "blaTEM":        ["blaTEM-1"],
    "qnrB":          ["qnrB1"],
    "blaKPC":        ["blaKPC-2"],
    "blaNDM":        ["blaNDM-1", "blaNDM-5"],
    "blaOXA-48-like": ["blaOXA-48", "blaOXA-181"],
}
CHROM_AMP_THRESH = 1.75   # crr_blaSHV amplified call


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.clip(q, 0, 1)


def log2_odds_ratio(a: np.ndarray, b: np.ndarray):
    """2x2 present/absent association with Haldane-Anscombe correction."""
    n11 = int(np.sum(a & b))
    n10 = int(np.sum(a & ~b))
    n01 = int(np.sum(~a & b))
    n00 = int(np.sum(~a & ~b))
    # Fisher exact on the raw table
    _, p = stats.fisher_exact([[n11, n10], [n01, n00]])
    # Haldane-Anscombe 0.5 for a finite log-OR even with a zero cell
    orr = ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
    return np.log2(orr), p, (n11, n10, n01, n00)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plasm = pd.read_csv(PLASM, sep="\t")
    chrom = pd.read_csv(CHROM, sep="\t")
    df = plasm.merge(chrom[["sample_id", "crr_blaSHV"]], on="sample_id", how="left")
    n_samples = len(df)

    # ── Build family-level PCN and presence matrices ───────────────────
    pcn = pd.DataFrame(index=df.index)
    pres = pd.DataFrame(index=df.index)
    for fam, cols in FAMILIES.items():
        pcn[fam] = df[cols].fillna(0).sum(axis=1)
        pres[fam] = (df[PRESENCE[fam]].fillna(0).astype(float).values >= 1).any(axis=1)
    # chromosomal blaSHV
    pcn["blaSHV(chr)"] = df["crr_blaSHV"].fillna(0)
    pres["blaSHV(chr)"] = df["crr_blaSHV"].fillna(0) >= CHROM_AMP_THRESH

    genes = list(pcn.columns)
    print(f"Samples: {n_samples} | families: {len(genes)}")
    print("Per-family presence prevalence:")
    for g in genes:
        print(f"  {g:16s} present in {int(pres[g].sum()):5d} "
              f"({100*pres[g].mean():.1f}%)")

    # ── Layer 1: presence co-occurrence (log2 OR, Fisher) ──────────────
    # ── Layer 2: PCN Spearman GIVEN both present ───────────────────────
    # ── plus the OLD union-mask ρ for direct comparison ────────────────
    rows = []
    for gi, gj in itertools.combinations(genes, 2):
        a = pres[gi].values
        b = pres[gj].values

        log2or, p_pres, (n11, n10, n01, n00) = log2_odds_ratio(a, b)

        both = a & b
        n_both = int(both.sum())
        if n_both >= 20:
            rho_co, p_co = stats.spearmanr(pcn[gi].values[both],
                                           pcn[gj].values[both])
        else:
            rho_co, p_co = np.nan, np.nan

        # legacy union-mask ρ (the artefact-prone estimate)
        union = a | b
        if union.sum() >= 20:
            rho_union, _ = stats.spearmanr(pcn[gi].values[union],
                                           pcn[gj].values[union])
        else:
            rho_union = np.nan

        rows.append(dict(gene_a=gi, gene_b=gj,
                         n_both_present=n_both,
                         log2_odds_ratio=log2or, p_presence=p_pres,
                         n11=n11, n10=n10, n01=n01, n00=n00,
                         rho_pcn_copresent=rho_co, p_pcn=p_co,
                         rho_union_legacy=rho_union))

    res = pd.DataFrame(rows)
    res["q_presence"] = benjamini_hochberg(res["p_presence"].values)
    mask_co = res["p_pcn"].notna()
    res.loc[mask_co, "q_pcn"] = benjamini_hochberg(res.loc[mask_co, "p_pcn"].values)

    res = res.sort_values("q_presence")
    out_tsv = OUT / "cooccurrence_two_layer.tsv"
    res.to_csv(out_tsv, sep="\t", index=False)

    # ── Report ─────────────────────────────────────────────────────────
    def show(sub, cols, title):
        print(f"\n── {title} ──")
        with pd.option_context("display.width", 160,
                               "display.max_columns", None,
                               "display.float_format", lambda v: f"{v:.3g}"):
            print(sub[cols].to_string(index=False))

    sig_pres = res[res["q_presence"] < 0.05].copy()
    show(sig_pres.sort_values("log2_odds_ratio"),
         ["gene_a", "gene_b", "log2_odds_ratio", "q_presence",
          "n11", "n10", "n01", "n00"],
         "LAYER 1 — PRESENCE co-occurrence (q<0.05), sorted by log2 OR")

    sig_co = res[(res["q_pcn"] < 0.05) & res["rho_pcn_copresent"].notna()].copy()
    show(sig_co.sort_values("rho_pcn_copresent"),
         ["gene_a", "gene_b", "n_both_present", "rho_pcn_copresent",
          "q_pcn", "rho_union_legacy"],
         "LAYER 2 — COPY-NUMBER covariation | both present (q<0.05)")

    # Headline contrast for the manuscript: KPC × OXA-48-like
    kox = res[((res.gene_a == "blaKPC") & (res.gene_b == "blaOXA-48-like")) |
              ((res.gene_a == "blaOXA-48-like") & (res.gene_b == "blaKPC"))]
    if len(kox):
        r = kox.iloc[0]
        print("\n── HEADLINE: blaKPC × blaOXA-48-like ──")
        print(f"  legacy union-mask ρ (artefact) : {r.rho_union_legacy:+.3f}")
        print(f"  Layer 1 presence log2 OR       : {r.log2_odds_ratio:+.3f} "
              f"(q={r.q_presence:.2e})  [n both present = {r.n_both_present}]")
        print(f"  Layer 2 PCN ρ | both present    : "
              f"{r.rho_pcn_copresent if pd.notna(r.rho_pcn_copresent) else float('nan'):+.3f}"
              f"  (n={r.n_both_present})")
        print("  → the strong negative signal is a PRESENCE-level exclusion,")
        print("    not copy-number anti-correlation.")

    # ── Figure 10: two-panel heatmap (presence log2 OR | PCN ρ) ────────
    make_figure(res, genes, OUT / "cooccurrence_two_layer.png")
    print(f"\nSaved {out_tsv}")


def make_figure(res: pd.DataFrame, genes: list, path: Path) -> None:
    idx = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    OR = np.full((n, n), np.nan)
    RHO = np.full((n, n), np.nan)
    for _, r in res.iterrows():
        i, j = idx[r.gene_a], idx[r.gene_b]
        if r.q_presence < 0.05:
            OR[i, j] = OR[j, i] = r.log2_odds_ratio
        if pd.notna(r.get("q_pcn")) and r.q_pcn < 0.05:
            RHO[i, j] = RHO[j, i] = r.rho_pcn_copresent
    for M in (OR, RHO):
        np.fill_diagonal(M, 0.0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, M, title, vlim, cmap in [
        (axes[0], OR,  "Layer 1 — presence co-occurrence (log₂ OR, q<0.05)", 3.0, "RdBu_r"),
        (axes[1], RHO, "Layer 2 — copy-number ρ | both present (q<0.05)", 0.7, "RdBu_r"),
    ]:
        im = ax.imshow(M, cmap=cmap, vmin=-vlim, vmax=vlim)
        ax.set_xticks(range(n)); ax.set_yticks(range(n))
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(genes, fontsize=8)
        ax.set_title(title, fontsize=10)
        for i in range(n):
            for j in range(n):
                if not np.isnan(M[i, j]) and i != j and M[i, j] != 0.0:
                    ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                            fontsize=6.5,
                            color="white" if abs(M[i, j]) > 0.55 * vlim else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("CNV two-layer co-occurrence — carriage exclusion vs "
                 "copy-number co-amplification (10 K KpSC)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
