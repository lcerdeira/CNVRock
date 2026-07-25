#!/usr/bin/env python3
"""
CNV co-occurrence network for plasmid resistance genes (exp33, 10K KpSC).

Builds a gene × gene undirected weighted graph where:
  - Nodes = resistance genes (plasmid PCN + chromosomal blaSHV CRR)
  - Edge weight = Spearman ρ between per-gene copy-number values
    (FDR-corrected, Benjamini-Hochberg q < 0.05)

Graph-theoretic analyses:
  1. Betweenness centrality — which gene is the "hub" of resistance networks?
  2. Community detection (Louvain) — which genes co-amplify?
  3. Positive vs negative edges — amplification synergy vs exclusion
     (negative ρ: genes on incompatible replicons?)

Run locally:
    python3 analysis/cnv_cooccurrence_network.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from scipy import stats
from pathlib import Path

REPO    = Path(__file__).resolve().parents[1]
PLASM   = REPO / "data/results/33_kpsc_expansion_10k/plasmid_gene_calls.tsv"
CHROM   = REPO / "data/results/33_kpsc_expansion_10k/gene_calls.tsv"
OUT     = REPO / "data/results/cnv_cooccurrence_network"


def benjamini_hochberg(pvals: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    qvals = pvals * n / ranks
    # enforce monotonicity from the right
    for i in range(n - 2, -1, -1):
        qvals[order[i]] = min(qvals[order[i]], qvals[order[i + 1]])
    return np.clip(qvals, 0, 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────
    plasm = pd.read_csv(PLASM, sep="\t")
    chrom = pd.read_csv(CHROM, sep="\t")

    # PCN columns (per-gene family-aggregated plasmid copy number)
    pcn_cols = [c for c in plasm.columns
                if c.startswith("pcn_") and not c.endswith("_fam")]
    # Rename for display
    pcn_labels = {c: c.replace("pcn_", "") for c in pcn_cols}

    # Add chromosomal blaSHV CRR (if available)
    if "crr_blaSHV" in chrom.columns:
        pcn_cols += ["crr_blaSHV"]
        plasm = plasm.merge(chrom[["sample_id", "crr_blaSHV"]],
                            on="sample_id", how="left")

    print(f"Genes in network: {len(pcn_cols)}")
    print(f"Samples: {len(plasm)}")

    X = plasm[pcn_cols].fillna(0).values.astype(float)

    # ── Spearman correlation matrix ────────────────────────────────────
    n = len(pcn_cols)
    rho_mat  = np.zeros((n, n))
    pval_mat = np.ones((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            # Only use samples where both genes are called (PCN > 0 threshold)
            mask = (X[:, i] > 0) | (X[:, j] > 0)
            if mask.sum() < 20:
                continue
            rho, p = stats.spearmanr(X[mask, i], X[mask, j])
            rho_mat[i, j] = rho_mat[j, i] = rho
            pval_mat[i, j] = pval_mat[j, i] = p

    np.fill_diagonal(rho_mat, 1.0)

    # FDR correction on upper triangle
    upper = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pvals_flat = np.array([pval_mat[i, j] for i, j in upper])
    qvals_flat = benjamini_hochberg(pvals_flat)
    qval_mat = np.ones((n, n))
    for (i, j), q in zip(upper, qvals_flat):
        qval_mat[i, j] = qval_mat[j, i] = q

    sig_mask = qval_mat < 0.05
    print(f"\nSignificant pairs (q < 0.05): {sig_mask.sum() // 2} / {len(upper)}")

    # ── Build NetworkX graph ───────────────────────────────────────────
    labels = [pcn_labels.get(c, c.replace("pcn_", "")) for c in pcn_cols]

    G = nx.Graph()
    for i, lab in enumerate(labels):
        G.add_node(i, label=lab)

    for i in range(n):
        for j in range(i + 1, n):
            if sig_mask[i, j]:
                G.add_edge(i, j,
                           weight=float(abs(rho_mat[i, j])),
                           rho=float(rho_mat[i, j]),
                           q=float(qval_mat[i, j]))

    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ── Graph-theoretic analysis ───────────────────────────────────────
    # Betweenness centrality
    bc = nx.betweenness_centrality(G, weight="weight")
    print("\n── Betweenness centrality (top 5) ───────────────────────────")
    for node, score in sorted(bc.items(), key=lambda x: -x[1])[:5]:
        print(f"  {labels[node]:20s}: {score:.4f}")

    # Degree centrality
    dc = nx.degree_centrality(G)
    print("\n── Degree centrality (top 5) ─────────────────────────────────")
    for node, score in sorted(dc.items(), key=lambda x: -x[1])[:5]:
        print(f"  {labels[node]:20s}: {score:.4f}")

    # Community detection (Louvain if available, else Girvan-Newman)
    try:
        from networkx.algorithms import community as nx_comm
        comms = list(nx_comm.louvain_communities(G, seed=42))
        print(f"\n── Louvain communities: {len(comms)} ─────────────────────────")
        for i, comm in enumerate(comms):
            mems = ", ".join(labels[n] for n in sorted(comm))
            print(f"  Community {i+1}: {mems}")
        node_comm = {}
        for ci, comm in enumerate(comms):
            for node in comm:
                node_comm[node] = ci
    except Exception:
        node_comm = {n: 0 for n in G.nodes()}

    # Positive vs negative correlation summary
    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d["rho"] > 0]
    neg_edges = [(u, v) for u, v, d in G.edges(data=True) if d["rho"] < 0]
    print(f"\nPositive correlations (co-amplification): {len(pos_edges)}")
    print(f"Negative correlations (mutual exclusion):  {len(neg_edges)}")
    for u, v, d in sorted(G.edges(data=True), key=lambda x: x[2]["rho"]):
        if d["rho"] < -0.05:
            print(f"  {labels[u]:15s} ↔ {labels[v]:15s}: ρ={d['rho']:+.3f}  q={d['q']:.3e}")

    # ── Save correlation matrix ────────────────────────────────────────
    rho_df = pd.DataFrame(rho_mat, index=labels, columns=labels)
    qval_df = pd.DataFrame(qval_mat, index=labels, columns=labels)
    rho_df.to_csv(OUT / "cnv_spearman_rho.tsv", sep="\t")
    qval_df.to_csv(OUT / "cnv_spearman_qval.tsv", sep="\t")

    bc_df = pd.DataFrame([{"gene": labels[n], "betweenness": v,
                            "degree": dc[n], "community": node_comm.get(n, -1)}
                           for n, v in bc.items()])
    bc_df.sort_values("betweenness", ascending=False, inplace=True)
    bc_df.to_csv(OUT / "cnv_network_centrality.tsv", sep="\t", index=False)

    # ── Figure ─────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: correlation heatmap
    ax = axes[0]
    im = ax.imshow(rho_mat, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45,
                                                  ha="right", fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=7)
    # Overlay significance stars
    for i in range(n):
        for j in range(n):
            if i != j and qval_mat[i, j] < 0.05:
                ax.text(j, i, "·", ha="center", va="center",
                        fontsize=10, color="black")
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title("Spearman ρ (· = q<0.05)", fontsize=9)

    # Right: network graph
    ax = axes[1]
    comm_colors = plt.cm.Set2(np.linspace(0, 1, max(node_comm.values()) + 2))
    node_colors = [comm_colors[node_comm.get(n, 0)] for n in G.nodes()]
    pos = nx.spring_layout(G, weight="weight", seed=42, k=1.5)

    # Edge colors: blue = positive rho, red = negative
    edge_colors = ["#2176ae" if G[u][v]["rho"] > 0 else "#e63946"
                   for u, v in G.edges()]
    edge_widths = [3 * abs(G[u][v]["rho"]) for u, v in G.edges()]
    node_sizes  = [300 + 3000 * bc[n] for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                           width=edge_widths, alpha=0.7)
    nx.draw_networkx_labels(G, pos, ax=ax,
                            labels={n: labels[n] for n in G.nodes()},
                            font_size=7)
    ax.set_title("CNV co-occurrence network\n"
                 "Node size ∝ betweenness | blue=co-amp | red=exclusive",
                 fontsize=9)
    ax.axis("off")

    legend_patches = [
        mpatches.Patch(color="#2176ae", label="Co-amplification (ρ>0)"),
        mpatches.Patch(color="#e63946", label="Mutual exclusion (ρ<0)"),
    ]
    ax.legend(handles=legend_patches, fontsize=7, loc="lower right")

    plt.tight_layout()
    out_fig = OUT / "cnv_cooccurrence_network.png"
    plt.savefig(out_fig, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved: {out_fig}")
    print(f"Network TSVs saved to {OUT}/")


if __name__ == "__main__":
    main()
