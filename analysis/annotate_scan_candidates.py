#!/usr/bin/env python3
"""
Enrich the genome-wide CNV scan candidate loci (§3.8) with metabolic-pathway
and functional annotation.

Annotation levels:
  1. Gene product + locus tag        — from HS11286.gff3 (already available)
  2. KEGG Orthology (KO) + pathway   — from the KEGG REST API; HS11286 is
     KEGG organism "kpm", so its genes are pre-annotated to KO / pathway /
     module / EC. No local database or install required.
  3. (separate script) mobile-element context — IntegronFinder / ISEScan.

Outputs
  data/results/cnv_scan_phase_e/candidate_annotation_per_gene.tsv
      one row per candidate gene: bin, locus, product, KO, KO_name,
      pathways, EC
  data/results/cnv_scan_phase_e/candidate_annotation_per_locus.tsv
      one row per candidate locus (contiguous bin cluster): genes,
      dominant pathway(s), n_STs, best_q — the enriched Table 6.

Run locally (needs internet for KEGG REST):
  python3 analysis/annotate_scan_candidates.py
"""
from __future__ import annotations
import io, re, time, json
from pathlib import Path
from urllib.request import urlopen
import pandas as pd

REPO   = Path(__file__).resolve().parents[1]
SCAN   = REPO / "data/results/cnv_scan_phase_e/scan_significant.tsv"
GFF    = REPO / "assets/HS11286.gff3"
OUTDIR = REPO / "data/results/cnv_scan_phase_e"
CACHE  = OUTDIR / "kegg_cache"
ORG    = "kpm"   # KEGG organism code for K. pneumoniae HS11286


def kegg_bulk(endpoint: str) -> str:
    """Fetch a KEGG REST endpoint, cached to disk."""
    CACHE.mkdir(parents=True, exist_ok=True)
    key = endpoint.replace("/", "_") + ".txt"
    fp = CACHE / key
    if fp.exists():
        return fp.read_text()
    url = f"https://rest.kegg.jp/{endpoint}"
    txt = urlopen(url, timeout=60).read().decode()
    fp.write_text(txt)
    time.sleep(0.4)
    return txt


def parse_link(txt: str) -> dict[str, list[str]]:
    """Parse a KEGG /link/ TSV into {gene: [targets]}."""
    d: dict[str, list[str]] = {}
    for line in txt.strip().split("\n"):
        if not line:
            continue
        a, b = line.split("\t")
        d.setdefault(a, []).append(b)
    return d


def parse_list(txt: str) -> dict[str, str]:
    """Parse a KEGG /list/ TSV into {id: description}."""
    d = {}
    for line in txt.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            d[parts[0]] = parts[1]
    return d


def gff_products(gff: Path) -> dict[str, dict]:
    """locus_tag -> {product, start, end, strand} from HS11286 GFF (CDS rows)."""
    out: dict[str, dict] = {}
    for line in gff.read_text().split("\n"):
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(re.findall(r'(\w+)=([^;]+)', f[8]))
        lt = attrs.get("locus_tag")
        if lt:
            out[lt] = {"product": attrs.get("product", "—"),
                       "start": int(f[3]), "end": int(f[4]), "strand": f[6]}
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    scan = pd.read_csv(SCAN, sep="\t")

    # Unique candidate genes + the bins / STs / q they belong to
    gene_rows = []
    for _, r in scan.iterrows():
        genes = str(r.get("gene", "")).split(";")
        for g in genes:
            g = g.strip()
            if g and g != "nan":
                gene_rows.append({"gene": g, "bin": r["bin"],
                                  "antibiotic": r["antibiotic"],
                                  "n_sts": r["n_sts_amplified"], "q": r["q"]})
    gdf = pd.DataFrame(gene_rows)
    genes = sorted(gdf["gene"].unique())
    print(f"Candidate genes: {len(genes)} across {gdf['bin'].nunique()} bins")

    # ── GFF products ──────────────────────────────────────────────────────
    prod = gff_products(GFF)

    # ── KEGG bulk annotation (organism-wide, then subset) ─────────────────
    print("Fetching KEGG annotation for organism", ORG, "…")
    g2ko   = parse_link(kegg_bulk(f"link/ko/{ORG}"))          # gene -> [ko]
    g2path = parse_link(kegg_bulk(f"link/pathway/{ORG}"))     # gene -> [path]
    g2ec   = parse_link(kegg_bulk(f"link/enzyme/{ORG}"))      # gene -> [ec]
    g2mod  = parse_link(kegg_bulk(f"link/module/{ORG}"))      # gene -> [module]
    path_names = parse_list(kegg_bulk(f"list/pathway/{ORG}")) # path -> name
    ko_names   = parse_list(kegg_bulk("list/ko"))             # ko -> name

    def kegg_key(g):  # KEGG keys are "kpm:KPHS_xxxxx"
        return f"{ORG}:{g}"

    # ── Per-gene annotation table ─────────────────────────────────────────
    rows = []
    for g in genes:
        k = kegg_key(g)
        kos  = g2ko.get(k, [])
        ko_id = kos[0].split(":")[-1] if kos else "—"
        ko_nm = ko_names.get(kos[0], "—") if kos else "—"
        paths = [path_names.get(p.replace("path:", ""), p).replace(
                     " - Klebsiella pneumoniae subsp. pneumoniae HS11286", "")
                 for p in g2path.get(k, [])]
        # drop the generic "Metabolic pathways" / "Biosynthesis" umbrellas
        paths_specific = [p for p in paths
                          if p not in ("Metabolic pathways",
                                       "Biosynthesis of secondary metabolites")]
        ecs  = [e.split(":")[-1] for e in g2ec.get(k, [])]
        mods = [m.split(":")[-1] for m in g2mod.get(k, [])]
        meta = gdf[gdf["gene"] == g]
        rows.append({
            "gene":      g,
            "product":   prod.get(g, {}).get("product", "—"),
            "KO":        ko_id,
            "KO_name":   ko_nm,
            "KEGG_pathways":  " | ".join(paths_specific) or "—",
            "KEGG_modules":   " | ".join(mods) or "—",
            "EC":        " | ".join(ecs) or "—",
            "bins":      ",".join(map(str, sorted(meta["bin"].unique()))),
            "max_n_STs": int(meta["n_sts"].max()),
            "best_q":    float(meta["q"].min()),
        })
    per_gene = pd.DataFrame(rows).sort_values(["max_n_STs", "best_q"],
                                              ascending=[False, True])
    per_gene.to_csv(OUTDIR / "candidate_annotation_per_gene.tsv",
                    sep="\t", index=False)
    print(f"wrote candidate_annotation_per_gene.tsv ({len(per_gene)} genes)")

    # ── Per-locus summary (cluster contiguous bins) ───────────────────────
    bins = sorted(gdf["bin"].unique())
    clusters, cur = [], [bins[0]]
    for b in bins[1:]:
        if b - cur[-1] <= 2:          # contiguous (≤2 bins apart)
            cur.append(b)
        else:
            clusters.append(cur); cur = [b]
    clusters.append(cur)

    loc_rows = []
    for cl in clusters:
        sub = gdf[gdf["bin"].isin(cl)]
        cl_genes = per_gene[per_gene["bins"].apply(
            lambda s: any(int(x) in cl for x in s.split(",")))]
        # dominant specific pathways across the cluster's genes
        paths = []
        for p in cl_genes["KEGG_pathways"]:
            paths += [x.strip() for x in p.split("|") if x.strip() != "—"]
        from collections import Counter
        top_paths = "; ".join(f"{p} ({c})" for p, c in
                              Counter(paths).most_common(3)) or "—"
        loc_rows.append({
            "bin_cluster":   f"{cl[0]}–{cl[-1]}" if len(cl) > 1 else str(cl[0]),
            "n_bins":        len(cl),
            "approx_pos_Mb": f"{cl[0]/1000:.2f}",
            "genes":         ", ".join(cl_genes["gene"].tolist()[:6]),
            "products":      " / ".join(
                                 dict.fromkeys(cl_genes["product"].tolist()[:4])),
            "KEGG_pathways": top_paths,
            "EC":            " | ".join(
                                 sorted({e for ec in cl_genes["EC"]
                                         for e in ec.split("|")
                                         if e.strip() != "—"}))[:120] or "—",
            "max_n_STs":     int(sub["n_sts"].max()),
            "best_q":        float(sub["q"].min()),
            "antibiotics":   ", ".join(sorted(sub["antibiotic"].unique())),
        })
    per_locus = pd.DataFrame(loc_rows).sort_values(
        ["max_n_STs", "best_q"], ascending=[False, True])
    per_locus.to_csv(OUTDIR / "candidate_annotation_per_locus.tsv",
                     sep="\t", index=False)
    print(f"wrote candidate_annotation_per_locus.tsv ({len(per_locus)} loci)")

    # ── Console summary of the headline loci ──────────────────────────────
    print("\n── Top candidate loci (enriched Table 6) ──────────────────────")
    for _, r in per_locus.head(8).iterrows():
        print(f"\n  bin {r['bin_cluster']}  ({r['max_n_STs']} STs, "
              f"q={r['best_q']:.2e})  [{r['antibiotics']}]")
        print(f"    genes:    {r['genes']}")
        print(f"    products: {r['products']}")
        print(f"    pathways: {r['KEGG_pathways']}")
        if r['EC'] != "—":
            print(f"    EC:       {r['EC']}")


if __name__ == "__main__":
    main()
