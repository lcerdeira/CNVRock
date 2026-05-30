#!/usr/bin/env python3
"""
Build the Phase E extended ground-truth file.

The AMRFinder GT (amrfinder_gt_expansion.tsv) only has columns for the
original 7 gene families. The 6 Phase E genes need ground truth too. We
derive it from the Kleborate v3 genotype table, which already reports the
relevant resistance classes per sample:

  Phase E gene   Kleborate column      match string     rationale
  -----------    -----------------     ------------     ---------
  sul1           Sul_acquired          "sul1"           clean, variant-specific
  sul2           Sul_acquired          "sul2"           clean, variant-specific
  dfrA12         Tmt_acquired          "dfrA12"         clean, variant-specific
  dfrA14         Tmt_acquired          "dfrA14"         clean, variant-specific
  aac3-IIa       AGly_acquired         "aac(3)-II"      FAMILY-level — the
                  reference contig cross-maps across the aac(3)-II subgroup,
                  so a variant-specific call would be unfair; matches the
                  gene-family-aggregation principle (Methods 2.4).
  blaOXA-232 is folded into the existing blaOXA-48-like family and is
  already covered by the AMRFinder blaOXA-48-like column; no new column.

Output: assets/amrfinder_gt_expansion_phaseE.tsv  (AMRFinder GT + 5 columns)
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
AMR_GT = REPO / "assets/amrfinder_gt_expansion.tsv"
KB_GT = REPO / "assets/kpsc_expansion_kleborate_gt_runlevel.tsv"
OUT = REPO / "assets/amrfinder_gt_expansion_phaseE.tsv"

# new GT column -> (kleborate column, match substring)
NEW_GT = {
    "sul1":     ("Sul_acquired", "sul1"),
    "sul2":     ("Sul_acquired", "sul2"),
    "dfrA12":   ("Tmt_acquired", "dfrA12"),
    "dfrA14":   ("Tmt_acquired", "dfrA14"),
    # aac(3)-II removed (2026-05-30): reference contig was aac(3)-IIa (zero
    # prevalence in KpSC); cross-mapping from IId/IIe made family-level GT
    # necessary but scientifically inaccurate. Removed from panel.
}


def main():
    amr = pd.read_csv(AMR_GT, sep="\t", dtype=str)
    kb = pd.read_csv(KB_GT, sep="\t", dtype=str)
    print(f"AMRFinder GT: {len(amr):,} rows, {len(amr.columns)} cols")
    print(f"Kleborate GT: {len(kb):,} rows")

    id_col = "sample_id"
    kb_idx = kb.set_index(id_col)

    for gene, (kb_col, needle) in NEW_GT.items():
        if kb_col not in kb_idx.columns:
            print(f"  WARNING: {kb_col} missing — skipping {gene}")
            continue
        present = kb_idx[kb_col].fillna("-").astype(str).str.contains(
            needle, case=False, regex=False)
        gt_series = present.astype(int)
        amr[gene] = amr[id_col].map(gt_series).fillna(0).astype(int)
        n_pos = int(amr[gene].sum())
        print(f"  {gene:10s} <- {kb_col} ~ '{needle}'   positives: {n_pos:,}"
              f" / {len(amr):,}")

    amr.to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT}  ({len(amr.columns)} columns)")


if __name__ == "__main__":
    main()
