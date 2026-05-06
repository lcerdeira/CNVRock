"""
Re-aggregate kpsc_amrfinder_ground_truth.tsv from existing per-sample AMRFinder
TSVs, adding columns for Phase C genes (qnrB1, blaOXA-48, aac6-Ib-cr).

The per-sample *.amrfinder.tsv files in data/amrfinder_raw/ are unchanged;
only the GENE_PATTERNS dict is extended.  Idempotent: safe to re-run.

Usage (from repo root):
    python3 data/setup/reaggregate_gt.py
"""

import sys
from pathlib import Path

# Extend path so we can reuse get_amrfinder_gt helpers
sys.path.insert(0, str(Path(__file__).parent))

from get_amrfinder_gt import (
    GENE_PATTERNS,
    OUT_TSV,
    RAW_DIR,
    build_sample_map,
    aggregate_results,
)


def main():
    print("Re-aggregating AMRFinder ground truth with updated GENE_PATTERNS...")
    print(f"Genes: {list(GENE_PATTERNS.keys())}")
    print(f"Source TSVs: {RAW_DIR}")
    print(f"Output: {OUT_TSV}")

    sample_to_run = build_sample_map()
    aggregate_results(sample_to_run, RAW_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
