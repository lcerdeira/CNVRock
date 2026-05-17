"""Helper: resolve the results root path with a demo fallback.

Production layout: full per-experiment outputs in `data/results/`.

When running on Streamlit Community Cloud or any fresh clone where
`data/results/` does not exist (it's gitignored — the full
reconstructions.npy is ~100 MB per experiment), we fall back to the small
demo bundle committed at `diagnostics/demo/`. That bundle contains the
first 200 samples of exp 32 so reviewers see a working UI immediately.
"""
import os

DATA_RESULTS = "../data/results"
DEMO_RESULTS = "demo"


def resolve_results_root() -> str:
    """First directory among (DATA_RESULTS, DEMO_RESULTS) that exists AND
    contains at least one experiment subdirectory."""
    for candidate in (DATA_RESULTS, DEMO_RESULTS):
        if not os.path.isdir(candidate):
            continue
        subdirs = [d for d in os.listdir(candidate)
                   if os.path.isdir(os.path.join(candidate, d))
                   and not d.startswith("__")]
        if subdirs:
            return candidate
    return DATA_RESULTS                       # error path, keeps message clean


def is_demo() -> bool:
    """True if we are currently serving from the demo bundle."""
    return resolve_results_root() == DEMO_RESULTS
