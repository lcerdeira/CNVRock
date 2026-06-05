"""Helper: resolve the results root path with a demo fallback.

Production layout: full per-experiment outputs in `data/results/`.

When running on Streamlit Community Cloud or any fresh clone where
`data/results/` does not exist (it's gitignored — the full
reconstructions.npy is ~100 MB per experiment), we fall back to the small
demo bundle committed at `diagnostics/demo/`. That bundle contains the
first 200 samples of exp 32 so reviewers see a working UI immediately.

All paths are resolved relative to THIS FILE's location so the app works
regardless of where streamlit is launched from (Streamlit Cloud runs from
the repo root; running locally with `cd diagnostics && streamlit run app.py`
runs from `diagnostics/`).
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))               # …/diagnostics/src
_DIAGNOSTICS = os.path.dirname(_HERE)                            # …/diagnostics
_REPO_ROOT = os.path.dirname(_DIAGNOSTICS)                       # …/CNVRock

DATA_RESULTS = os.path.join(_REPO_ROOT, "data", "results")       # production
DEMO_RESULTS = os.path.join(_DIAGNOSTICS, "demo")                # bundled


_EXPERIMENT_MARKERS = ("latents.npy", "reconstructions.npy", "training_log.json")


def _looks_like_experiment(path: str) -> bool:
    """A directory is a real experiment output only if it carries at least
    one model artefact. This guards against stray committed sub-directories
    (e.g. `data/results/cnv_scan_phase_e/` holding only annotation TSVs),
    which must NOT flip the app out of demo mode on a fresh clone / cloud."""
    return any(os.path.isfile(os.path.join(path, m)) for m in _EXPERIMENT_MARKERS)


def resolve_results_root() -> str:
    """First directory among (DATA_RESULTS, DEMO_RESULTS) that exists AND
    contains at least one real experiment subdirectory (model artefacts —
    not just any committed folder)."""
    for candidate in (DATA_RESULTS, DEMO_RESULTS):
        if not os.path.isdir(candidate):
            continue
        subdirs = [d for d in os.listdir(candidate)
                   if os.path.isdir(os.path.join(candidate, d))
                   and not d.startswith("__")
                   and _looks_like_experiment(os.path.join(candidate, d))]
        if subdirs:
            return candidate
    return DATA_RESULTS                       # error path, keeps message clean


def is_demo() -> bool:
    """True if we are currently serving from the demo bundle."""
    return resolve_results_root() == DEMO_RESULTS


def list_demo_bundles() -> dict:
    """Return {bundle_name: abs_path} for every demo bundle that has the
    minimum files for the sample viewer (latents + reconstructions)."""
    out = {}
    if not os.path.isdir(DEMO_RESULTS):
        return out
    for d in sorted(os.listdir(DEMO_RESULTS)):
        p = os.path.join(DEMO_RESULTS, d)
        if (os.path.isdir(p)
                and os.path.isfile(os.path.join(p, "latents.npy"))
                and os.path.isfile(os.path.join(p, "reconstructions.npy"))):
            out[d] = p
    return out
