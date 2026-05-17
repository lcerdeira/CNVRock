import os
import streamlit as st
import streamlit.components.v1 as components

from src.utils import load_results, load_meta, compute_coverage, plot_coverage
from src.results_root import resolve_results_root


def page2():
    st.title("Latent space coverage")

    RESULTS_ROOT = resolve_results_root()
    if not os.path.isdir(RESULTS_ROOT) or not os.listdir(RESULTS_ROOT):
        st.warning(f"No results directories found in {RESULTS_ROOT}.")
        st.stop()
    RESULTS_DIR = st.selectbox("Select results directory",
                               options=sorted(os.listdir(RESULTS_ROOT)), index=0)
    if not RESULTS_DIR:
        st.stop("Please select a results directory.")

    results    = load_results(os.path.join(RESULTS_ROOT, RESULTS_DIR))
    meta, _gff = load_meta()  # uses default KpSC path or Pf9 fallback

    with st.expander("Settings", expanded=False):
        n_void      = st.slider("Void probes",    1_000, 50_000,  5_000, step=1_000)
        n_neighbors = st.slider("UMAP n_neighbors", 5, 100, 15)
        min_dist    = st.slider("UMAP min_dist",    0.0, 0.5, 0.05, step=0.05)

    with st.spinner("Computing coverage (cached after first run)…"):
        coverage_df = compute_coverage(
            results["latents"], n_void=n_void,
            umap_n_neighbors=n_neighbors, umap_min_dist=min_dist,
        )

    n_void_    = int((coverage_df["type"] == "void_probe").sum())
    n_samples_ = int((coverage_df["type"] == "sample").sum())
    st.caption(f"{n_void_:,} void probes + {n_samples_:,} real samples")

    html = plot_coverage(coverage_df, results["latents"], meta)
    components.html(html, height=700)
