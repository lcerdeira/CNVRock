import os
import streamlit as st
import streamlit.components.v1 as components

from src.utils import load_results, load_meta, compute_coverage, plot_coverage


def page2():
    st.title("Latent space coverage")

    RESULTS_DIR = st.selectbox("Select results directory", options=os.listdir("../data/results/"), index=0)
    if not RESULTS_DIR:
        st.stop("Please select a results directory.")

    results    = load_results(os.path.join("../data/results/", RESULTS_DIR))
    meta, _gff = load_meta()

    with st.expander("Settings", expanded=False):
        n_void      = st.slider("Void probes",    5_000, 50_000, 20_000, step=5_000)
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
