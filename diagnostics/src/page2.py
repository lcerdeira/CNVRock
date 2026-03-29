import os
import streamlit as st
import streamlit.components.v1 as components

from bokeh.embed import file_html
from bokeh.resources import CDN

from src.utils import load_results, compute_umap_coverage, plot_umap_coverage


def page2():
    st.title("Latent space coverage")

    RESULTS_DIR = st.selectbox("Select results directory", options=os.listdir("../data/results/"), index=0)
    if not RESULTS_DIR:
        st.stop("Please select a results directory.")

    results = load_results(os.path.join("../data/results/", RESULTS_DIR))

    with st.expander("Settings", expanded=False):
        n_probes    = st.slider("Random probes",          50_000, 500_000, 300_000, step=50_000)
        empty_pct   = st.slider("Empty threshold (percentile)", 80, 99, 95)
        n_neighbors = st.slider("UMAP n_neighbors",       5, 100, 30)
        min_dist    = st.slider("UMAP min_dist",          0.0, 0.5, 0.1, step=0.05)

    with st.spinner("Running UMAP on empty probes (cached after first run)…"):
        umap_df = compute_umap_coverage(
            results["latents"],
            n_probes=n_probes,
            empty_pct=empty_pct,
            umap_n_neighbors=n_neighbors,
            umap_min_dist=min_dist,
        )

    st.caption(f"{len(umap_df):,} empty probes embedded")

    layout = plot_umap_coverage(umap_df, results["latents"])
    components.html(file_html(layout, CDN), height=660)
