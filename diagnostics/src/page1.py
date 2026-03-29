import os

import streamlit as st
import streamlit.components.v1 as components

from src.utils import (load_meta, load_results, load_inputs, process_sample,
                       compute_pca, compute_pca_contours, plot_page1_html)

def page1():
    st.title("First page")

    RESULTS_DIR = st.selectbox("Select results directory", options = os.listdir("../data/results/"), index=0)
    INPUTS_DIR  = st.selectbox("Select inputs directory", options = os.listdir("../data/inputs/"), index=0)

    if not RESULTS_DIR or not INPUTS_DIR:
        st.stop("Please select both a results directory and an inputs directory to proceed.")

    results   = load_results(os.path.join("../data/results/", RESULTS_DIR))
    inputs    = load_inputs(os.path.join("../data/inputs/", INPUTS_DIR))
    meta, gff = load_meta()

    st.dataframe(meta)

    sample_options = list(results["latents"].index)
    SAMPLE_ID = st.selectbox(
        "Coverage sample", options=sample_options, key="sample_select",
        help="Changes which sample's coverage is loaded. "
             "Sample selection inside the plot updates PCA and latents instantly.",
    )

    pca_df, variance = compute_pca(results["latents"])
    contours = compute_pca_contours(pca_df, meta)

    data = process_sample(
        inputs["contigs"], inputs["counts"].loc[SAMPLE_ID],
        results["reconstructions"].loc[SAMPLE_ID]
    )

    html = plot_page1_html(pca_df, variance, contours, results["latents"], data, SAMPLE_ID)
    components.html(html, height=620)

    st.dataframe(gff, hide_index=True)
