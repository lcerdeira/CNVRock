import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
import random

import matplotlib.pyplot as plt

from src.utils import (load_meta, load_results, load_inputs, process_sample,
                       compute_pca, compute_pca_contours, plot_latents, plot_pca,
                       plot_copy_number, fit_hmm_sample, call_all_genes)
from bokeh.embed import file_html
from bokeh.resources import CDN

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

    def on_lucky_click():
        st.session_state["sample_select"] = random.choice(sample_options)
        st.session_state["lucky_chrom"] = "__random__"

    def on_random_sample_click():
        st.session_state["sample_select"] = random.choice(sample_options)

    col1, col2, col3 = st.columns([4, 1, 1], vertical_alignment="bottom")
    with col1:
        SAMPLE_ID = st.selectbox("Select sample ID", options=sample_options, key="sample_select")
    with col2:
        st.button("Random sample", on_click=on_random_sample_click, width="stretch")
    with col3:
        st.button("I'm Feeling Lucky", on_click=on_lucky_click, width="stretch")

    pca_df, variance = compute_pca(results["latents"])
    contours = compute_pca_contours(pca_df, meta)

    data = process_sample(
        inputs["contigs"], inputs["counts"].loc[SAMPLE_ID],
        results["reconstructions"].loc[SAMPLE_ID]
    )

    col_pca, col_cn = st.columns([1, 3])
    with col_pca:
        lat_fig = plot_latents(results["latents"].loc[SAMPLE_ID])
        st.pyplot(lat_fig, width="stretch")
        plt.close(lat_fig)

        fig = plot_pca(pca_df, variance, contours, SAMPLE_ID)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
    with col_cn:
        precomputed = results["segments"]
        if precomputed is not None:
            sample_segs = precomputed[precomputed["sample_id"] == SAMPLE_ID]
        else:
            with st.expander("HMM parameters", expanded=False):
                n_states          = st.slider("CN states",              3, 8,    6)
                self_transition   = st.slider("Self-transition prob",   0.80, 0.999, 0.95,
                                              step=0.005, format="%.3f")
                low_cov_threshold = st.slider("Low-coverage threshold", 0, 100, 10)
            with st.spinner("Fitting HMM…"):
                sample_segs = fit_hmm_sample(
                    data,
                    n_states=n_states,
                    self_transition=self_transition,
                    low_cov_threshold=low_cov_threshold,
                )
        cn_layout = plot_copy_number(data, sample_segs)
        components.html(file_html(cn_layout, CDN), height=520)

    gene_calls = call_all_genes(data, sample_segs)
    st.dataframe(pd.DataFrame(gene_calls), hide_index=True, width="stretch")

    @st.dialog("Gene annotations", width="large")
    def _show_gff(chrom):
        st.dataframe(gff[gff["seqid"] == chrom], hide_index=True, width="stretch")

    if st.button("Gene annotations"):
        _show_gff(st.session_state.get("chrom_slider", data["chrom"].iloc[0]))
