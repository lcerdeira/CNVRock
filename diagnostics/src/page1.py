import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import random

import matplotlib.pyplot as plt

from src.utils import (load_meta, load_results, load_inputs, process_sample,
                       compute_pca, compute_pca_contours, plot_latents, plot_pca,
                       plot_copy_number,
                       list_experiments, load_experiment_config,
                       fit_hmm_sample_versioned, call_all_genes_versioned)
from bokeh.embed import file_html
from bokeh.resources import CDN

def page1():
    st.title("First page")

    experiments = list_experiments()
    EXPERIMENT  = st.selectbox("Experiment", options=experiments, index=0)

    if not EXPERIMENT:
        st.stop("Please select an experiment to proceed.")

    cfg = load_experiment_config(EXPERIMENT)

    results = load_results(cfg["out_dir"])
    inputs  = load_inputs(cfg["store_path"])
    meta, gff = load_meta()

    # --- Sample filter ---------------------------------------------------------
    filter_key = f"meta_filter_{EXPERIMENT}"

    with st.expander("Filter samples", expanded=False):
        filter_text = st.text_area(
            "One pandas query condition per line (lines are AND-ed together)",
            value=st.session_state.get(filter_key, ""),
            height=120,
            placeholder="Sample_type == 'aAMP'\nCountry == 'Ghana'",
            key=f"filter_text_{EXPERIMENT}",
        )
        col_apply, col_clear = st.columns([1, 1])
        if col_apply.button("Apply filter", key=f"apply_{EXPERIMENT}"):
            st.session_state[filter_key] = filter_text
        if col_clear.button("Clear filter", key=f"clear_{EXPERIMENT}"):
            st.session_state[filter_key] = ""
            st.rerun()

    active_filter = st.session_state.get(filter_key, "")
    filtered_meta = meta
    if active_filter.strip():
        conditions = [ln.strip() for ln in active_filter.splitlines() if ln.strip()]
        combined   = " & ".join(f"({c})" for c in conditions)
        try:
            filtered_meta = meta.query(combined)
        except Exception as e:
            st.error(f"Filter error: {e}")

    st.dataframe(filtered_meta)

    # Restrict sample options to those present in the filtered meta
    filtered_ids   = set(filtered_meta.index)
    sample_options = [s for s in results["latents"].index if s in filtered_ids]
    if not sample_options:
        st.warning("No samples match the current filter — showing all samples.")
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
            with st.spinner("Fitting HMM…"):
                sample_segs = fit_hmm_sample_versioned(
                    cfg["hmm"], data,
                    n_states          = cfg["hmm_n_states"],
                    self_transition   = cfg["hmm_self_transition"],
                    low_cov_threshold = cfg["hmm_low_cov_threshold"],
                )
        cn_layout = plot_copy_number(data, sample_segs)
        components.html(file_html(cn_layout, CDN), height=520)

    precomputed_calls = results["gene_calls"]
    if precomputed_calls is not None and SAMPLE_ID in precomputed_calls.index:
        gene_calls = precomputed_calls.loc[[SAMPLE_ID]].to_dict(orient="records")
    else:
        gene_calls = call_all_genes_versioned(
            cfg["cnv"], data, sample_segs,
            min_cn1_proportion    = cfg["cnv_min_cn1_proportion"],
            min_confidence        = cfg["cnv_min_confidence"],
            flank_padding         = cfg["cnv_flank_padding"],
            crr_amp_threshold     = cfg["cnv_crr_amp_threshold"],
            crr_min_bins_fallback = cfg["cnv_crr_min_bins_fallback"],
        )
    st.dataframe(pd.DataFrame(gene_calls), hide_index=True, width="stretch")

    @st.dialog("Gene annotations", width="large")
    def _show_gff(chrom):
        st.dataframe(gff[gff["seqid"] == chrom], hide_index=True, width="stretch")

    if st.button("Gene annotations"):
        _show_gff(st.session_state.get("chrom_slider", data["chrom"].iloc[0]))
