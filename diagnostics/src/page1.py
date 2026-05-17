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
    if not experiments:
        st.warning(
            "No experiments found. Either you're on Streamlit Cloud demo mode "
            "and only the bundled subsample exists (use the Monitor page), or "
            "the `models/experiments/` directory hasn't been populated."
        )
        st.stop()
    EXPERIMENT  = st.selectbox("Experiment", options=experiments, index=0)

    if not EXPERIMENT:
        st.warning("Please select an experiment to proceed.")
        st.stop()

    cfg = load_experiment_config(EXPERIMENT)

    results = load_results(cfg["out_dir"])
    inputs  = load_inputs(cfg["store_path"])
    meta, gff = load_meta(cfg.get("kpsc_meta_path"))

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
        sample_options = list(results["latents"].index)

    # --- Gene call filter (plasmid + chromosomal) ----------------------------
    _CN_STATE_LABELS = {-1: "uncallable", 0: "absent / deletion", 1: "normal / present", 2: "amplified"}

    with st.expander("Filter by gene calls", expanded=False):
        gene_filtered_ids = set(sample_options)

        plasmid_calls = results.get("plasmid_calls")
        if plasmid_calls is not None:
            st.markdown("**Plasmid genes** (cn: 0 = absent, 1 = present, 2 = amplified)")
            p_genes = [c for c in plasmid_calls.columns if not c.startswith("pcn_")]
            cols = st.columns(min(4, len(p_genes)))
            for i, gene in enumerate(p_genes):
                with cols[i % len(cols)]:
                    available_states = sorted(plasmid_calls[gene].dropna().unique().astype(int).tolist())
                    state_labels = [_CN_STATE_LABELS.get(s, str(s)) for s in available_states]
                    selected = st.multiselect(
                        gene, options=available_states,
                        format_func=lambda s: _CN_STATE_LABELS.get(s, str(s)),
                        key=f"pfilter_{EXPERIMENT}_{gene}",
                        placeholder="any",
                    )
                    if selected:
                        keep = set(plasmid_calls.index[plasmid_calls[gene].isin(selected)].tolist())
                        gene_filtered_ids &= keep

        chrom_calls = results.get("gene_calls")
        if chrom_calls is not None:
            st.markdown("**Chromosomal genes** (cn: 0 = deletion, 1 = normal, 2 = amplified, -1 = uncallable)")
            c_genes = [c for c in chrom_calls.columns if not c.startswith("crr_")]
            cols = st.columns(min(4, len(c_genes)))
            for i, gene in enumerate(c_genes):
                with cols[i % len(cols)]:
                    available_states = sorted(chrom_calls[gene].dropna().unique().astype(int).tolist())
                    selected = st.multiselect(
                        gene, options=available_states,
                        format_func=lambda s: _CN_STATE_LABELS.get(s, str(s)),
                        key=f"cfilter_{EXPERIMENT}_{gene}",
                        placeholder="any",
                    )
                    if selected:
                        keep = set(chrom_calls.index[chrom_calls[gene].isin(selected)].tolist())
                        gene_filtered_ids &= keep

        n_before = len(sample_options)
        sample_options = [s for s in sample_options if s in gene_filtered_ids]
        if len(sample_options) < n_before:
            st.caption(f"{len(sample_options)} / {n_before} samples match gene filter")
        if not sample_options:
            st.warning("No samples match — gene filter cleared.")
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

    # Sample metadata badge
    if SAMPLE_ID and not meta.empty and SAMPLE_ID in meta.index:
        row    = meta.loc[SAMPLE_ID]
        badges = []
        if "Species" in meta.columns:
            sp = str(row.get("Species", ""))
            if sp and sp != "nan":
                badges.append(f"**Species:** {sp}")
        if "ST" in meta.columns:
            st_val = str(row.get("ST", ""))
            if st_val and st_val not in ("nan", "Unknown", "-"):
                badges.append(f"**ST:** {st_val}")
        if badges:
            st.caption("  ·  ".join(badges))

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

        fig = plot_pca(pca_df, variance, contours, SAMPLE_ID, meta=meta)
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

    if not gff.empty:
        @st.dialog("Gene annotations", width="large")
        def _show_gff(chrom):
            st.dataframe(gff[gff["seqid"] == chrom], hide_index=True, width="stretch")

        if st.button("Gene annotations"):
            _show_gff(st.session_state.get("chrom_slider", data["chrom"].iloc[0]))
