import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
import random

from src.utils import load_meta, load_results, load_inputs, process_sample, plot_copy_number
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

    col1, col2 = st.columns([4, 1], vertical_alignment="bottom")
    with col1:
        SAMPLE_ID = st.selectbox("Select sample ID", options=sample_options, key="sample_select")
    with col2:
        st.button("I'm Feeling Lucky", on_click=on_lucky_click)

    data = process_sample(
        inputs["contigs"], inputs["counts"].loc[SAMPLE_ID],
        results["reconstructions"].loc[SAMPLE_ID]
    )
    
    layout = plot_copy_number(data)
    html = file_html(layout, CDN)
    components.html(html, height=500)
    
    st.dataframe(gff, hide_index = True)
