import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os

from src.utils import load_results, load_inputs, load_sample_meta, process_sample, plot_copy_number
from bokeh.embed import file_html
from bokeh.resources import CDN

def page1():
    st.title("First page")

    RESULTS_DIR = st.selectbox("Select results directory", options = os.listdir("../data/results/"), index=0)
    INPUTS_DIR  = st.selectbox("Select inputs directory", options = os.listdir("../data/inputs/"), index=0)

    if not RESULTS_DIR or not INPUTS_DIR:
        st.stop("Please select both a results directory and an inputs directory to proceed.")

    results = load_results(os.path.join("../data/results/", RESULTS_DIR))
    inputs  = load_inputs(os.path.join("../data/inputs/", INPUTS_DIR))
    meta    = load_sample_meta()

    st.dataframe(meta)

    SAMPLE_ID = st.selectbox("Select sample ID", options = results["latents"].index, index=0)

    data = process_sample(
        inputs["contigs"], inputs["counts"].loc[SAMPLE_ID],
        results["reconstructions"].loc[SAMPLE_ID]
    )
    
    p = plot_copy_number(data)
    p.sizing_mode = "stretch_width"
    html = file_html(p, CDN)
    components.html(html, height=300)
