import os
import random
import zarr

import numpy as np
import pandas as pd
import streamlit as st
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import ColumnDataSource
from bokeh.models import NumeralTickFormatter, PanTool, WheelZoomTool

@st.cache_data
def load_results(results_dir):
    latents         = np.load(os.path.join(results_dir, "latents.npy"))
    reconstructions = np.load(os.path.join(results_dir, "reconstructions.npy"))
    sample_ids      = np.load(os.path.join(results_dir, "sample_ids.npy"), allow_pickle=True)

    latents_df = pd.DataFrame(latents, index=sample_ids)
    latents_df.columns = [f"latent_{i+1}" for i in range(latents_df.shape[1])]

    reconstructions_df = pd.DataFrame(reconstructions, index=sample_ids)

    return {
        "latents": latents_df,
        "reconstructions": reconstructions_df
    }

@st.cache_data
def load_meta():
    meta_df = pd.read_csv(
        "../assets/Pf_9_samples_20260227.txt", index_col=0, sep = "\t",
        usecols = [
            "Sample", "Study", "Country", "Admin level 1", "Year", "Population", "% callable",
            "QC pass", "Exclusion reason", "Sample type"]
        )
    cnv_calls = pd.read_csv(
        "../assets/20260313_full_cnv_data_pf9.tsv", sep = "\t", index_col=0,
        usecols = [
            "Sample",
            "CRT_uncurated_coverage_only", "CRT_curated_coverage_only", "CRT_faceaway_only",
            "GCH1_uncurated_coverage_only", "GCH1_curated_coverage_only", "GCH1_faceaway_only",
            "MDR1_uncurated_coverage_only", "MDR1_curated_coverage_only", "MDR1_faceaway_only",
            "PM2_PM3_uncurated_coverage_only", "PM2_PM3_curated_coverage_only", "PM2_PM3_faceaway_only",
            "HRP2_uncurated_coverage_only", "HRP2_final_deletion_call",
            "HRP3_uncurated_coverage_only", "HRP3_final_deletion_call",
        ]
    )
    gff = pd.read_csv(
        "../assets/PlasmoDB-54_Pfalciparum3D7.gff",
        sep="\t", comment="#", header=None,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )

    def parse_attributes(attr_string):
        attrs = {}
        for item in attr_string.split(';'):
            if "=" in item:
                key, value = item.split('=', 1)
                attrs[key] = value
        return attrs
    
    attr_df = gff['attributes'].apply(parse_attributes).apply(pd.Series)
    gff = pd.concat([gff, attr_df], axis=1)
    gff = gff.loc[gff["type"].isin(["protein_coding_gene", "ncRNA"])].reset_index(drop = True).drop(
        columns = [
            "source", "attributes", "Note", "score", "protein_source_id",
            "strand", "phase", "Parent", "gene_id", "type"]
    ).sort_values("ID")

    return meta_df.merge(cnv_calls, left_index=True, right_index=True), gff

@st.cache_data
def load_inputs(inputs_path):
    contigs    = np.load(os.path.join(inputs_path, "contigs.npy"), allow_pickle=True)
    counts     = np.load(os.path.join(inputs_path, "counts.npy"))
    sample_ids = np.load(os.path.join(inputs_path, "sample_ids.npy"), allow_pickle=True)

    counts_df = pd.DataFrame(counts, index=sample_ids)
    contigs_df = pd.DataFrame(contigs)

    return {
        "contigs": contigs_df,
        "counts": counts_df
    }

@st.cache_data
def process_sample(contigs, sample_inputs, sample_reconstruction):
    copy_ratio = pd.DataFrame({
        "input": sample_inputs.values.flatten(),
        "reconstruction": sample_reconstruction.values.flatten()
    }, index=sample_inputs.index)
    copy_ratio["copy_ratio"] = copy_ratio["input"] / (copy_ratio["reconstruction"] + 1e-6)
    
    return pd.concat([contigs, copy_ratio], axis=1)

def plot_copy_number(data):
    chrom_options = data["chrom"].unique().tolist()

    if st.session_state.get("lucky_chrom") == "__random__":
        st.session_state["chrom_slider"] = random.choice(chrom_options)
        st.session_state["lucky_chrom"] = st.session_state["chrom_slider"]

    selected_chrom = st.select_slider("Chromosome", chrom_options, key="chrom_slider", label_visibility="hidden")

    filtered = data[data["chrom"] == selected_chrom].reset_index(drop=True)
    x       = filtered.start.values.astype(float)
    y_ratio = filtered["copy_ratio"].values.astype(float)
    y_input = filtered["input"].values.astype(float)
    y_recon = filtered["reconstruction"].values.astype(float)

    gaps    = np.where(np.diff(x) > 1000)[0] + 1
    x       = np.insert(x,       gaps, np.nan)
    y_ratio = np.insert(y_ratio, gaps, np.nan)
    y_input = np.insert(y_input, gaps, np.nan)
    y_recon = np.insert(y_recon, gaps, np.nan)

    tools = "pan,wheel_zoom,reset"
    fig_kwargs = dict(tools=tools, sizing_mode="stretch_width", output_backend="webgl")
    p1 = figure(height=250, width=700, **fig_kwargs)
    p2 = figure(height=200, width=700, x_range=p1.x_range, **fig_kwargs)

    for p in [p1, p2]:
        for tool in p.toolbar.tools:
            if isinstance(tool, PanTool):
                tool.dimensions = "width"
            elif isinstance(tool, WheelZoomTool):
                tool.dimensions = "width"

    s1 = ColumnDataSource(data=dict(x=x, y=y_ratio))
    s2 = ColumnDataSource(data=dict(x=x, input=y_input, reconstruction=y_recon))

    p1.line('x', 'y', source=s1, line_width=2)
    p1.y_range.start = -0.1
    p1.y_range.end = 5.1
    p1.xaxis.formatter = NumeralTickFormatter(format="0,0")

    p2.line('x', 'input',          source=s2, line_width=2, color="green", legend_label="input")
    p2.line('x', 'reconstruction', source=s2, line_width=2, color="red",
            legend_label="reconstruction")
    p2.xaxis.formatter = NumeralTickFormatter(format="0,0")

    return column([p1, p2], sizing_mode="stretch_width")