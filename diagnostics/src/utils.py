import os
import zarr

import numpy as np
import pandas as pd
import streamlit as st
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.resources import CDN
from bokeh.embed import components

@st.cache_data
def load_results(results_dir):
    latents         = np.load(os.path.join(results_dir, "latents.npy"))
    reconstructions = np.load(os.path.join(results_dir, "reconstructions.npy"))
    sample_ids      = np.load(os.path.join(results_dir, "sample_ids.npy"))

    latents_df = pd.DataFrame(latents, index=sample_ids)
    latents_df.columns = [f"latent_{i+1}" for i in range(latents_df.shape[1])]

    reconstructions_df = pd.DataFrame(reconstructions, index=sample_ids)

    return {
        "latents": latents_df,
        "reconstructions": reconstructions_df
    }

@st.cache_data
def load_sample_meta():
    meta_df = pd.read_csv(
        "../assets/Pf_9_samples_20260227.txt", index_col=0, sep = "\t",
        usecols = [
            "Sample", "Study", "Country", "Admin level 1", "Year", "Population", "% callable",
            "QC pass", "Exclusion reason", "Sample type"]
        )
    cnv_calls = pd.read_csv(
        "../assets/20260311_full_cnv_data_pf9.tsv", sep = "\t", index_col=0,
        usecols = [
            "CRT_curated_coverage_only", "CRT_faceaway_only",
            "GCH1_curated_coverage_only", "GCH1_faceaway_only",
            "MDR1_curated_coverage_only", "MDR1_faceaway_only",
            "PM2_PM3_curated_coverage_only", "PM2_PM3_faceaway_only",
            "HRP2_uncurated_coverage_only", "HRP2_final_deletion_call",
            "HRP3_uncurated_coverage_only", "HRP3_final_deletion_call",
        ]
    )
    return meta_df

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
    selected_chrom = st.select_slider("Chromosome", chrom_options, label_visibility = "hidden")
    
    filtered = data[data["chrom"] == selected_chrom].reset_index(drop=True)
    x = np.arange(len(filtered))
    source = ColumnDataSource(data=dict(x=x, y=filtered["copy_ratio"].values))
    
    p = figure(height=250, width=700, title=f"Copy Number Profile — {selected_chrom}",
               x_axis_label="Bin", y_axis_label="Copy Ratio",
               tools="pan,wheel_zoom,box_zoom,reset,save")
    p.line('x', 'y', source=source, line_width=2)
    p.circle('x', 'y', source=source, size=3, color="navy", alpha=0.5)
    p.y_range.start = -0.1
    p.y_range.end = 5.1
    return p
