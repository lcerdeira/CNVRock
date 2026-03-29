import os
import random

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
import umap
from matplotlib.lines import Line2D
from bokeh.plotting import figure
from bokeh.layouts import column, row as bk_row
from bokeh.models import (ColumnDataSource, HoverTool, ColorBar,
                           LinearColorMapper, CustomJS, Span)
from bokeh.models import NumeralTickFormatter, PanTool, WheelZoomTool
from bokeh.palettes import Plasma256

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

@st.cache_data
def compute_pca(latents_df):
    pca = PCA(n_components=2)
    coords = pca.fit_transform(latents_df.values)
    df = pd.DataFrame(coords, index=latents_df.index, columns=["PC1", "PC2"])
    variance = pca.explained_variance_ratio_ * 100
    return df, variance

_KDE_CLASSES = {"gDNA": "green", "sWGA": "blue"}

@st.cache_data
def compute_pca_contours(pca_df, meta):
    """Compute KDE contour paths per sample class. Cached — independent of selected sample."""
    pca_with_type = pca_df.join(meta[["Sample type"]], how="left")

    pad_x = (pca_df["PC1"].max() - pca_df["PC1"].min()) * 0.15
    pad_y = (pca_df["PC2"].max() - pca_df["PC2"].min()) * 0.15
    gx = np.linspace(pca_df["PC1"].min() - pad_x, pca_df["PC1"].max() + pad_x, 100)
    gy = np.linspace(pca_df["PC2"].min() - pad_y, pca_df["PC2"].max() + pad_y, 100)
    xx, yy = np.meshgrid(gx, gy)
    grid_pts = np.vstack([xx.ravel(), yy.ravel()])

    result = {}
    for stype, color in _KDE_CLASSES.items():
        subset = pca_with_type[pca_with_type["Sample type"] == stype]
        if len(subset) < 5:
            continue
        kde = gaussian_kde(subset[["PC1", "PC2"]].values.T)
        zz  = kde(grid_pts).reshape(xx.shape)
        # Use percentiles of the grid density so levels are guaranteed within zz range
        levels = np.percentile(zz, [96, 98, 99.8]).tolist()

        fig, ax = plt.subplots()
        cs = ax.contour(gx, gy, zz, levels=sorted(levels))
        paths = [path.vertices for path in cs.get_paths()]
        plt.close(fig)

        result[stype] = (color, paths)
    return result

def plot_latents(latent_values):
    vals   = latent_values.values.astype(float)
    labels = [f"z{i+1}" for i in range(len(vals))]
    abs_max = np.abs(vals).max() or 1.0

    norm   = plt.Normalize(vmin=0, vmax=abs_max)
    colors = plt.cm.Reds(norm(np.abs(vals)))

    fig, ax = plt.subplots(figsize=(8, 1.8))
    bars = ax.bar(labels, vals, color=colors, edgecolor="none")
    ax.axhline(0, color="black", linewidth=0.6)

    pad = abs_max * 0.04
    for bar, val in zip(bars, vals):
        x  = bar.get_x() + bar.get_width() / 2
        if val >= 0:
            ax.text(x, val + pad, f"{val:.1f}", ha="center", va="bottom", fontsize=7)
        else:
            ax.text(x, val - pad, f"{val:.1f}", ha="center", va="top",    fontsize=7)

    ax.set_xlim(-0.5, len(vals) - 0.5)
    ax.set_ylim(-abs_max * 1.25, abs_max * 1.25)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("latent value", fontsize=8)
    fig.tight_layout()
    return fig

def plot_pca(pca_df, variance, contours, selected_sample):
    fig, ax = plt.subplots(figsize=(4, 4))

    other    = pca_df[pca_df.index != selected_sample]
    selected = pca_df[pca_df.index == selected_sample]

    for color, paths in contours.values():
        for verts in paths:
            ax.plot(verts[:, 0], verts[:, 1], color=color, alpha=0.7, linewidth=1.5)

    ax.scatter(other["PC1"],    other["PC2"],    c="grey", s=8,  alpha=0.4, zorder=-1)
    ax.scatter(selected["PC1"], selected["PC2"], c="red",  s=40, alpha=1.0, zorder=3)

    handles = [
        Line2D([0], [0], color=color, linewidth=1.5, label=stype)
        for stype, (color, _) in contours.items()
    ]
    ax.legend(handles=handles, framealpha=0.7, fontsize=8)

    ax.set_xlabel(f"PC1 ({variance[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1]:.1f}%)")
    fig.tight_layout()

    return fig

@st.cache_data
def compute_umap_coverage(latents_df, n_probes=300_000, empty_pct=95, umap_n_neighbors=30, umap_min_dist=0.1):
    z   = latents_df.values
    rng = np.random.default_rng(42)

    lo, hi = z.min(axis=0), z.max(axis=0)
    probes = rng.uniform(lo, hi, size=(n_probes, z.shape[1]))

    tree = cKDTree(z)
    dists, _ = tree.query(probes, k=1, workers=-1)
    threshold    = np.percentile(dists, empty_pct)
    empty_probes = probes[dists > threshold]
    empty_dists  = dists[dists > threshold]

    reducer   = umap.UMAP(n_components=2, n_neighbors=umap_n_neighbors,
                          min_dist=umap_min_dist, random_state=42, low_memory=True)
    embedding = reducer.fit_transform(empty_probes)

    df             = pd.DataFrame(embedding, columns=["u1", "u2"])
    df["nn_dist"]  = empty_dists
    df[list(latents_df.columns)] = empty_probes
    return df


def plot_umap_coverage(umap_df, latents_df):
    lat_cols = list(latents_df.columns)
    n_lats   = len(lat_cols)
    x_labels = [f"z{i+1}" for i in range(n_lats)]

    y_abs_max = float(np.abs(umap_df[lat_cols].values).max()) * 1.15

    src = ColumnDataSource({
        "u1":      umap_df["u1"].values,
        "u2":      umap_df["u2"].values,
        "nn_dist": umap_df["nn_dist"].values,
        **{c: umap_df[c].values for c in lat_cols},
    })

    bar_src = ColumnDataSource({
        "x":     x_labels,
        "top":   [0.0] * n_lats,
        "color": ["#eeeeee"] * n_lats,
    })

    # --- Scatter ---
    color_mapper = LinearColorMapper(
        palette=Plasma256,
        low=float(umap_df["nn_dist"].min()),
        high=float(umap_df["nn_dist"].max()),
    )
    p_scatter = figure(
        width=720, height=640,
        output_backend="webgl",
        tools="pan,wheel_zoom,box_zoom,reset",
        title="UMAP of empty latent regions — hover to inspect",
    )
    p_scatter.scatter("u1", "u2", source=src, size=4, alpha=0.55,
                      color={"field": "nn_dist", "transform": color_mapper})
    p_scatter.add_layout(
        ColorBar(color_mapper=color_mapper, label_standoff=8, width=12,
                 title="dist to nearest sample"),
        "right"
    )
    p_scatter.xaxis.axis_label = "UMAP 1"
    p_scatter.yaxis.axis_label = "UMAP 2"

    # --- Bar chart ---
    p_bar = figure(
        width=360, height=640,
        x_range=x_labels,
        y_range=(-y_abs_max, y_abs_max),
        title="Latent values at cursor",
        tools="",
    )
    p_bar.vbar(x="x", top="top", bottom=0, width=0.75, color="color", source=bar_src)
    p_bar.add_layout(Span(location=0, dimension="width", line_color="black", line_width=1))
    p_bar.xaxis.major_label_orientation = 0.5
    p_bar.yaxis.axis_label = "latent value"
    p_bar.xgrid.grid_line_color = None

    # --- CustomJS hover callback (pure client-side — no round trips) ---
    callback = CustomJS(args=dict(source=src, bar_source=bar_src, lat_cols=lat_cols), code="""
        const indices = cb_data.index.indices;
        if (!indices.length) return;
        const i = indices[0];

        const vals = lat_cols.map(c => source.data[c][i]);
        const abs_max = Math.max(...vals.map(v => Math.abs(v)));

        const colors = vals.map(v => {
            if (abs_max === 0) return '#eeeeee';
            const t = Math.abs(v) / abs_max;
            // white (255,255,255) -> crimson (220,20,60)
            const r = Math.round(255 - t * 35);
            const g = Math.round(255 - t * 235);
            const b = Math.round(255 - t * 195);
            return `rgb(${r},${g},${b})`;
        });

        bar_source.data['top']   = vals;
        bar_source.data['color'] = colors;
        bar_source.change.emit();
    """)
    p_scatter.add_tools(HoverTool(callback=callback, tooltips=None, mode="mouse"))

    return bk_row([p_scatter, p_bar])


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