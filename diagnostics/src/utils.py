import importlib
import os
import random
import sys

import numpy as np
import yaml
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
                           LinearColorMapper, CustomJS, Span, BoxAnnotation)
from bokeh.models import NumeralTickFormatter, PanTool, WheelZoomTool, LabelSet
from bokeh.palettes import Plasma256

@st.cache_data
def load_results(results_dir):
    latents         = np.load(os.path.join(results_dir, "latents.npy"))
    reconstructions = np.load(os.path.join(results_dir, "reconstructions.npy"))
    sample_ids      = np.load(os.path.join(results_dir, "sample_ids.npy"), allow_pickle=True)

    latents_df = pd.DataFrame(latents, index=sample_ids)
    latents_df.columns = [f"latent_{i+1}" for i in range(latents_df.shape[1])]

    reconstructions_df = pd.DataFrame(reconstructions, index=sample_ids)

    segs_path = os.path.join(results_dir, "segments.parquet")
    segments  = pd.read_parquet(segs_path) if os.path.exists(segs_path) else None

    return {
        "latents": latents_df,
        "reconstructions": reconstructions_df,
        "segments": segments,
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

    return meta_df.merge(cnv_calls, left_index=True, right_index=True, how="left"), gff

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
def compute_coverage(latents_df, n_void=20_000, umap_n_neighbors=15, umap_min_dist=0.05):
    z        = latents_df.values
    lat_cols = list(latents_df.columns)
    n_dims   = z.shape[1]
    rng      = np.random.default_rng(42)

    sample_tree   = cKDTree(z)
    sample_dists  = sample_tree.query(z, k=2,   p=1, workers=-1)[0][:, 1]
    cluster_dists = sample_tree.query(z, k=501, p=1, workers=-1)[0][:, -1]
    cluster_radius = float(np.median(cluster_dists))  # ~7 — inter-cluster scale

    # Perturb real samples with Laplace noise targeting cluster-scale displacement
    # E[L1 norm] = n_dims * noise_scale = 2 * cluster_radius
    noise_scale = 2.0 * cluster_radius / n_dims
    K           = 4  # perturbations per sample → ~215k candidates
    anchors     = np.repeat(z, K, axis=0)
    noise       = rng.laplace(loc=0.0, scale=noise_scale, size=anchors.shape)
    candidates  = anchors + noise

    cand_dists, _ = sample_tree.query(candidates, k=1, p=1, workers=-1)

    # Keep probes in the inter-cluster void zone
    mask        = (cand_dists > cluster_radius) & (cand_dists < 3.0 * cluster_radius)
    probes      = candidates[mask]
    probe_dists = cand_dists[mask]

    if len(probes) < n_void:
        raise ValueError(f"Only {len(probes)} probe candidates survived filtering (need {n_void})")

    # Take the n_void closest to the boundary
    top_idx     = np.argpartition(probe_dists, n_void)[:n_void]
    probes      = probes[top_idx]
    probe_dists = probe_dists[top_idx]

    # PCA to 6D first, then UMAP
    combined     = np.vstack([probes, z])
    combined_pca = PCA(n_components=6).fit_transform(combined)

    reducer   = umap.UMAP(n_components=2, n_neighbors=umap_n_neighbors,
                          min_dist=umap_min_dist, random_state=42, low_memory=True)
    embedding = reducer.fit_transform(combined_pca)

    n_probes_ = len(probes)

    probe_df            = pd.DataFrame(embedding[:n_probes_], columns=["u1", "u2"])
    probe_df["nn_dist"] = probe_dists
    probe_df["type"]    = "void_probe"
    probe_df[lat_cols]  = probes

    sample_df            = pd.DataFrame(embedding[n_probes_:], columns=["u1", "u2"])
    sample_df["nn_dist"] = sample_dists
    sample_df["type"]    = "sample"
    sample_df[lat_cols]  = z
    sample_df.index      = latents_df.index

    return pd.concat([probe_df, sample_df])


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


_SAMPLE_TYPE_COLORS = {
    "gDNA":    "rgba(80,210,120,0.6)",
    "sWGA":    "rgba(80,170,255,0.6)",
    "aAMP":    "rgba(255,160,60,0.6)",
    "Unknown": "rgba(180,180,180,0.5)",
}

def plot_coverage(coverage_df, latents_df, meta):
    import json

    lat_cols  = list(latents_df.columns)
    n_lats    = len(lat_cols)
    x_labels  = [f"z{i+1}" for i in range(n_lats)]
    y_abs_max = float(np.abs(coverage_df[lat_cols].values).max()) * 1.15

    probe_df  = coverage_df[coverage_df["type"] == "void_probe"]
    sample_df = coverage_df[coverage_df["type"] == "sample"].copy()
    sample_df = sample_df.join(meta[["Sample type"]], how="left")
    sample_df["Sample type"] = sample_df["Sample type"].fillna("Unknown")

    probe_payload = json.dumps({
        "x":          probe_df["u1"].tolist(),
        "y":          probe_df["u2"].tolist(),
        "nn_dist":    probe_df["nn_dist"].tolist(),
        "customdata": probe_df[lat_cols].values.tolist(),
    })

    ordered_types = [t for t in _SAMPLE_TYPE_COLORS if t in sample_df["Sample type"].values]
    other_types   = [t for t in sample_df["Sample type"].unique() if t not in _SAMPLE_TYPE_COLORS]
    type_payloads = {}
    for stype in ordered_types + other_types:
        sub = sample_df[sample_df["Sample type"] == stype]
        type_payloads[stype] = {
            "x":          sub["u1"].tolist(),
            "y":          sub["u2"].tolist(),
            "customdata": sub[lat_cols].values.tolist(),
            "color":      _SAMPLE_TYPE_COLORS.get(stype, "rgba(180,180,180,0.5)"),
        }

    return f"""<!DOCTYPE html>
<html>
<head>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
  html, body {{ margin: 0; height: 100%; background: #0e1117; }}
  #wrap {{ display: flex; height: 100vh; gap: 0; }}
  #scatter {{ flex: 3; min-width: 0; }}
  #bar     {{ flex: 1; min-width: 0; }}
</style>
</head>
<body>
<div id="wrap"><div id="scatter"></div><div id="bar"></div></div>
<script>
const probes      = {probe_payload};
const typeData    = {json.dumps(type_payloads)};
const X_LABELS    = {json.dumps(x_labels)};
const N_LATS      = {n_lats};
const Y_ABS_MAX   = {y_abs_max};

const DARK  = '#0e1117';
const PANEL = '#1a1f2e';
const FONT  = {{color: '#e0e0e0'}};

const sampleTraces = Object.entries(typeData).map(([name, d]) => ({{
    type: 'scattergl', mode: 'markers',
    name,
    x: d.x, y: d.y,
    customdata: d.customdata,
    marker: {{size: 3, color: d.color}},
    hovertemplate: `<extra>${{name}}</extra>`,
}}));

const probeTrace = {{
    type: 'scattergl', mode: 'markers',
    name: 'Void probes',
    x: probes.x, y: probes.y,
    customdata: probes.customdata,
    marker: {{
        size: 3, opacity: 0.5,
        color: probes.nn_dist, colorscale: 'Plasma',
        colorbar: {{
            title: {{text: 'dist to<br>nearest<br>sample', font: FONT}},
            tickfont: FONT, thickness: 14, len: 0.5, x: 1.02,
        }},
    }},
    hovertemplate: 'dist: %{{marker.color:.3f}}<extra>Void probe</extra>',
}};

Plotly.newPlot('scatter', [...sampleTraces, probeTrace], {{
    paper_bgcolor: DARK, plot_bgcolor: PANEL,
    xaxis: {{title: 'UMAP 1', color: '#aaa', gridcolor: '#333', zerolinecolor: '#444'}},
    yaxis: {{title: 'UMAP 2', color: '#aaa', gridcolor: '#333', zerolinecolor: '#444'}},
    legend: {{font: FONT, bgcolor: 'rgba(30,35,50,0.85)', bordercolor: '#444', borderwidth: 1}},
    margin: {{l:50, r:80, t:40, b:50}},
    title: {{text: 'Latent space coverage — click legend to toggle', font: {{...FONT, size: 13}}}},
    hovermode: 'closest',
}}, {{responsive: true}});

Plotly.newPlot('bar', [{{
    type: 'bar',
    x: X_LABELS,
    y: new Array(N_LATS).fill(0),
    marker: {{color: new Array(N_LATS).fill('#333')}},
}}], {{
    paper_bgcolor: DARK, plot_bgcolor: PANEL,
    font: FONT,
    yaxis: {{range: [-Y_ABS_MAX, Y_ABS_MAX], gridcolor: '#333', zerolinecolor: '#888',
             title: {{text: 'latent value', font: FONT}}}},
    xaxis: {{gridcolor: '#333', tickfont: {{size: 11}}}},
    margin: {{l:55, r:20, t:50, b:50}},
    title: {{text: 'Latent profile', font: {{...FONT, size: 14}}}},
    bargap: 0.25,
}}, {{responsive: true}});

function updateBar(vals) {{
    const amax = Math.max(...vals.map(v => Math.abs(v)));
    const colors = vals.map(v => {{
        if (amax === 0) return '#555';
        const t = Math.abs(v) / amax;
        return `rgb(${{Math.round(255-t*35)}},${{Math.round(255-t*235)}},${{Math.round(255-t*195)}})`;
    }});
    Plotly.restyle('bar', {{y: [vals], 'marker.color': [colors]}}, 0);
}}

document.getElementById('scatter').on('plotly_hover', function(ev) {{
    updateBar(ev.points[0].customdata);
}});
</script>
</body>
</html>"""


# CN 0 → 5: deep blue, light blue, grey, orange, red, dark red
_CN_COLORS = ["#313695", "#74add1", "#888888", "#fdae61", "#d73027", "#a50026"]


# ---------------------------------------------------------------------------
# Versioned model loading
# ---------------------------------------------------------------------------

# models/ lives two levels up from diagnostics/src/
_MODELS_DIR     = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../models"))
_EXPERIMENTS_DIR = os.path.join(_MODELS_DIR, "experiments")

if _MODELS_DIR not in sys.path:
    sys.path.insert(0, _MODELS_DIR)


def list_experiments():
    """Return sorted list of experiment folder names that contain a config.yaml."""
    return sorted([
        d for d in os.listdir(_EXPERIMENTS_DIR)
        if d[0].isdigit()
        and os.path.isfile(os.path.join(_EXPERIMENTS_DIR, d, "config.yaml"))
    ])


@st.cache_data
def load_experiment_config(experiment_id):
    """Load the config.yaml for the given experiment and resolve all paths to absolute."""
    config_dir  = os.path.join(_EXPERIMENTS_DIR, experiment_id)
    config_path = os.path.join(config_dir, "config.yaml")
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    # Resolve relative paths against the config file's own directory
    for key in ("store_path", "out_dir", "pf9_gt_path", "pf9_meta_path"):
        if key in cfg and not os.path.isabs(cfg[key]):
            cfg[key] = os.path.normpath(os.path.join(config_dir, cfg[key]))
    return cfg


@st.cache_data
def fit_hmm_sample_versioned(hmm_version, data, n_states, self_transition, low_cov_threshold):
    """Load fit_hmm_sample from the named HMM version and run it.

    Results are cached by Streamlit — changing any argument invalidates the cache.
    """
    fn = importlib.import_module(f"hmm.{hmm_version}").fit_hmm_sample
    return fn(data, n_states, self_transition, low_cov_threshold)


def call_all_genes_versioned(cnv_version, data, segments,
                             min_cn1_proportion, min_confidence, flank_padding):
    """Load call_all_genes from the named CNV version and run it."""
    fn = importlib.import_module(f"cnv.{cnv_version}").call_all_genes
    return fn(data, segments, min_cn1_proportion, min_confidence, flank_padding)


def _confidence_color(conf):
    """Red (low confidence) → black (high confidence)."""
    r = int(215 * (1 - conf))
    g = int(48  * (1 - conf))
    b = int(39  * (1 - conf))
    return f"#{r:02x}{g:02x}{b:02x}"


def _render_segments(segs):
    """Add rendering columns to a pre-computed segments DataFrame."""
    segs = segs.copy()
    segs["color"]      = segs["cn"].map(lambda cn: _CN_COLORS[min(cn, len(_CN_COLORS) - 1)])
    segs["conf_color"] = segs["confidence"].map(_confidence_color)
    segs["label_x"]    = (segs["x0"] + segs["x1"]) / 2
    segs["label"]      = segs["confidence"].map(lambda c: f"{c:.2f}")
    return segs


def plot_copy_number(data, segments=None):
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

    # Load pre-computed HMM segments for this chromosome
    if segments is not None and len(segments) > 0:
        chrom_segs = segments[segments["chrom"] == selected_chrom]
        segs = _render_segments(chrom_segs) if len(chrom_segs) > 0 else pd.DataFrame()
    else:
        segs = pd.DataFrame()

    # Low-coverage vrects: bins where input or reconstruction < 10
    low_cov = (filtered["input"].values < 10) | (filtered["reconstruction"].values < 10)
    low_starts, low_ends = [], []
    in_run = False
    for i, flag in enumerate(low_cov):
        if flag and not in_run:
            low_starts.append(float(filtered["start"].iloc[i]))
            in_run = True
        elif not flag and in_run:
            low_ends.append(float(filtered["start"].iloc[i]))
            in_run = False
    if in_run:
        low_ends.append(float(filtered["start"].iloc[-1]))
    s1 = ColumnDataSource(data=dict(x=x, y=y_ratio))
    s2 = ColumnDataSource(data=dict(x=x, input=y_input, reconstruction=y_recon))

    for left, right in zip(low_starts, low_ends):
        for p in [p1, p2]:
            p.add_layout(BoxAnnotation(left=left, right=right,
                                       fill_color="red", fill_alpha=0.10, line_color=None))

    p1.line('x', 'y', source=s1, line_width=1, color="#aaaaaa", alpha=0.7)
    p1.y_range.start = -0.1
    p1.y_range.end = 5.1
    p1.xaxis.formatter = NumeralTickFormatter(format="0,0")

    if len(segs) > 0:
        seg_src = ColumnDataSource(dict(
            x0=segs["x0"].tolist(),
            x1=segs["x1"].tolist(),
            y0=segs["cn"].tolist(),
            y1=segs["cn"].tolist(),
            cn=segs["cn"].tolist(),
            color=segs["color"].tolist(),
            label_x=segs["label_x"].tolist(),
            label=segs["label"].tolist(),
            conf_color=segs["conf_color"].tolist(),
        ))
        seg_line = p1.segment("x0", "y0", "x1", "y1", source=seg_src,
                              line_color="color", line_width=3, line_alpha=0.95)
        p1.add_layout(LabelSet(
            x="label_x", y="y0", text="label", source=seg_src,
            text_color="conf_color", text_font_size="9px",
            text_align="center", y_offset=4,
        ))
        p1.add_tools(HoverTool(renderers=[seg_line], tooltips=[("CN", "@cn")]))

    p2.line('x', 'input',          source=s2, line_width=2, color="green", legend_label="input")
    p2.line('x', 'reconstruction', source=s2, line_width=2, color="red",
            legend_label="reconstruction")
    p2.xaxis.formatter = NumeralTickFormatter(format="0,0")

    return column([p1, p2], sizing_mode="stretch_width")