import json
import os

import pandas as pd
import streamlit as st


RESULTS_ROOT = "../data/results"


def _read_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return None


@st.fragment(run_every=3)
def _training_monitor(log_path):
    data = _read_json(log_path)
    if data is None:
        st.caption("No training log yet.")
        return

    status = data.get("status", "unknown")
    epoch  = data.get("epoch", 0)
    total  = data.get("total_epochs", 1)

    col_stat, col_bar = st.columns([1, 4])
    with col_stat:
        if status == "done":
            st.success(f"Done  ({epoch} epochs)")
        else:
            st.info(f"Epoch {epoch} / {total}")
    with col_bar:
        st.progress(epoch / max(total, 1))

    history = data.get("history", [])
    if history:
        df = pd.DataFrame(history).set_index("epoch")
        st.line_chart(df[["recon", "kl"]], height=260)
        # Beta ramp on a separate (much smaller) scale
        with st.expander("Beta schedule", expanded=False):
            st.line_chart(df[["beta"]], height=160)


def page_monitor():
    st.title("Monitor")

    dirs = sorted(os.listdir(RESULTS_ROOT)) if os.path.isdir(RESULTS_ROOT) else []
    if not dirs:
        st.warning(f"No results directories found in {RESULTS_ROOT}")
        return

    results_dir = st.selectbox("Results directory", dirs)
    out_dir     = os.path.join(RESULTS_ROOT, results_dir)

    st.subheader("Training")
    _training_monitor(os.path.join(out_dir, "training_log.json"))
