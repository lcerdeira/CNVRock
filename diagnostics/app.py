import streamlit as st

from src.page1 import page1
from src.page2 import page2
from src.page_monitor import page_monitor

st.set_page_config(layout="wide")

pg = st.navigation([
    st.Page(page_monitor, title="Monitor",          icon="📡"),
    st.Page(page1,        title="Sample viewer",    icon="🔬"),
    st.Page(page2,        title="Latent coverage",  icon="🗺️"),
])

pg.run()