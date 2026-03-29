import streamlit as st

from src.page1 import page1
from src.page2 import page2

st.set_page_config(layout="wide")

pg = st.navigation([
    st.Page(page1, title="First page"),
    st.Page(page2, title="Latent coverage"),
])

pg.run()