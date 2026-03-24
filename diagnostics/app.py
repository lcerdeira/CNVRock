import streamlit as st

from src.page1 import page1

st.set_page_config(layout="wide")

pg = st.navigation([
    st.Page(page1, title="First page")
])

pg.run()