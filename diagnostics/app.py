import streamlit as st

from src.page1 import page1

pg = st.navigation([
    st.Page(page1, title="First page")
])

pg.run()