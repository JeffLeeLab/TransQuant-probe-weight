"""TransQuant probe weight factor (W) — Streamlit / stlite app."""

import streamlit as st

import transquant_w.ui

st.set_page_config(page_title="TransQuant probe weight W", page_icon="🧬", layout="centered")

transquant_w.ui.render()
