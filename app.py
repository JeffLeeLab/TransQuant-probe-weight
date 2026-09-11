"""TransQuant probe weight factor (W) — Streamlit / stlite app."""

import streamlit as st

from transquant_w import core

st.set_page_config(page_title="TransQuant probe weight W", page_icon="🧬", layout="centered")

st.title("TransQuant: Compute Probe Weight Factor")
st.markdown(
    "Computes the probe weight factor **W**, gene length **L** and probe localisation "
    "profile **N** for an smFISH probe set, as defined in "
    "[Halpern & Itzkovitz, *Methods* 2016](https://doi.org/10.1016/j.ymeth.2015.11.015) "
)

# ------------------------------------------------------------------ inputs --

st.subheader("1 · Target gene sequence")
# st.info(
#     "**Paste the full genomic locus, in the direction of transcription.**\n\n"
#     "- Include *all* exons, introns and UTRs, from transcription start to end. "
#     "Do **not** use the spliced mRNA/cDNA: intron lengths set how long Pol II spends between probes.\n"
#     "- The sequence must read 5' → 3' along the RNA (the coding/sense strand). "
#     "For a minus-strand gene in the UCSC Genome Browser, tick *Reverse complement* when downloading DNA. "
#     "The app does not read strand information from FASTA headers.\n"
#     "- FASTA or plain sequence; `N` bases are fine. A multi-entry FASTA (UCSC exon/intron download) "
#     "is concatenated in order."
# )
st.caption(
    "**Paste the full genomic locus, in the direction of transcription.**\n"
    "- Include *all* exons, introns and UTRs, from transcription start to end. "
    "Do **not** use the spliced mRNA/cDNA: intron lengths set how long Pol II spends between probes.\n"
    "- The sequence must read 5' → 3' along the RNA (the coding/sense strand).\n"
    "- FASTA or plain sequence; `N` bases are fine. A multi-entry FASTA "
    "is concatenated in order."
)
tcol1, tcol2 = st.columns([3, 2])
with tcol1:
    target_text = st.text_area("Target sequence (paste)", height=180, placeholder=">gene chr:start-end\nACGT...",
                               label_visibility="collapsed")
with tcol2:
    target_file = st.file_uploader("…or upload FASTA", type=["fa", "fasta", "txt"], key="target_file")
    concat = st.toggle("Concatenate multi-entry FASTA", value=True)

st.subheader("2 · Probe sequences")
st.caption(
    "One probe per line, or paste rows straight from Excel / Google Sheets / the Stellaris designer. "
    "Extra columns and a header row are detected and ignored. Probes are the antisense oligos as ordered."
)
pcol1, pcol2 = st.columns([3, 2])
with pcol1:
    probes_text = st.text_area("Probe sequences (paste)", height=180, placeholder="1\tACGT...\t50\n2\t...",
                               label_visibility="collapsed")
with pcol2:
    probes_file = st.file_uploader("…or upload tab-delimited .txt", type=["txt", "tsv", "csv"], key="probes_file")
    revcomp = st.toggle("Reverse complement probes", value=True,
                        help="On: probes are antisense oligos (normal). Off: you pasted target-sense sequences.")

run = st.button("Compute W", type="primary", use_container_width=True)
st.markdown("<div style='height:1.5rem'></div>", unsafe_allow_html=True)

# ----------------------------------------------------------------- compute --

if run:
    t_raw = target_file.getvalue().decode("utf-8", "replace") if target_file else target_text
    p_raw = probes_file.getvalue().decode("utf-8", "replace") if probes_file else probes_text
    try:
        if not t_raw.strip():
            raise core.InputError("Paste or upload a target sequence.")
        if not p_raw.strip():
            raise core.InputError("Paste or upload probe sequences.")
        target, w1 = core.parse_target(t_raw, concatenate=concat)
        probes, w2 = core.parse_probes(p_raw)
        res = core.compute(target, probes, revcomp_probes=revcomp, warnings=w1 + w2)
    except core.InputError as e:
        st.error(str(e))
        st.stop()

    st.divider()
    st.subheader("Results")
    for w in res.warnings:
        st.warning(w)

    m1, m2, m3 = st.columns(3)
    m1.metric("Probe weight factor W", f"{res.W:.4f}")
    m2.metric("Gene length L (bp)", f"{res.L:,}")
    m3.metric("Probes matched / supplied", f"{res.n_matched_probes} / {res.n_probes}",
              delta=None if res.n_sites == res.n_matched_probes else f"{res.n_sites} binding sites",
              delta_color="off")

    fig = core.make_plot(res)
    st.pyplot(fig, use_container_width=True)
    d1, d2 = st.columns(2)
    d1.download_button("Download plot (PNG)", core.figure_bytes(fig, "png"),
                       file_name=f"probe_profile_W{res.W:.3f}.png", mime="image/png", use_container_width=True)
    d2.download_button("Download plot (SVG)", core.figure_bytes(fig, "svg"),
                       file_name=f"probe_profile_W{res.W:.3f}.svg", mime="image/svg+xml", use_container_width=True)

    st.markdown("**Methods sentence**")
    st.code(
        f"The probe weight factor W = {res.W:.3f} was computed for {res.n_matched_probes} probes "
        f"({res.n_sites} binding sites) over the {res.L:,} bp genomic locus using TransQuant-probe-weight, "
        "following Bahar Halpern & Itzkovitz (Methods, 2016).",
        language=None,
    )

    with st.expander("Probe positions", expanded=False):
        rows = [{"probe": h.probe, "start": h.start, "end": h.end, "midpoint": h.mid,
                 "sites": res.multi_site.get(h.probe, 1)} for h in res.hits]
        rows += [{"probe": p, "start": None, "end": None, "midpoint": None, "sites": 0} for p in res.unmatched]
        st.dataframe(rows, use_container_width=True, hide_index=True)

# ------------------------------------------------------------ how to use W --

with st.expander("How to use W, L and your smFISH spot count measurements", expanded=not run):
    st.markdown(
        """<style>
        [data-testid="stExpander"] .katex-display { text-align: left !important; margin: 0.3rem 0 1rem 0; overflow-x: auto; }
        [data-testid="stExpander"] .katex-display > .katex { text-align: left !important; font-size: 0.9rem; }
        [data-testid="stExpander"] .katex-display > .katex > .katex-html { text-align: left !important; }
        /* Streamlit shrink-wraps and centres the st.latex block itself; force it to full width, left edge. */
        [data-testid="stExpander"] div:has(> .katex-display),
        [data-testid="stExpander"] div:has(> div > .katex-display) { width: 100% !important; margin-left: 0 !important; margin-right: 0 !important; align-items: flex-start !important; }
        </style>""",
        unsafe_allow_html=True,
    )
    st.markdown("**1. Calculate RNA *Transcription rate***")
    st.latex(
        r"\text{Transcription rate (mRNA} \cdot \text{hour}^{-1}\text{)} = "
        r"\frac{\dfrac{\text{nascent transcript number}}{\text{probe weight factor } W} \times \text{elongation rate}}"
        r"{\text{gene length } L}"
    )
    st.markdown("**2. Calculate RNA *Decay rate***")
    st.latex(
        r"\text{Decay rate (hour}^{-1}\text{)} = "
        r"\frac{\text{chromosome fraction} \times \text{transcription rate} \times \text{number of chromosome copies}}"
        r"{\text{transcripts in the cell}}"
    )
    st.markdown("**3. Calculate RNA *Half-life***")
    st.latex(r"\text{Half-life (min)} = \frac{\ln(2)}{\text{decay rate}} \times 60")
    st.caption("See Eqs. 3–8 in the original publication for the full treatment.")

st.caption(
    "Source and CLI: github.com/JeffLeeLab/TransQuant-probe-weight · "
    "Algorithm reimplemented from the TransQuant MATLAB code (compute_weight_factor.m)."
)
