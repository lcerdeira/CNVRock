"""
Validation & robustness page.

Renders the post-scaling-study analyses that test whether the copy-number
*values* are correct, rather than whether events are merely detected:

    long-read depth agreement, per-gene threshold calibration, per-call
    depth-bootstrap confidence intervals, multi-reference chromosomal
    calling for K. variicola, and the calendar-time amplification trend.

Every panel is a pre-rendered PNG committed under ``assets/validation/`` so
this page works on Streamlit Cloud, where the gitignored ``data/results/``
tree is absent. Each figure load is guarded: a missing asset degrades to an
inline note instead of taking the whole app down.
"""
import os

import streamlit as st

_HERE = os.path.dirname(os.path.abspath(__file__))          # …/diagnostics/src
_REPO = os.path.dirname(os.path.dirname(_HERE))             # …/CNVRock
ASSETS = os.path.join(_REPO, "assets", "validation")


def _figure(filename: str, caption: str) -> None:
    """Render a validation figure, or a note if the asset is missing."""
    path = os.path.join(ASSETS, filename)
    if os.path.isfile(path):
        st.image(path, caption=caption, use_container_width=True)
    else:
        st.info(f"Figure not bundled in this deployment (`{filename}`).")


def page_validation():
    st.title("Validation & robustness")
    st.markdown(
        "The scaling study shows CNVRock **detects** copy-number events "
        "reliably. These analyses ask whether the copy-number **values** are "
        "correct, how sensitive the calls are to our choices, and whether the "
        "stated limitations hold up."
    )

    tab_lr, tab_thr, tab_ci, tab_kv, tab_time = st.tabs(
        ["Long-read agreement", "Threshold calibration", "Per-call uncertainty",
         "K. variicola reference", "Temporal trend"])

    # ── orthogonal validation against long-read depth ───────────────────────
    with tab_lr:
        st.subheader("Orthogonal validation against long-read depth")
        c1, c2, c3 = st.columns(3)
        c1.metric("Paired isolates", "255")
        c2.metric("Pearson r", "0.960")
        c3.metric("Amplified recovered", "9 / 10")
        _figure("longread_depth_validation.png",
                "Short-read CRR vs long-read-depth CRR, chromosomal blaSHV "
                "(n = 255).")
        st.markdown(
            "Copy-ratio is computed from **long-read depth, not a long-read "
            "assembly** — assemblers collapse tandem arrays and would erase "
            "the signal under test. The comparison uses a different "
            "chemistry, read-length regime and alignment path, which is what "
            "makes it a genuine external check."
        )
        st.caption(
            "The Pearson correlation is carried largely by the amplified "
            "isolates, which span the widest range; among the mostly "
            "single-copy majority the rank correlation is much weaker. The "
            "operationally meaningful claim is the 9/10 amplified-call "
            "recovery, not the coefficient alone."
        )

    # ── detection robustness to threshold choice ────────────────────────────
    with tab_thr:
        st.subheader("Detection robustness to threshold choice")
        c1, c2 = st.columns(2)
        c1.metric("Mean ROC AUC", "0.976")
        c2.metric("Gene families", "6")
        _figure("threshold_sensitivity.png",
                "MCC as a function of the PCN presence threshold, per family, "
                "with per-family ROC AUC.")
        st.markdown(
            "| Family | *n* positive | ROC AUC | MCC at operating | Best threshold |\n"
            "|---|---|---|---|---|\n"
            "| blaKPC | 865 | 1.000 | 0.985 | 0.10 |\n"
            "| qnrB | 1 527 | 0.999 | 0.974 | 0.30 |\n"
            "| blaNDM | 863 | 0.985 | 0.765 | 0.60 |\n"
            "| blaOXA-48-like | 713 | 0.985 | 0.907 | 0.60 |\n"
            "| aac6-Ib-cr | 1 729 | 0.973 | 0.857 | 0.15 |\n"
            "| blaCTX-M | 2 724 | 0.915 | 0.606 | 0.45 |\n"
        )
        st.markdown(
            "The MCC curves are flat around the operating threshold, so "
            "detection rests on genuine separability rather than a finely "
            "tuned cut-off. blaCTX-M is weakest, for the documented "
            "cross-mapping reason."
        )

    # ── per-call uncertainty ────────────────────────────────────────────────
    with tab_ci:
        st.subheader("Per-call uncertainty")
        c1, c2, c3 = st.columns(3)
        c1.metric("Isolates in scope", "6 078")
        c2.metric("Amplified calls", "162")
        c3.metric("CI excludes 1.0", "162 / 162")
        _figure("percall_uncertainty.png",
                "Top blaSHV amplification calls with 95 % depth-bootstrap "
                "confidence intervals.")
        st.markdown(
            "Each chromosomal blaSHV copy-ratio carries a 95 % interval "
            "obtained by bootstrapping **the estimator the caller actually "
            "reports**: the VAE-normalised gene-bin copy-ratio divided by its "
            "chromosomal flank mean. Poisson resampling of the gene-bin "
            "counts is exact (1 000 replicates); the flank mean, spanning "
            "5 133 bins, is propagated analytically under the central limit "
            "theorem. Every amplification call is statistically separated "
            "from single-copy, so none is a depth-sampling artefact. Median "
            "CI width is 0.168 across all calls and 0.313 among amplified "
            "ones."
        )
        st.caption(
            "Scope matches the main text: chromosomal blaSHV is called only "
            "for K. pneumoniae and K. quasipneumoniae, since K. variicola "
            "carries blaLEN at the syntenic locus and cross-maps. The script "
            "reproduces the caller's own copy-ratio to floating-point "
            "precision (max discrepancy 1.8e-15)."
        )

    # ── multi-reference chromosomal calling ─────────────────────────────────
    with tab_kv:
        st.subheader("Multi-reference chromosomal calling (K. variicola)")
        c1, c2, c3 = st.columns(3)
        c1.metric("Isolates re-aligned", "223")
        c2.metric("blaLEN CRR median", "0.94")
        c3.metric("Spearman ρ vs blaSHV", "0.45")
        _figure("kvariicola_multiref.png",
                "blaLEN copy-ratio on a K. variicola reference (A) vs the "
                "blaSHV value the same isolates receive on HS11286 (B).")
        st.markdown(
            "Chromosomal blaSHV is called only for K. pneumoniae and "
            "K. quasipneumoniae, because K. variicola carries the LEN-family "
            "homolog at the syntenic locus and its reads cross-map. "
            "Re-aligning the 223 K. variicola isolates to NC_011283.1 makes "
            "blaLEN callable and centred on single-copy — so the restriction "
            "is a **reference-choice limitation, not an architectural one**. "
            "The ρ = 0.45 also shows the HS11286 value is a distorted proxy "
            "rather than pure noise, which is why these isolates are excluded "
            "rather than reinterpreted."
        )

    # ── temporal trend ──────────────────────────────────────────────────────
    with tab_time:
        st.subheader("Temporal trend in amplification")
        c1, c2 = st.columns(2)
        c1.metric("Isolates with a usable year", "7 448")
        c2.metric("Genes with a significant trend", "0 / 3")
        _figure("temporal_amplification.png",
                "Amplification prevalence among carriers by collection year "
                "(2013-2024).")
        st.markdown(
            "| Gene | Carriers | Amplified | OR per year | 95 % CI |\n"
            "|---|---|---|---|---|\n"
            "| blaKPC | 539 | 123 | 0.947 | 0.881-1.008 |\n"
            "| blaCTX-M | 2 401 | 420 | 0.983 | 0.946-1.018 |\n"
            "| blaSHV (chr) | 3 731 | 98 | 0.983 | 0.946-1.024 |\n"
        )
        st.markdown(
            "No gene shows a significant calendar-time trend, supporting the "
            "steady-state assumption behind the per-gene λ/μ estimates."
        )
        st.warning(
            "Unadjusted, first-pass trend. Which clones and countries were "
            "sequenced shifts across years, so a trend — or its absence — may "
            "reflect ascertainment as much as biology. A confounder-adjusted "
            "hierarchical model remains future work."
        )
