import math
import streamlit as st
import plotly.express as px
import pandas as pd

from gwas_catalog import fetch_gwas_variants
from alphagenome_scorer import score_variants

DISEASES = {
    "Alzheimer's Disease": "alzheimers",
    "Type 2 Diabetes": "t2d",
    "Schizophrenia": "schizophrenia",
    "Parkinson's Disease": "parkinsons",
}

st.set_page_config(page_title="GWAS Fine-Mapping Dashboard", layout="wide")
st.title("GWAS Fine-Mapping Dashboard")
st.caption("Powered by NHGRI-EBI GWAS Catalog + AlphaGenome")

with st.sidebar:
    disease_display = st.selectbox("Disease", list(DISEASES.keys()))
    disease_slug = DISEASES[disease_display]
    max_variants = st.slider("Max variants to fetch", 10, 200, 50)
    run = st.button("Run analysis", type="primary")

if run:
    with st.spinner(f"Fetching GWAS variants for {disease_display}…"):
        try:
            variants = fetch_gwas_variants(disease_slug, max_variants=max_variants)
        except Exception as e:
            st.error(f"GWAS Catalog error: {e}")
            st.stop()

    if not variants:
        st.warning("No variants returned from GWAS Catalog for this query.")
        st.stop()

    n_with_alleles = sum(1 for v in variants if v.get("ref") and v.get("alt"))
    st.info(
        f"Fetched {len(variants)} variants ({n_with_alleles} with allele info). "
        "Scoring with AlphaGenome…"
    )

    if n_with_alleles == 0:
        st.warning(
            "No variants have ref/alt allele information — the GWAS Catalog v2 API "
            "does not always return alleles. Try a smaller batch or check the trait name."
        )
        st.stop()

    with st.spinner("Scoring variants (this takes ~10–30 s per variant)…"):
        try:
            df = score_variants(variants, disease=disease_slug)
        except EnvironmentError as e:
            st.error(str(e))
            st.stop()
        except Exception as e:
            st.error(f"AlphaGenome scoring error: {e}")
            st.stop()

    if df.empty:
        st.warning("No variants could be scored.")
        st.stop()

    # ── Main results table ────────────────────────────────────────────────────
    st.subheader("Ranked variants by causal impact")

    display_cols = [
        "rank", "rsid", "chrom", "pos", "ref", "alt",
        "composite_score", "pip_weighted_score", "pip",
        "p_value", "rare_variant_caution", "error",
    ]
    show_df = df[[c for c in display_cols if c in df.columns]].copy()

    # Format scores to 4 decimal places for readability
    for col in ("composite_score", "pip_weighted_score"):
        if col in show_df.columns:
            show_df[col] = show_df[col].apply(
                lambda x: f"{x:.4f}" if pd.notna(x) else ""
            )

    st.dataframe(show_df, use_container_width=True, hide_index=True)

    # ── Scatter: p-value vs composite score ──────────────────────────────────
    plot_df = df.copy()
    if "p_value" in plot_df.columns and plot_df["p_value"].notna().any():
        plot_df["neg_log10_p"] = plot_df["p_value"].apply(
            lambda p: -math.log10(p) if p and p > 0 else None
        )
        fig = px.scatter(
            plot_df.dropna(subset=["neg_log10_p", "composite_score"]),
            x="neg_log10_p",
            y="composite_score",
            hover_data=["rsid", "chrom", "pos", "ref", "alt"],
            color="rare_variant_caution" if "rare_variant_caution" in plot_df.columns else None,
            labels={
                "neg_log10_p": "−log₁₀(p-value)",
                "composite_score": "AlphaGenome Composite Score",
                "rare_variant_caution": "Rare variant",
            },
            title=f"GWAS p-value vs. AlphaGenome Composite Score — {disease_display}",
        )
        st.plotly_chart(fig, use_container_width=True)

    # ── Per-variant modality breakdown ────────────────────────────────────────
    with st.expander("Per-variant modality breakdown"):
        modality_breakdowns = getattr(df, "modality_breakdowns", {})
        for _, row in df.iterrows():
            rsid = row.get("rsid", "")
            composite = row.get("composite_score")
            score_str = f"{composite:.4f}" if pd.notna(composite) else "N/A"
            rare_flag = " ⚠ rare variant" if row.get("rare_variant_caution") else ""
            err = row.get("error")

            st.markdown(
                f"**{rsid}** — chr{row.get('chrom')}:{row.get('pos')} "
                f"{row.get('ref', '?')}→{row.get('alt', '?')}  "
                f"| composite: `{score_str}`{rare_flag}"
            )

            if err:
                st.caption(f"  Error: {err}")
            elif rsid in modality_breakdowns and not modality_breakdowns[rsid].empty:
                breakdown = modality_breakdowns[rsid][
                    ["modality", "max_abs_score", "frac_above_threshold"]
                ].copy()
                breakdown["max_abs_score"] = breakdown["max_abs_score"].apply(
                    lambda x: f"{x:.4f}"
                )
                breakdown["frac_above_threshold"] = breakdown["frac_above_threshold"].apply(
                    lambda x: f"{x:.2%}"
                )
                st.dataframe(breakdown, use_container_width=True, hide_index=True)
