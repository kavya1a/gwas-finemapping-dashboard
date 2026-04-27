"""Background distribution calibration for variant scores.

## How AlphaGenome quantile scores are calibrated

The AlphaGenome API returns a ``quantile_score`` alongside each ``raw_score``.
Per the SDK documentation and the AlphaGenome Nature 2026 paper (Extended Data),
quantile scores are pre-computed by ranking the query variant's raw score
within a background distribution of approximately 300,000 common variants
(minor allele frequency > 0.01) drawn genome-wide from gnomAD / 1000 Genomes.

Key properties of this calibration:
  - Scale: signed [-1, 1] for directional scorers (DIFF_*);
           unsigned [0, 1] for magnitude-only scorers (L2_DIFF*, contact maps).
  - Saturation: scores saturate at ±0.9999 due to the finite (~300K) reference
    set. Multiple variants at the ceiling cannot be ranked by quantile alone;
    use raw_score as a tiebreaker.
  - NOT allele-frequency-matched: the background pool is all common variants
    (MAF > 0.01), not stratified by the query variant's own allele frequency.

## Implication for rare variants

Because the background is genome-wide common variants, a rare variant
(MAF < 0.01) compared against this background will tend to show inflated
quantile scores if rare variants genuinely have larger effect sizes (as
expected under purifying selection). In practice, this affects variants like
TREM2 rs75932628 (MAF ~0.003) and LRRK2 G2019S rs34637584 (MAF < 0.001).

We flag variants with MAF < 0.01 in the output with ``rare_variant_caution=True``
so downstream users know the quantile score may overstate relative effect.

## Effort to implement allele-frequency-matched calibration

High-effort (~3–4 days):
  1. Collect a large set of background variants (e.g., 10,000) stratified into
     MAF bins (e.g., [0.001–0.01], [0.01–0.05], [0.05–0.2], [0.2–0.5]).
  2. Score all background variants via AlphaGenome API (~10,000 calls × $cost).
  3. For each query variant, find its MAF bin and compute within-bin percentile.
  4. Return bin-calibrated quantile score.
The SDK does not expose the raw background distribution, so this must be built
entirely client-side. The benefit is modest for common variants (MAF > 0.05)
but meaningful for rare/low-frequency variants.

Current workaround: use ``rare_variant_caution`` flag and downweight quantile
scores for rare variants in the composite (see ``RARE_VARIANT_DISCOUNT`` in
composite.py).

## Two scoring modes

  1. API quantile scores (default, free): Use the SDK-provided quantile_score
     directly. These are the authoritative per-track calibration.

  2. Local null sampling (optional, expensive): Score N synthetic nearby
     variants at non-functional positions to build a local background. Returns
     z-scores relative to that distribution. Useful for controlling local
     sequence composition biases but requires N extra API calls per locus.
     Disabled by default (n_null=0).

The primary output of this module is a per-modality summary with:
  - max_abs_quantile: max |quantile_score| across tissue-filtered tracks
  - mean_abs_quantile: mean |quantile_score| across tissue-filtered tracks
  - frac_above_threshold: fraction of tracks with |quantile_score| > threshold
  - top_track: biosample_name of the single highest-scoring track
  - raw_score_max: raw score of the highest-scoring track (for sign)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Fraction of tracks (within a modality/tissue group) that must exceed
# this quantile threshold to count towards the credibility bonus.
QUANTILE_THRESHOLD: float = 0.90

# Output type strings that belong to each modality group, matching the
# `output_type` column of tidy_scores().
OUTPUT_TYPE_TO_MODALITY: dict[str, str] = {
    "RNA_SEQ": "expression",
    "CAGE": "expression",
    "PROCAP": "expression",
    "ATAC": "chromatin",
    "DNASE": "chromatin",
    "CHIP_HISTONE": "chromatin",
    "CHIP_TF": "tf_binding",
    "SPLICE_SITES": "splice_sites",        # kept separate, combined later
    "SPLICE_SITE_USAGE": "splice_site_usage",
    "SPLICE_JUNCTIONS": "splice_junctions",
    "CONTACT_MAPS": "contact",
    "POLYADENYLATION": "polyadenylation",
}


def _score_col(df: pd.DataFrame) -> pd.Series:
    """Returns quantile_score if present, else falls back to raw_score."""
    if "quantile_score" in df.columns and df["quantile_score"].notna().any():
        return df["quantile_score"].fillna(0.0)
    # Normalize raw_score to [-1, 1] via tanh as a rough quantile proxy
    raw = df["raw_score"].fillna(0.0)
    return pd.Series(np.tanh(raw.values / (raw.abs().quantile(0.95) + 1e-9)), index=df.index)


def summarize_modality_scores(
    tidy_df: pd.DataFrame,
    quantile_threshold: float = QUANTILE_THRESHOLD,
) -> pd.DataFrame:
    """Compute per-modality summary stats for a single variant's tidy DataFrame.

    Args:
        tidy_df: Output of tidy_scores() for one variant, optionally tissue-
            filtered via tissue_config.filter_tracks().
        quantile_threshold: Fraction threshold for frac_above_threshold.

    Returns:
        DataFrame with one row per modality group, columns:
            modality, max_abs_score, mean_abs_score, frac_above_threshold,
            signed_max_score, top_track_biosample, n_tracks.
    """
    if tidy_df is None or tidy_df.empty:
        return pd.DataFrame(columns=[
            "modality", "max_abs_score", "mean_abs_score",
            "frac_above_threshold", "signed_max_score",
            "top_track_biosample", "n_tracks",
        ])

    scores_col = _score_col(tidy_df)
    tidy_df = tidy_df.copy()
    tidy_df["_score"] = scores_col

    # Map output_type to modality group
    tidy_df["_modality"] = (
        tidy_df["output_type"]
        .astype(str)
        .str.replace("OutputType.", "", regex=False)
        .map(OUTPUT_TYPE_TO_MODALITY)
    )

    rows = []
    for modality, group_df in tidy_df.groupby("_modality", dropna=True):
        s = group_df["_score"]
        abs_s = s.abs()
        idx_max = abs_s.idxmax()
        signed_max = float(s.loc[idx_max])

        top_bio = (
            group_df.loc[idx_max, "biosample_name"]
            if "biosample_name" in group_df.columns
            else None
        )

        rows.append({
            "modality": modality,
            "max_abs_score": float(abs_s.max()),
            "mean_abs_score": float(abs_s.mean()),
            "frac_above_threshold": float((abs_s > quantile_threshold).mean()),
            "signed_max_score": signed_max,
            "top_track_biosample": top_bio,
            "n_tracks": len(group_df),
        })

    return pd.DataFrame(rows)
