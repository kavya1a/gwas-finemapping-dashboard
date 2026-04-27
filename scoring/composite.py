"""Composite causal-likelihood scorer.

Pipeline per variant:
  1. Call AlphaGenome score_variant() → list[AnnData]
  2. tidy_scores() → long DataFrame
  3. tissue_config.filter_tracks() → tissue-relevant subset
  4. null_distribution.summarize_modality_scores() → per-modality stats
  5. Weighted sum across modality groups → composite_score in [0, 1]
  6. Optional PIP multiplication → pip_weighted_score

Modality weights are defined in MODALITY_WEIGHTS. The splice sub-scores
(splice_sites, splice_site_usage, splice_junctions) are averaged first into
a single "splicing" score before weighting, so the group weight applies once.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers as vs

from scoring.diff_functions import get_scorers_for_organism
from scoring.null_distribution import summarize_modality_scores
from scoring.tissue_config import TissueProfile, filter_tracks, get_profile

load_dotenv(Path(__file__).parent.parent / ".env")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Weights sum to 1.0. Splice sub-scorers are averaged first, so this
# "splicing" weight applies to the combined splice score.
MODALITY_WEIGHTS: dict[str, float] = {
    "expression": 0.30,   # RNA-seq + CAGE + PRO-cap
    "chromatin": 0.25,    # ATAC + DNase + ChIP-Histone
    "tf_binding": 0.15,   # ChIP-TF
    "splicing": 0.20,     # splice_sites + splice_site_usage + splice_junctions (averaged)
    "contact": 0.05,      # contact maps
    "polyadenylation": 0.05,  # PAS usage
}

# Splice sub-modalities that are averaged into one "splicing" composite
SPLICE_SUBMODALITIES = {"splice_sites", "splice_site_usage", "splice_junctions"}

# Interval half-width for AlphaGenome (must be exactly 524288 bp = 1 Mb / 2)
INTERVAL_HALF: int = 524288

# AlphaGenome sequence length options: 16384, 98304, 524288, 1048576
SEQ_LEN: int = 1_048_576


# Quantile score discount applied to rare variants (MAF < 0.01) to account
# for the fact that the SDK's background distribution uses only common variants.
# A value of 0.8 means rare variant scores are reported at 80% face value.
RARE_VARIANT_DISCOUNT: float = 0.80
RARE_VARIANT_MAF_THRESHOLD: float = 0.01


@dataclass
class VariantInput:
    rsid: str
    chrom: str       # e.g. "chr19"
    pos: int         # 1-based VCF position
    ref: str
    alt: str
    pip: float | None = None   # Posterior Inclusion Probability (optional)
    p_value: float | None = None
    maf: float | None = None   # Minor allele frequency (for rare-variant flag)


def _build_interval(chrom: str, pos_1based: int) -> genome.Interval:
    center = pos_1based - 1  # 0-based
    return genome.Interval(
        chromosome=chrom,
        start=max(0, center - INTERVAL_HALF),
        end=center + INTERVAL_HALF,
    )


def _compute_composite(modality_df: pd.DataFrame) -> float:
    """Weighted sum of per-modality max_abs_scores → composite in [0, 1]."""
    if modality_df.empty:
        return 0.0

    score_by_modality: dict[str, float] = {}
    for _, row in modality_df.iterrows():
        score_by_modality[row["modality"]] = row["max_abs_score"]

    # Average splice sub-scores into single "splicing" value
    splice_vals = [
        score_by_modality[m]
        for m in SPLICE_SUBMODALITIES
        if m in score_by_modality
    ]
    if splice_vals:
        score_by_modality["splicing"] = sum(splice_vals) / len(splice_vals)
    for m in SPLICE_SUBMODALITIES:
        score_by_modality.pop(m, None)

    composite = 0.0
    weight_used = 0.0
    for modality, weight in MODALITY_WEIGHTS.items():
        val = score_by_modality.get(modality, 0.0)
        # Cap at 1.0 — quantile scores saturate at ~0.9999
        composite += weight * min(val, 1.0)
        weight_used += weight

    return composite / weight_used if weight_used > 0 else 0.0


def _process_single_result(
    variant: VariantInput,
    scores: list,
    tissue_profile: TissueProfile,
) -> dict:
    """Convert raw SDK scores for one variant into a result dict."""
    tidy_df = vs.tidy_scores(scores)
    if tidy_df is None or tidy_df.empty:
        return _error_result(variant, "empty scores")

    filtered_df = filter_tracks(tidy_df, tissue_profile)
    modality_df = summarize_modality_scores(filtered_df)
    composite = _compute_composite(modality_df)

    # Discount rare variants — their quantile scores are calibrated against a
    # common-variant background, potentially inflating the score.
    rare = variant.maf is not None and variant.maf < RARE_VARIANT_MAF_THRESHOLD
    if rare:
        composite *= RARE_VARIANT_DISCOUNT

    pip_score = composite * variant.pip if variant.pip is not None else composite

    return {
        "rsid": variant.rsid,
        "chrom": variant.chrom,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "maf": variant.maf,
        "composite_score": composite,
        "pip_weighted_score": pip_score,
        "pip": variant.pip,
        "rare_variant_caution": rare,
        "error": None,
        "modality_breakdown": modality_df,
        "tidy_df": tidy_df,
    }


def _error_result(variant: VariantInput, msg: str) -> dict:
    return {
        "rsid": variant.rsid,
        "chrom": variant.chrom,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "maf": variant.maf,
        "composite_score": None,
        "pip_weighted_score": None,
        "pip": variant.pip,
        "rare_variant_caution": False,
        "error": msg,
        "modality_breakdown": pd.DataFrame(),
        "tidy_df": pd.DataFrame(),
    }


def score_single_variant(
    model: dna_client.DnaClient,
    variant: VariantInput,
    tissue_profile: TissueProfile,
) -> dict:
    """Score one variant; returns a result dict."""
    interval = _build_interval(variant.chrom, variant.pos)
    sdk_variant = genome.Variant(
        chromosome=variant.chrom,
        position=variant.pos,
        reference_bases=variant.ref,
        alternate_bases=variant.alt,
    )
    try:
        raw_scores = model.score_variant(
            interval, sdk_variant, variant_scorers=get_scorers_for_organism()
        )
    except Exception as e:
        return _error_result(variant, str(e))

    return _process_single_result(variant, raw_scores, tissue_profile)


def score_variants_composite(
    variants: list[VariantInput],
    disease: str = "",
    tissue_profile: TissueProfile | None = None,
    api_key: str | None = None,
    max_workers: int = 4,
) -> pd.DataFrame:
    """Score a list of variants using the SDK's parallel batch endpoint.

    Uses model.score_variants() (plural) which runs up to max_workers requests
    concurrently. Significantly faster than the sequential single-variant path
    for batches of 5+ variants.

    Args:
        variants: List of VariantInput objects.
        disease: Disease slug (e.g. "alzheimers") used to auto-select tissue
            profile when tissue_profile is None.
        tissue_profile: Explicit tissue profile; overrides disease lookup.
        api_key: AlphaGenome API key; reads from ALPHAGENOME_API_KEY env if None.
        max_workers: Parallel requests to AlphaGenome API.

    Returns:
        DataFrame sorted by pip_weighted_score (if PIPs provided) or
        composite_score, with modality_breakdowns and tidy_dfs dicts attached.
    """
    key = api_key or os.environ.get("ALPHAGENOME_API_KEY", "")
    if not key:
        raise EnvironmentError("ALPHAGENOME_API_KEY not set")

    profile = tissue_profile or get_profile(disease)
    model = dna_client.create(key)
    scorers = get_scorers_for_organism()

    intervals = [_build_interval(v.chrom, v.pos) for v in variants]
    sdk_variants = [
        genome.Variant(
            chromosome=v.chrom,
            position=v.pos,
            reference_bases=v.ref,
            alternate_bases=v.alt,
        )
        for v in variants
    ]

    print(f"  Batch-scoring {len(variants)} variants (max_workers={max_workers})...")
    try:
        all_scores = model.score_variants(
            intervals,
            sdk_variants,
            variant_scorers=scorers,
            max_workers=max_workers,
            progress_bar=True,
        )
    except Exception as e:
        # Fall back to sequential scoring if batch fails
        print(f"  Batch failed ({e}), falling back to sequential...")
        all_scores = []
        for i, (iv, sv) in enumerate(zip(intervals, sdk_variants)):
            print(f"    [{i+1}/{len(variants)}] {variants[i].rsid}...")
            try:
                all_scores.append(model.score_variant(iv, sv, variant_scorers=scorers))
            except Exception as e2:
                all_scores.append(None)
                print(f"    ERROR: {e2}")

    results = []
    for v, scores in zip(variants, all_scores):
        if scores is None:
            results.append(_error_result(v, "scoring failed"))
        else:
            res = _process_single_result(v, scores, profile)
            flag = " [RARE]" if res["rare_variant_caution"] else ""
            err = f" ERROR: {res['error']}" if res["error"] else f" composite={res['composite_score']:.4f}"
            print(f"  {v.rsid}{flag}:{err}")
            results.append(res)

    sort_col = (
        "pip_weighted_score"
        if any(v.pip is not None for v in variants)
        else "composite_score"
    )
    rows = [{k: v for k, v in r.items() if k not in ("modality_breakdown", "tidy_df")}
            for r in results]

    df = pd.DataFrame(rows)
    df = df.sort_values(sort_col, ascending=False, na_position="last").reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    df._metadata = ["modality_breakdowns", "tidy_dfs"]
    df.modality_breakdowns = {r["rsid"]: r["modality_breakdown"] for r in results}
    df.tidy_dfs = {r["rsid"]: r["tidy_df"] for r in results}

    return df
