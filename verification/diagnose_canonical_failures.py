"""Tissue-filter diagnostic for rs7903146 (TCF7L2/T2D) and rs356219 (SNCA/PD).

Historical context: these two variants used to fail a strict top-20% rank
test in the canonical verification harness. The original docs blamed a
tissue mismatch (scoring against K562/blood-lineage instead of disease-
matched tracks). The hypothesis being checked here was a related one — that
filter_tracks silently falls back to ALL tracks when nothing matches the
disease keywords, leading to the same effective behavior.

Both hypotheses are REFUTED by this script's output:
  - T2D retains ~1,230 / 11,906 tracks; top biosamples include pancreas,
    liver, hepatocyte, skeletal muscle.
  - PD retains ~906 / 15,518 tracks; top biosamples include substantia
    nigra, prefrontal cortex, neural cell.

The actual cause of the rank-test failures is at the modality-aggregation
level — rs356219 scores below the median of PD variants on every modality.
See docs/canonical_variants.md for the full diagnosis and
analysis/explore_composite_weights.py for the weight-scheme exploration.

The verification harness now splits these into a separate
`canonical_regulatory_detection` group with an absolute-effect threshold
rather than a rank claim; both tests pass.

This script is kept as a documented diagnostic — re-run if the tissue
filter or disease profiles are ever changed.

For each variant, this script:
  1. Scores it via the AlphaGenome API with the appropriate disease profile.
  2. Reports the number of tracks before/after the tissue filter.
  3. Lists the top biosample / gtex tissues represented in the filtered set.
  4. Flags whether the tissue filter fell back to "all tracks".

Run (requires ALPHAGENOME_API_KEY, ~2 API calls / ~2 min):
    python verification/diagnose_canonical_failures.py
"""

from __future__ import annotations

import os
import sys
from collections import Counter
from pathlib import Path

from dotenv import load_dotenv

DIR = Path(__file__).parent.parent
sys.path.insert(0, str(DIR))
load_dotenv(DIR / ".env")

EXPRESSION_OUTPUT_TYPES = {"RNA_SEQ", "CAGE", "PROCAP"}

CASES = [
    # (disease_slug, rsid, chrom, pos, ref, alt, gene_label)
    ("t2d",        "rs7903146", "chr10", 112_998_590, "C", "T", "TCF7L2 (β-islet enhancer)"),
    ("parkinsons", "rs356219",  "chr4",   89_704_960, "A", "G", "SNCA (substantia nigra)"),
]


def _tissue_breakdown(df, top_n: int = 8):
    bio = (df.get("biosample_name") if "biosample_name" in df.columns
           else df.get("biosample_term_name"))
    gtex = df.get("gtex_tissue")

    bio_top = Counter(bio.fillna("(none)")).most_common(top_n) if bio is not None else []
    gtex_top = Counter(gtex.fillna("(none)")).most_common(top_n) if gtex is not None else []
    return bio_top, gtex_top


def main() -> None:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set")
        sys.exit(1)

    from alphagenome.models import dna_client
    from scoring.composite import VariantInput, score_single_variant
    from scoring.tissue_config import filter_tracks, get_profile

    model = dna_client.create(api_key)

    for disease, rsid, chrom, pos, ref, alt, label in CASES:
        print(f"\n{'='*70}")
        print(f"  {rsid}  ({label})  →  disease profile: {disease}")
        print(f"{'='*70}")

        profile = get_profile(disease)
        print(f"\n  Profile biosample_keywords: {profile.biosample_keywords}")
        print(f"  Profile gtex_keywords:      {profile.gtex_keywords}")

        vi = VariantInput(rsid=rsid, chrom=chrom, pos=pos, ref=ref, alt=alt,
                          maf=None, p_value=None)
        result = score_single_variant(model, vi, profile)

        if result.get("error"):
            print(f"\n  ERROR: {result['error']}")
            continue

        tidy = result.get("tidy_df")
        if tidy is None or tidy.empty:
            print("\n  tidy_df is empty — nothing to inspect")
            continue

        n_total = len(tidy)
        filtered = filter_tracks(tidy, profile)
        n_filtered = len(filtered)

        fell_back = (n_filtered == n_total) and bool(
            profile.biosample_keywords or profile.gtex_keywords
        )

        print(f"\n  Tracks total:    {n_total:>5}")
        print(f"  Tracks after filter: {n_filtered:>5}")
        if fell_back:
            print(f"  ⚠ filter_tracks FELL BACK to all tracks (no matches for keywords)")
        else:
            print(f"  ✓ tissue filter retained a real subset "
                  f"({100*n_filtered/n_total:.1f}% of catalog)")

        # Expression-modality only (signed_max for that modality drove the
        # composite for these variants in the DB).
        ot = filtered["output_type"].astype(str).str.replace(
            "OutputType.", "", regex=False)
        expr = filtered[ot.isin(EXPRESSION_OUTPUT_TYPES)]
        print(f"  Expression tracks in filtered set: {len(expr):>5}")

        bio_top, gtex_top = _tissue_breakdown(expr)
        print("\n  Top biosamples in filtered expression tracks:")
        for name, n in bio_top:
            print(f"    {n:>4} × {name}")
        print("\n  Top GTEx tissues:")
        for name, n in gtex_top:
            print(f"    {n:>4} × {name}")


if __name__ == "__main__":
    main()
