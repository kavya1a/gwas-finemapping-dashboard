"""End-to-end test: score 10 Alzheimer's disease GWAS variants.

Variants are drawn from published GWAS meta-analyses (Lambert et al. 2013,
Kunkle et al. 2019) with hg38 positions from Ensembl/dbSNP.
APOE rs429358 is variant #1 and must appear in the top 5 for the ranking
to be biologically plausible.

Run with:
    /opt/homebrew/bin/python3.11 run_ad_test.py
"""

import sys
import os
from pathlib import Path

# Ensure project root is on path
sys.path.insert(0, str(Path(__file__).parent))

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent / ".env")

import pandas as pd

from scoring.composite import VariantInput, score_variants_composite

# ---------------------------------------------------------------------------
# Curated Alzheimer's GWAS variants (hg38 / GRCh38)
# Positions and alleles from dbSNP / Ensembl / NHGRI-EBI GWAS Catalog.
# ---------------------------------------------------------------------------
AD_VARIANTS: list[VariantInput] = [
    # APOE locus — by far the largest common-variant AD signal
    VariantInput("rs429358", "chr19", 44908684, "T", "C",   # APOE ε4, OR ~3.5
                 p_value=1e-700),
    VariantInput("rs7412",   "chr19", 44908822, "C", "T",   # APOE ε2, protective
                 p_value=1e-50),

    # TREM2 — rare coding variant with large effect (R47H)
    VariantInput("rs75932628", "chr6", 41163199, "G", "A",  # TREM2 R47H
                 p_value=3.7e-25),

    # Genome-wide significant non-APOE loci
    VariantInput("rs6656401",  "chr1",  207692049, "G", "A",  # CR1
                 p_value=3.7e-21),
    VariantInput("rs744373",   "chr2",  127892810, "C", "T",  # BIN1
                 p_value=1.6e-11),
    VariantInput("rs9331896",  "chr8",   27467686, "C", "T",  # CLU
                 p_value=8.5e-10),
    VariantInput("rs3865444",  "chr19",  51728299, "C", "A",  # CD33
                 p_value=1.6e-9),
    VariantInput("rs10948363", "chr6",   47571792, "G", "A",  # CD2AP
                 p_value=8.6e-9),
    VariantInput("rs3764650",  "chr19",   1050228, "G", "T",  # ABCA7
                 p_value=4.5e-17),
    VariantInput("rs10993776", "chr11",  85868640, "A", "G",  # PICALM
                 p_value=3.2e-9),
]


def print_full_breakdown(rsid: str, breakdown_df: pd.DataFrame) -> None:
    if breakdown_df.empty:
        print("    (no modality breakdown available)")
        return
    print(f"\n  Modality breakdown for {rsid}:")
    cols = ["modality", "max_abs_score", "mean_abs_score",
            "frac_above_threshold", "signed_max_score", "top_track_biosample"]
    for _, row in breakdown_df.sort_values("max_abs_score", ascending=False).iterrows():
        sign = "▲" if row["signed_max_score"] >= 0 else "▼"
        bio = row["top_track_biosample"] or "—"
        print(
            f"    {row['modality']:<18}  max={row['max_abs_score']:.3f}  "
            f"mean={row['mean_abs_score']:.3f}  "
            f"frac>{0.90:.0%}={row['frac_above_threshold']:.2f}  "
            f"dir={sign}  top_tissue={bio}"
        )


def main():
    print("=" * 70)
    print("AlphaGenome end-to-end test: 10 Alzheimer's Disease variants")
    print("Tissue profile: brain/neuronal (Alzheimer's)")
    print("=" * 70)

    print(f"\nScoring {len(AD_VARIANTS)} variants...\n")
    ranked = score_variants_composite(AD_VARIANTS, disease="alzheimers")

    print("\n" + "=" * 70)
    print("RANKED RESULTS")
    print("=" * 70)
    display_cols = ["rank", "rsid", "chrom", "pos", "ref", "alt",
                    "composite_score", "error"]
    print(ranked[display_cols].to_string(index=False))

    print("\n" + "=" * 70)
    print("TOP 3 — FULL MODALITY/TISSUE BREAKDOWN")
    print("=" * 70)
    top3 = ranked.head(3)
    for _, row in top3.iterrows():
        print(f"\n#{int(row['rank'])}  {row['rsid']}  "
              f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}  "
              f"composite={row['composite_score']:.4f}")
        bd = ranked.modality_breakdowns.get(row["rsid"], pd.DataFrame())
        print_full_breakdown(row["rsid"], bd)

    print("\n" + "=" * 70)
    print("APOE rs429358 DIAGNOSIS")
    print("=" * 70)
    apoe_rows = ranked[ranked["rsid"] == "rs429358"]
    if apoe_rows.empty:
        print("rs429358 was NOT scored (likely missing ref/alt).")
    else:
        apoe = apoe_rows.iloc[0]
        rank = int(apoe["rank"])
        score = apoe["composite_score"]
        print(f"rs429358 rank: #{rank}  composite_score: {score:.4f}")
        if rank <= 5:
            print("PASS — APOE ε4 is in the top 5. Ranking is biologically plausible.")
        else:
            print(f"INVESTIGATE — APOE ε4 ranked #{rank}. Modality breakdown:")
            bd = ranked.modality_breakdowns.get("rs429358", pd.DataFrame())
            print_full_breakdown("rs429358", bd)
            print("\nPossible causes:")
            print("  - rs429358 is an APOE coding variant (chr19 APOE exon 4);")
            print("    its primary mechanism is protein-level (ApoE isoforms),")
            print("    which AlphaGenome's regulatory predictions won't capture fully.")
            print("  - Expected top signals: RNA_SEQ (expression change), CAGE")
            print("  - If expression/chromatin scores are near zero, it confirms")
            print("    the variant acts post-translationally and is expected to rank low")
            print("    on purely regulatory metrics.")


if __name__ == "__main__":
    main()
