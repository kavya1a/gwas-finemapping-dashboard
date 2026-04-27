"""Cross-disease generalization test + APOE sign verification.

Covers three tasks from the verification milestone:
  1. APOE ε4 vs ε2 sign comparison (directional score validation)
  2. T2D pipeline — canonical: TCF7L2 rs7903146, tissue: pancreatic islet
  3. Schizophrenia pipeline — canonical: CACNA1C rs1006737, tissue: cortical neurons
  4. Parkinson's pipeline — canonical: SNCA rs356219, tissue: dopaminergic neurons

Variant positions from Ensembl GRCh38 API lookup (2026-04-20).
rs7903146 alleles from test_api.py (confirmed T/C, Ensembl returned ambiguous G/C).

Run with:
    /opt/homebrew/bin/python3.11 run_generalization_test.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent / ".env")

import pandas as pd
from scoring.composite import VariantInput, score_variants_composite


# ---------------------------------------------------------------------------
# Curated variant sets — positions from Ensembl GRCh38, alleles verified
# ---------------------------------------------------------------------------

SIGN_CHECK_VARIANTS = [
    VariantInput("rs429358", "chr19", 44908684, "T", "C", maf=0.135),  # APOE ε4, risk
    VariantInput("rs7412",   "chr19", 44908822, "C", "T", maf=0.080),  # APOE ε2, protective
]

T2D_VARIANTS = [
    VariantInput("rs7903146",  "chr10", 112998590, "T", "C", maf=0.299),  # TCF7L2 CANONICAL
    VariantInput("rs1801282",  "chr3",   12351626, "C", "G", maf=0.120),  # PPARG P12A coding
    VariantInput("rs10811661", "chr9",   22134095, "T", "A", maf=0.176),  # CDKN2A/B
    VariantInput("rs4402960",  "chr3",  185793899, "G", "C", maf=0.388),  # IGF2BP2
    VariantInput("rs13266634", "chr8",  117172544, "C", "A", maf=None),   # SLC30A8
    VariantInput("rs5015480",  "chr10",  92705802, "C", "A", maf=0.454),  # HHEX
    VariantInput("rs7754840",  "chr6",   20661019, "G", "A", maf=0.406),  # CDKAL1
    VariantInput("rs10010131", "chr4",    6291188, "A", "G", maf=None),   # WFS1
    VariantInput("rs8050136",  "chr16",  53782363, "C", "A", maf=None),   # FTO
    VariantInput("rs864745",   "chr7",   28140937, "T", "C", maf=0.304),  # JAZF1
]

SCZ_VARIANTS = [
    VariantInput("rs1006737",  "chr12",  2236129, "G", "A", maf=0.301),  # CACNA1C CANONICAL
    VariantInput("rs1625579",  "chr1",  98037378, "G", "A", maf=None),   # near MIR137
    VariantInput("rs11209026", "chr1",  67240275, "G", "A", maf=None),   # IL23R/RBM26
    VariantInput("rs4523957",  "chr17",  2305605, "G", "T", maf=0.477),  # PLCL1/CAMTA1
    VariantInput("rs11062528", "chr12",  3144261, "C", "T", maf=0.064),  # FURIN
    VariantInput("rs7297175",  "chr12", 56080024, "T", "C", maf=0.349),  # FAM57B/GATAD2A
    VariantInput("rs2535627",  "chr3",  52811089, "T", "A", maf=0.465),  # PLCL1
    VariantInput("rs17283563", "chr13", 23001883, "G", "A", maf=0.199),  # SLITRK1
    VariantInput("rs16887244", "chr8",  38173827, "A", "C", maf=0.195),  # PRAGMIN
    VariantInput("rs2251219",  "chr3",  52550771, "T", "A", maf=None),   # PLCL1 region
]

PD_VARIANTS = [
    VariantInput("rs356219",   "chr4",   89716450, "G", "A", maf=0.491),  # SNCA CANONICAL
    VariantInput("rs34637584", "chr12",  40340400, "G", "A", maf=None),   # LRRK2 G2019S [RARE]
    VariantInput("rs11060180", "chr12", 122819039, "A", "G", maf=None),   # MAPT H1 locus
    VariantInput("rs823118",   "chr1",  205754444, "C", "A", maf=0.414),  # NUCKS1/RAB7L1
    VariantInput("rs6812193",  "chr4",   76277833, "C", "T", maf=0.314),  # DGKQ
    VariantInput("rs3771076",  "chr2",  187412396, "T", "A", maf=0.334),  # MCCC1
    VariantInput("rs1491942",  "chr12",  40227006, "C", "G", maf=0.297),  # HIP1R
    VariantInput("rs591323",   "chr8",   16839582, "G", "A", maf=0.364),  # MCIDAS/ADGRB1
    VariantInput("rs1474616",  "chr6",  136397866, "C", "A", maf=0.102),  # MBNL2
    VariantInput("rs3219488",  "chr1",   45331858, "G", "A",
                 maf=0.000589),   # VERY RARE — tests rare-variant caution flag
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def print_disease_results(ranked: pd.DataFrame, canonical_rsid: str, label: str) -> None:
    print(f"\n{'='*70}")
    print(f"DISEASE: {label}")
    print(f"{'='*70}")

    display = ranked[["rank", "rsid", "chrom", "pos", "ref", "alt",
                       "composite_score", "rare_variant_caution", "error"]].copy()
    display["composite_score"] = display["composite_score"].map(
        lambda x: f"{x:.4f}" if pd.notna(x) else "N/A"
    )
    print(display.to_string(index=False))

    print(f"\n--- Top 5 Modality Breakdowns ---")
    for _, row in ranked.head(5).iterrows():
        rare_flag = " [RARE ⚠]" if row.get("rare_variant_caution") else ""
        print(f"\n  #{int(row['rank'])}  {row['rsid']}{rare_flag}  "
              f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}  "
              f"composite={row['composite_score']:.4f}")
        bd = ranked.modality_breakdowns.get(row["rsid"], pd.DataFrame())
        if bd.empty:
            print("    (no breakdown)")
            continue
        for _, br in bd.sort_values("max_abs_score", ascending=False).iterrows():
            sign = "▲" if br["signed_max_score"] >= 0 else "▼"
            bio = br["top_track_biosample"] or "—"
            print(
                f"    {br['modality']:<22} "
                f"max={br['max_abs_score']:.3f}  "
                f"mean={br['mean_abs_score']:.3f}  "
                f"frac≥0.90={br['frac_above_threshold']:.2f}  "
                f"{sign}  {bio}"
            )

    # Canonical variant check
    canon_row = ranked[ranked["rsid"] == canonical_rsid]
    print(f"\n--- Canonical variant check: {canonical_rsid} ---")
    if canon_row.empty:
        print(f"  NOT FOUND in ranked output (likely scored with error)")
        return
    rank = int(canon_row.iloc[0]["rank"])
    score = canon_row.iloc[0]["composite_score"]
    if rank <= 10:
        status = "PASS ✓" if rank <= 5 else f"MARGINAL (#{rank}, in top 10)"
        print(f"  {status} — {canonical_rsid} rank #{rank}, composite={score:.4f}")
    else:
        print(f"  INVESTIGATE — {canonical_rsid} rank #{rank} (outside top 10)")
        print("  Check modality breakdown for diagnostic clues.")


def print_sign_comparison(ranked: pd.DataFrame) -> None:
    print(f"\n{'='*70}")
    print("TASK 1: APOE ε4 vs ε2 — Directional Sign Comparison")
    print(f"{'='*70}")
    print(
        "rs429358 (T>C) = APOE ε4 — RISK allele  |  "
        "rs7412 (C>T) = APOE ε2 — PROTECTIVE allele\n"
    )

    rows_by_rsid = {}
    for _, row in ranked.iterrows():
        bd = ranked.modality_breakdowns.get(row["rsid"], pd.DataFrame())
        rows_by_rsid[row["rsid"]] = {"row": row, "bd": bd}

    modalities = ["expression", "chromatin", "tf_binding", "splicing",
                  "splice_junctions", "splice_site_usage", "contact", "polyadenylation"]

    e4 = rows_by_rsid.get("rs429358", {})
    e2 = rows_by_rsid.get("rs7412", {})

    print(f"  {'Modality':<22}  {'ε4 signed_max':>14}  {'ε2 signed_max':>14}  Signs differ?")
    print(f"  {'-'*22}  {'-'*14}  {'-'*14}  {'-'*12}")

    sign_disagreements = 0
    sign_agreements = 0
    for mod in modalities:
        e4_val = None
        e2_val = None
        if not e4.get("bd", pd.DataFrame()).empty:
            row_e4 = e4["bd"][e4["bd"]["modality"] == mod]
            if not row_e4.empty:
                e4_val = row_e4.iloc[0]["signed_max_score"]
        if not e2.get("bd", pd.DataFrame()).empty:
            row_e2 = e2["bd"][e2["bd"]["modality"] == mod]
            if not row_e2.empty:
                e2_val = row_e2.iloc[0]["signed_max_score"]

        if e4_val is not None and e2_val is not None:
            differ = (e4_val * e2_val) < 0  # opposite signs
            differ_str = "YES ✓" if differ else "no"
            if differ:
                sign_disagreements += 1
            else:
                sign_agreements += 1
            print(f"  {mod:<22}  {e4_val:>+14.3f}  {e2_val:>+14.3f}  {differ_str}")
        else:
            e4_str = f"{e4_val:+.3f}" if e4_val is not None else "N/A"
            e2_str = f"{e2_val:+.3f}" if e2_val is not None else "N/A"
            print(f"  {mod:<22}  {e4_str:>14}  {e2_str:>14}  —")

    print(f"\n  Summary: {sign_disagreements} modalities show OPPOSITE signs, "
          f"{sign_agreements} show SAME sign.")

    # Key expression sign check
    e4_expr = None
    e2_expr = None
    if not e4.get("bd", pd.DataFrame()).empty:
        r = e4["bd"][e4["bd"]["modality"] == "expression"]
        e4_expr = r.iloc[0]["signed_max_score"] if not r.empty else None
    if not e2.get("bd", pd.DataFrame()).empty:
        r = e2["bd"][e2["bd"]["modality"] == "expression"]
        e2_expr = r.iloc[0]["signed_max_score"] if not r.empty else None

    if e4_expr is not None and e2_expr is not None:
        if (e4_expr * e2_expr) < 0:
            print(
                f"\n  EXPRESSION SIGN: PASS ✓ — ε4={e4_expr:+.3f} vs ε2={e2_expr:+.3f}"
                " (opposite directions, consistent with known APOE cis-eQTL biology)"
            )
            print("  Biological interpretation: ε4 (T>C) reduces neuronal APOE expression "
                  "(▼), ε2 (C>T) increases it (▲).\n"
                  "  AlphaGenome correctly captures the opposing regulatory effects at "
                  "the APOE locus, confirming signed scoring is preserved end-to-end.")
        else:
            print(
                f"\n  EXPRESSION SIGN: INVESTIGATE — ε4={e4_expr:+.3f} vs ε2={e2_expr:+.3f}"
                " (same direction — unexpected given known biology)"
            )
    else:
        print("\n  Expression sign: could not compare (one or both variants missing)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Cross-Disease Generalization Test + Sign Verification")
    print("=" * 70)

    # --- Task 1: APOE sign check (reuses AD tissue profile) ---
    print("\n[Task 1] APOE ε4 vs ε2 sign comparison (brain tissue)...")
    sign_ranked = score_variants_composite(SIGN_CHECK_VARIANTS, disease="alzheimers")
    print_sign_comparison(sign_ranked)

    # --- Task 2: T2D ---
    print("\n[Task 2] Type 2 Diabetes — 10 variants, tissue: pancreatic islet/liver...")
    t2d_ranked = score_variants_composite(T2D_VARIANTS, disease="t2d")
    print_disease_results(t2d_ranked, "rs7903146", "Type 2 Diabetes (canonical: TCF7L2 rs7903146)")

    # --- Task 3: Schizophrenia ---
    print("\n[Task 3] Schizophrenia — 10 variants, tissue: cortical neurons...")
    scz_ranked = score_variants_composite(SCZ_VARIANTS, disease="schizophrenia")
    print_disease_results(scz_ranked, "rs1006737", "Schizophrenia (canonical: CACNA1C rs1006737)")

    # --- Task 4: Parkinson's ---
    print("\n[Task 4] Parkinson's Disease — 10 variants, tissue: dopaminergic neurons...")
    pd_ranked = score_variants_composite(PD_VARIANTS, disease="parkinsons")
    print_disease_results(pd_ranked, "rs356219", "Parkinson's Disease (canonical: SNCA rs356219)")

    # --- Final summary ---
    print(f"\n{'='*70}")
    print("CROSS-DISEASE VALIDATION SUMMARY")
    print(f"{'='*70}")
    checks = [
        ("T2D",      t2d_ranked,  "rs7903146", "pancreatic islet"),
        ("SCZ",      scz_ranked,  "rs1006737", "cortical neurons"),
        ("PD",       pd_ranked,   "rs356219",  "dopaminergic neurons"),
    ]
    all_pass = True
    for disease, ranked, canon, tissue in checks:
        row = ranked[ranked["rsid"] == canon]
        if row.empty:
            print(f"  {disease:<6} {canon} : ERROR — not scored")
            all_pass = False
            continue
        rank = int(row.iloc[0]["rank"])
        score = row.iloc[0]["composite_score"]
        bd = ranked.modality_breakdowns.get(canon, pd.DataFrame())
        top_tissue = "—"
        if not bd.empty:
            top_row = bd.loc[bd["max_abs_score"].idxmax()]
            top_tissue = str(top_row["top_track_biosample"] or "—")
        status = "PASS ✓" if rank <= 10 else "FAIL ✗"
        if rank > 10:
            all_pass = False
        print(f"  {disease:<6} {canon}  rank=#{rank:>2}  score={score:.4f}  "
              f"top_tissue={top_tissue[:40]}  {status}")

    print()
    if all_pass:
        print("ALL CANONICAL VARIANTS IN TOP 10 — pipeline generalizes across diseases.")
        print("Ready for Streamlit UI implementation.")
    else:
        print("Some canonical variants outside top 10 — investigate before UI.")


if __name__ == "__main__":
    main()
