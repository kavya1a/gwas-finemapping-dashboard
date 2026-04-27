"""Benchmark orchestrator.

Runs three tasks:
  1. Positive benchmark: Score 100 pathogenic ClinVar variants → recall@K metrics.
  2. Baseline comparisons: CADD PHRED and GWAS p-value rankings on same set.
  3. Negative controls: Score 50 benign variants, report specificity.

All AlphaGenome scores cached in benchmark/cache/ag_scores.json to avoid
re-running the expensive API calls.

Usage:
    /opt/homebrew/bin/python3.11 benchmark/runner.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

import pandas as pd
import numpy as np

from benchmark.clinvar_fetcher import sample_benchmark_set
from benchmark.external_scores import fetch_all_external_scores
from benchmark.metrics import compute_all_metrics, scores_from_df
from scoring.composite import VariantInput, score_variants_composite
from scoring.tissue_config import get_profile

AG_CACHE = Path(__file__).parent / "cache" / "ag_scores.json"


# ---------------------------------------------------------------------------
# AlphaGenome scoring with per-variant caching
# ---------------------------------------------------------------------------

def _load_ag_cache() -> dict:
    if AG_CACHE.exists():
        with open(AG_CACHE) as f:
            return json.load(f)
    return {}


def _save_ag_cache(cache: dict) -> None:
    AG_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(AG_CACHE, "w") as f:
        json.dump(cache, f)


def score_with_ag(
    variants: list[dict],
    batch_size: int = 10,
) -> dict[str, dict]:
    """Score variants with AlphaGenome, using cache to avoid re-running.

    Returns dict: rsid → {"composite_score": float, "disease": str}.
    Variants are batched by disease so tissue profiles are applied correctly.
    """
    cache = _load_ag_cache()

    # Split into disease groups
    by_disease: dict[str, list[dict]] = {}
    for v in variants:
        disease = v.get("disease") or "alzheimers"  # default for benign
        by_disease.setdefault(disease, []).append(v)

    for disease, dvariants in by_disease.items():
        uncached = [v for v in dvariants if v["rsid"] not in cache]
        if not uncached:
            continue

        print(f"\n  Scoring {len(uncached)} uncached {disease} variants...")
        inputs = [
            VariantInput(
                rsid=v["rsid"],
                chrom=v["chrom"],
                pos=v["pos"],
                ref=v["ref"],
                alt=v["alt"],
                maf=v.get("maf"),
            )
            for v in uncached
        ]

        # Process in batches to limit memory
        for start in range(0, len(inputs), batch_size):
            batch = inputs[start : start + batch_size]
            try:
                ranked = score_variants_composite(batch, disease=disease, max_workers=4)
                for _, row in ranked.iterrows():
                    rsid = row["rsid"]
                    cache[rsid] = {
                        "composite_score": (
                            float(row["composite_score"])
                            if pd.notna(row["composite_score"])
                            else None
                        ),
                        "error": row.get("error"),
                        "disease": disease,
                    }
            except Exception as e:
                print(f"    Batch error: {e}")
                for inp in batch:
                    cache[inp.rsid] = {"composite_score": None, "error": str(e), "disease": disease}

        _save_ag_cache(cache)

    return cache


# ---------------------------------------------------------------------------
# Assemble scored DataFrame
# ---------------------------------------------------------------------------

def build_results_df(
    pathogenic: list[dict],
    benign: list[dict],
    ag_cache: dict,
    ext_scores: dict,
) -> pd.DataFrame:
    """Merge all scores into a single analysis DataFrame."""
    rows = []
    for v in pathogenic + benign:
        is_path = v.get("disease") is not None
        rsid = v["rsid"]
        ag = ag_cache.get(rsid, {})
        ext = ext_scores.get(rsid, {})

        cadd = ext.get("cadd_phred")
        gwas_p = ext.get("gwas_pvalue")
        gwas_score = -math.log10(gwas_p) if gwas_p and gwas_p > 0 else 0.0

        rows.append({
            "rsid": rsid,
            "chrom": v["chrom"],
            "pos": v["pos"],
            "gene": v.get("gene", ""),
            "phenotype": v.get("phenotype", "")[:80],
            "clinsig": v.get("clinsig", ""),
            "disease": v.get("disease") or "benign_control",
            "is_pathogenic": 1 if is_path else 0,
            "ag_composite": ag.get("composite_score"),
            "cadd_phred": cadd if cadd is not None else 0.0,
            "gwas_neg_log10p": gwas_score,
            "ag_error": ag.get("error"),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_benchmark(
    n_pathogenic: int = 100,
    n_benign: int = 50,
    seed: int = 42,
) -> dict:
    print("=" * 70)
    print("Benchmark: AlphaGenome vs CADD vs GWAS p-value")
    print(f"Target: {n_pathogenic} pathogenic, {n_benign} benign ClinVar variants")
    print("=" * 70)

    # 1. Load/sample benchmark set
    print("\n[1/4] Loading ClinVar benchmark variants...")
    pathogenic, benign = sample_benchmark_set(n_pathogenic, n_benign, seed=seed)
    all_variants = pathogenic + benign

    print(f"  Pathogenic by disease:")
    from collections import Counter
    for d, cnt in Counter(v["disease"] for v in pathogenic).items():
        print(f"    {d}: {cnt}")

    # 2. External scores (CADD + GWAS p-value)
    print(f"\n[2/4] Fetching CADD + GWAS p-values for {len(all_variants)} variants...")
    ext_scores = fetch_all_external_scores(all_variants)
    cadd_available = sum(1 for e in ext_scores.values() if e.get("cadd_phred") is not None)
    gwas_available = sum(1 for e in ext_scores.values() if e.get("gwas_pvalue") is not None)
    print(f"  CADD scores available: {cadd_available}/{len(all_variants)}")
    print(f"  GWAS p-values available: {gwas_available}/{len(all_variants)}")

    # 3. AlphaGenome scoring
    print(f"\n[3/4] AlphaGenome scoring ({len(all_variants)} variants)...")
    ag_cache = score_with_ag(all_variants)
    ag_scored = sum(1 for r in ag_cache.values() if r.get("composite_score") is not None)
    print(f"  AlphaGenome scores available: {ag_scored}/{len(all_variants)}")

    # 4. Build combined DataFrame
    print("\n[4/4] Computing metrics...")
    df = build_results_df(pathogenic, benign, ag_cache, ext_scores)

    # Filter to variants with AlphaGenome scores (drop errors)
    df_ag = df[df["ag_composite"].notna()].copy()
    path_ag, benign_ag = scores_from_df(df_ag, "ag_composite")

    # CADD metrics (all variants; use 0.0 for missing = conservative)
    path_cadd, benign_cadd = scores_from_df(df, "cadd_phred")

    # GWAS metrics (all variants; use 0.0 for no GWAS signal = not associated)
    path_gwas, benign_gwas = scores_from_df(df, "gwas_neg_log10p")

    metrics = {
        "alphagenome_composite": compute_all_metrics(
            path_ag, benign_ag, label="AlphaGenome Composite"
        ),
        "cadd_phred": compute_all_metrics(
            path_cadd, benign_cadd, label="CADD PHRED v1.7"
        ),
        "gwas_neg_log10p": compute_all_metrics(
            path_gwas, benign_gwas, label="GWAS -log10(p-value)"
        ),
    }

    # Negative control specificity breakdown per disease
    nc_by_disease = {}
    for disease in ["alzheimers", "parkinsons", "t2d", "schizophrenia"]:
        sub_path = df_ag[
            (df_ag["is_pathogenic"] == 1) & (df_ag["disease"] == disease)
        ]["ag_composite"].dropna().tolist()
        sub_benign = df_ag[df_ag["is_pathogenic"] == 0]["ag_composite"].dropna().tolist()
        if sub_path and sub_benign:
            from benchmark.metrics import negative_control_specificity, recall_at_k, auroc
            nc_by_disease[disease] = {
                "n_path": len(sub_path),
                "median_path_score": round(float(np.median(sub_path)), 4),
                "auroc": round(auroc(sub_path, sub_benign), 4),
                **recall_at_k(sub_path, sub_benign),
            }

    results = {
        "summary": {
            "n_pathogenic_scored": len(path_ag),
            "n_benign_scored": len(benign_ag),
            "n_pathogenic_total": n_pathogenic,
            "n_benign_total": n_benign,
            "seed": seed,
        },
        "metrics": metrics,
        "per_disease": nc_by_disease,
        "variant_scores": df.to_dict(orient="records"),
    }

    # Save results
    out_path = Path(__file__).parent.parent / "benchmark_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    results = run_benchmark()

    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    for method, m in results["metrics"].items():
        print(f"\n{m['label']}:")
        print(f"  Recall@1 / @5 / @10 = "
              f"{m.get('recall_at_1', 'N/A')} / "
              f"{m.get('recall_at_5', 'N/A')} / "
              f"{m.get('recall_at_10', 'N/A')}")
        print(f"  auROC = {m['auROC']}  AUPRC = {m['AUPRC']}")
        print(f"  Neg ctrl frac > P75 = {m['neg_ctrl_frac_above_p75']}")
        print(f"  Median pathogenic={m['median_pathogenic_score']}  "
              f"benign={m['median_benign_score']}")

    print("\nPer-disease AlphaGenome auROC:")
    for d, dm in results["per_disease"].items():
        print(f"  {d:<14} n={dm['n_path']}  "
              f"auROC={dm['auroc']}  "
              f"R@1={dm.get('recall_at_1','?')}  "
              f"R@5={dm.get('recall_at_5','?')}")
