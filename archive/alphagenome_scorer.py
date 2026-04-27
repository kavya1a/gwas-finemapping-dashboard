"""Public interface for AlphaGenome variant scoring.

Thin adapter between app.py and scoring/composite.py.
Checks scored_variants.db cache first; only calls the API for uncached variants.
"""

import os
import sqlite3
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from scoring.composite import VariantInput, score_variants_composite

load_dotenv(Path(__file__).parent / ".env")

SCORED_DB = Path(__file__).parent / "scored_variants.db"


def _load_cached(rsids: list[str], disease: str) -> dict[str, dict]:
    """Return cached score rows keyed by rsid."""
    if not SCORED_DB.exists():
        return {}
    conn = sqlite3.connect(SCORED_DB)
    try:
        placeholders = ",".join("?" * len(rsids))
        rows = conn.execute(
            f"SELECT rsid,chrom,pos,ref,alt,maf,composite_score,pip_weighted_score,"
            f"pip,rare_variant_caution,error,rank "
            f"FROM scores WHERE disease=? AND rsid IN ({placeholders})",
            [disease, *rsids],
        ).fetchall()
        cols = [
            "rsid", "chrom", "pos", "ref", "alt", "maf",
            "composite_score", "pip_weighted_score", "pip",
            "rare_variant_caution", "error", "rank",
        ]
        cached = {}
        for row in rows:
            d = dict(zip(cols, row))
            d["rare_variant_caution"] = bool(d["rare_variant_caution"])
            cached[d["rsid"]] = d

        # Load modality breakdowns for cached variants
        if cached:
            rs_list = list(cached.keys())
            ph2 = ",".join("?" * len(rs_list))
            mod_rows = conn.execute(
                f"SELECT rsid,modality,max_abs_score,frac_above_threshold "
                f"FROM modality_scores WHERE disease=? AND rsid IN ({ph2})",
                [disease, *rs_list],
            ).fetchall()
            breakdowns: dict[str, list] = {r: [] for r in rs_list}
            for rsid, modality, max_score, frac in mod_rows:
                breakdowns[rsid].append({
                    "modality": modality,
                    "max_abs_score": max_score,
                    "frac_above_threshold": frac,
                })
            for rsid, rows_list in breakdowns.items():
                cached[rsid]["_modality_breakdown"] = pd.DataFrame(rows_list)

        return cached
    finally:
        conn.close()


def score_variants(
    variants: list[dict],
    disease: str = "",
    tissue: str = "",
) -> pd.DataFrame:
    """Score variants, using scored_variants.db cache where available.

    Args:
        variants: List of dicts with keys rsid, chrom, pos, ref, alt,
            p_value. ref/alt may be None.
        disease: Disease slug for tissue selection and cache lookup.
        tissue: Ignored legacy parameter.

    Returns:
        DataFrame sorted by composite_score.
    """
    scoreable = [v for v in variants if v.get("ref") and v.get("alt")]
    if not scoreable:
        return pd.DataFrame()

    rsids = [v["rsid"] for v in scoreable]
    cached = _load_cached(rsids, disease)

    need_scoring = [v for v in scoreable if v["rsid"] not in cached]

    live_df = pd.DataFrame()
    if need_scoring:
        inputs = []
        for v in need_scoring:
            inputs.append(VariantInput(
                rsid=v.get("rsid", ""),
                chrom=str(v["chrom"]) if not str(v["chrom"]).startswith("chr") else v["chrom"],
                pos=int(v["pos"]),
                ref=v["ref"],
                alt=v["alt"],
                pip=v.get("pip"),
                p_value=v.get("p_value"),
                maf=v.get("maf"),
            ))
        live_df = score_variants_composite(inputs, disease=disease)

    # Build combined rows from cache + live
    all_rows = []
    mod_breakdowns = {}

    for rsid in rsids:
        if rsid in cached:
            row = {k: v for k, v in cached[rsid].items() if k != "_modality_breakdown"}
            all_rows.append(row)
            mod_breakdowns[rsid] = cached[rsid].get("_modality_breakdown", pd.DataFrame())
        elif not live_df.empty:
            match = live_df[live_df["rsid"] == rsid]
            if not match.empty:
                r = match.iloc[0].to_dict()
                all_rows.append({k: v for k, v in r.items()})
                if hasattr(live_df, "modality_breakdowns"):
                    mod_breakdowns[rsid] = live_df.modality_breakdowns.get(rsid, pd.DataFrame())

    if not all_rows:
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)
    sort_col = "pip_weighted_score" if df.get("pip") is not None and df["pip"].notna().any() else "composite_score"
    if sort_col in df.columns:
        df = df.sort_values(sort_col, ascending=False, na_position="last").reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    df._metadata = ["modality_breakdowns"]
    df.modality_breakdowns = mod_breakdowns

    return df
