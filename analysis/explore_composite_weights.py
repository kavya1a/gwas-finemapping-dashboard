"""Explore alternative composite-score weightings using existing modality_scores.

Recomputes per-variant composite scores under several weighting schemes and
reports the canonical variants' ranks under each scheme. No API calls — uses
scored_variants.db only.

Schemes:
  baseline               current weights in scoring/composite.py
  expression_heavy       global re-balance: up-weight expression, down chromatin
  splice_heavy           global re-balance: emphasize splicing
  brain_emphasis         per-disease: AD/SCZ/PD up-weight expression + splicing
  saturation_discount    discount each modality by (1 − its disease-median),
                         so saturated modalities contribute less
  zscore                 per-disease z-score each modality before weighting

Decision rule: a scheme is "promising" if all 4 rank-based canonical tests
(rs429358 AD, rs7903146 T2D, rs1006737 SCZ, rs356219 PD) land in top 20%.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

DIR = Path(__file__).parent.parent
SCORED_DB = DIR / "scored_variants.db"

CANONICAL = [
    ("alzheimers",    "rs429358"),
    ("t2d",           "rs7903146"),
    ("schizophrenia", "rs1006737"),
    ("parkinsons",    "rs356219"),
]

TOP_PCT = 0.20

# ---------------------------------------------------------------------------
# Weighting schemes
# ---------------------------------------------------------------------------

BASELINE = {
    "expression":  0.30,
    "chromatin":   0.25,
    "tf_binding":  0.15,
    "splicing":    0.20,
    "contact":     0.05,
    "polyadenylation": 0.05,
}

EXPRESSION_HEAVY = {
    "expression":  0.45,
    "chromatin":   0.15,
    "tf_binding":  0.15,
    "splicing":    0.20,
    "contact":     0.025,
    "polyadenylation": 0.025,
}

SPLICE_HEAVY = {
    "expression":  0.25,
    "chromatin":   0.20,
    "tf_binding":  0.15,
    "splicing":    0.30,
    "contact":     0.05,
    "polyadenylation": 0.05,
}

BRAIN_EMPHASIS = {
    "alzheimers":    {"expression": 0.40, "chromatin": 0.15, "tf_binding": 0.10,
                      "splicing": 0.25, "contact": 0.05, "polyadenylation": 0.05},
    "schizophrenia": {"expression": 0.40, "chromatin": 0.15, "tf_binding": 0.10,
                      "splicing": 0.25, "contact": 0.05, "polyadenylation": 0.05},
    "parkinsons":    {"expression": 0.40, "chromatin": 0.15, "tf_binding": 0.10,
                      "splicing": 0.25, "contact": 0.05, "polyadenylation": 0.05},
    "t2d":           {"expression": 0.25, "chromatin": 0.35, "tf_binding": 0.15,
                      "splicing": 0.15, "contact": 0.05, "polyadenylation": 0.05},
}

SPLICE_SUBMODALITIES = {"splice_sites", "splice_site_usage", "splice_junctions"}

RARE_VARIANT_DISCOUNT = 0.80  # matches scoring/composite.py:62


# ---------------------------------------------------------------------------
# Composite recomputation
# ---------------------------------------------------------------------------

def _build_wide(conn) -> pd.DataFrame:
    """Wide DF: rows = (disease, rsid); cols = modality scores + rare flag."""
    ms = pd.read_sql_query(
        "SELECT disease, rsid, modality, max_abs_score FROM modality_scores "
        "WHERE max_abs_score IS NOT NULL",
        conn,
    )
    wide = ms.pivot_table(index=["disease", "rsid"], columns="modality",
                          values="max_abs_score", aggfunc="first")
    # Average the splice sub-modalities into a single "splicing" column.
    splice_cols = [c for c in SPLICE_SUBMODALITIES if c in wide.columns]
    if splice_cols:
        wide["splicing"] = wide[splice_cols].mean(axis=1, skipna=True)
        wide = wide.drop(columns=splice_cols)
    wide = wide.reset_index()

    # Rare-variant flag — composite.py multiplies the final score by 0.80 for
    # variants with MAF < 0.01 (RARE_VARIANT_DISCOUNT). Pull the stored flag.
    rare = pd.read_sql_query(
        "SELECT disease, rsid, rare_variant_caution FROM scores "
        "WHERE error IS NULL AND composite_score IS NOT NULL",
        conn,
    )
    wide = wide.merge(rare, on=["disease", "rsid"], how="left")
    wide["rare_variant_caution"] = wide["rare_variant_caution"].fillna(0).astype(int)
    return wide


def _apply_weights(wide: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    """Weighted sum / total configured weight — matches composite.py:107-115.

    Missing modalities contribute 0 to the numerator but the total denominator
    is fixed at sum(weights), so variants with fewer modalities present score
    lower (intentional: it reflects partial coverage)."""
    total_w = sum(weights.values()) or 1.0
    composite = pd.Series(0.0, index=wide.index)
    for m, w in weights.items():
        if m not in wide.columns:
            continue
        s = wide[m].clip(upper=1.0).fillna(0.0)
        composite = composite + w * s
    return composite / total_w


def _apply_rare_discount(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["composite"] = df["composite"] * np.where(
        df["rare_variant_caution"] == 1, RARE_VARIANT_DISCOUNT, 1.0
    )
    return df


def _compute_per_disease(wide: pd.DataFrame,
                         scheme: dict | dict[str, dict]) -> pd.DataFrame:
    """Compute composite per (disease, rsid) for either a single weight dict
    or a per-disease dict of dicts."""
    out = wide.copy()
    if isinstance(next(iter(scheme.values())), dict):
        composites = []
        for disease in out["disease"].unique():
            sub = out[out["disease"] == disease].copy()
            sub["composite"] = _apply_weights(sub, scheme[disease])
            composites.append(sub)
        out = pd.concat(composites, ignore_index=True)
    else:
        out["composite"] = _apply_weights(out, scheme)
    return _apply_rare_discount(out)


def _saturation_discount_weights(wide: pd.DataFrame, disease: str,
                                 base: dict[str, float]) -> dict[str, float]:
    """For one disease, scale each modality weight by (1 − median of that
    modality across the disease set), so saturated modalities contribute less."""
    sub = wide[wide["disease"] == disease]
    out = {}
    for m, w in base.items():
        if m in sub.columns:
            med = sub[m].median(skipna=True)
            scale = max(0.05, 1.0 - (med if pd.notna(med) else 0.0))
            out[m] = w * scale
        else:
            out[m] = w
    return out


def _zscore_composite(wide: pd.DataFrame, weights: dict[str, float]) -> pd.DataFrame:
    """Per-disease z-score each modality, then weighted-sum (fixed denominator)."""
    out = wide.copy()
    total_w = sum(weights.values()) or 1.0
    composites = []
    for disease in out["disease"].unique():
        sub = out[out["disease"] == disease].copy()
        composite = pd.Series(0.0, index=sub.index)
        for m, w in weights.items():
            if m not in sub.columns:
                continue
            s = sub[m]
            mu, sigma = s.mean(skipna=True), s.std(skipna=True)
            if not (sigma > 0):
                continue
            z = ((s - mu) / sigma).fillna(0.0)
            composite = composite + w * z
        sub["composite"] = composite / total_w
        composites.append(sub)
    return pd.concat(composites, ignore_index=True)


# ---------------------------------------------------------------------------
# Rank evaluation
# ---------------------------------------------------------------------------

def _rank_pct(composites: pd.DataFrame, disease: str, rsid: str) -> tuple[int | None, int, float | None]:
    sub = composites[(composites["disease"] == disease) & composites["composite"].notna()]
    sub = sub.sort_values("composite", ascending=False).reset_index(drop=True)
    n = len(sub)
    hit = sub.index[sub["rsid"] == rsid].tolist()
    if not hit:
        return None, n, None
    rank = hit[0] + 1
    return rank, n, rank / n


def _summary(scheme_name: str, composites: pd.DataFrame) -> dict:
    row = {"scheme": scheme_name}
    n_pass = 0
    for disease, rsid in CANONICAL:
        rank, total, pct = _rank_pct(composites, disease, rsid)
        if rank is None:
            row[f"{disease}/{rsid}"] = "missing"
            continue
        flag = "PASS" if pct <= TOP_PCT else "FAIL"
        row[f"{disease}/{rsid}"] = f"{rank}/{total} ({100*pct:.1f}%) {flag}"
        if pct <= TOP_PCT:
            n_pass += 1
    row["n_pass"] = f"{n_pass}/4"
    return row


def main() -> None:
    if not SCORED_DB.exists():
        raise SystemExit(f"{SCORED_DB} not found")

    conn = sqlite3.connect(SCORED_DB)
    try:
        wide = _build_wide(conn)
    finally:
        conn.close()

    print(f"Loaded {len(wide)} (disease, rsid) rows across {wide['disease'].nunique()} diseases")
    print(f"Modality columns: {[c for c in wide.columns if c not in ('disease','rsid')]}")

    rows = []

    # 1. baseline
    rows.append(_summary("baseline", _compute_per_disease(wide, BASELINE)))

    # 2. expression_heavy
    rows.append(_summary("expression_heavy", _compute_per_disease(wide, EXPRESSION_HEAVY)))

    # 3. splice_heavy
    rows.append(_summary("splice_heavy", _compute_per_disease(wide, SPLICE_HEAVY)))

    # 4. brain_emphasis (per-disease)
    rows.append(_summary("brain_emphasis", _compute_per_disease(wide, BRAIN_EMPHASIS)))

    # 5. saturation_discount (per-disease, derived from baseline)
    composites = []
    for disease in wide["disease"].unique():
        w = _saturation_discount_weights(wide, disease, BASELINE)
        sub = wide[wide["disease"] == disease].copy()
        sub["composite"] = _apply_weights(sub, w)
        composites.append(sub)
    rows.append(_summary("saturation_discount",
                         _apply_rare_discount(pd.concat(composites, ignore_index=True))))

    # 6. zscore (per-disease)
    rows.append(_summary("zscore",
                         _apply_rare_discount(_zscore_composite(wide, BASELINE))))

    # 7. t2d_chromatin_heavy: T2D-specific up-weight of chromatin (where rs7903146
    # is strong, 24% beat it), other diseases unchanged.
    scheme_t2d_chrom = {
        "alzheimers":    BASELINE,
        "schizophrenia": BASELINE,
        "parkinsons":    BASELINE,
        "t2d":           {"expression": 0.20, "chromatin": 0.40, "tf_binding": 0.15,
                          "splicing":   0.15, "contact":   0.05, "polyadenylation": 0.05},
    }
    rows.append(_summary("t2d_chromatin_heavy",
                         _compute_per_disease(wide, scheme_t2d_chrom)))

    # 8. bestmod: rank by max(modality_score) — outlier-friendly. Use the
    # rare-variant discount on top for parity with the other schemes.
    out = wide.copy()
    mod_cols = [c for c in ["expression", "chromatin", "tf_binding", "splicing"]
                if c in out.columns]
    out["composite"] = out[mod_cols].max(axis=1, skipna=True)
    rows.append(_summary("bestmod_max", _apply_rare_discount(out)))

    out = pd.DataFrame(rows)
    cols = ["scheme", "n_pass"] + [f"{d}/{r}" for d, r in CANONICAL]
    out = out[cols]
    print("\n" + out.to_string(index=False))


if __name__ == "__main__":
    main()
