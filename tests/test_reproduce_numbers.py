"""Smoke suite: every headline number in the README, recomputed from committed cache.

Self-contained — reads only the committed parquet/SQLite artifacts, needs no
AlphaGenome API key and does not import the SDK. If any of these assertions fail,
a number in the README no longer matches the data that is supposed to produce it.

Run:  python -m pytest tests/test_reproduce_numbers.py -q
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent

# Tolerances
RHO_TOL = 2e-3      # Spearman ρ absolute tolerance
PCT_TOL = 0.6       # saturation percentage-point tolerance


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def null_df() -> pd.DataFrame:
    return pd.read_parquet(ROOT / "matched_calibration_null.parquet")


@pytest.fixture(scope="module")
def tewhey() -> pd.DataFrame:
    parq = pd.read_parquet(ROOT / "tewhey_mpra.parquet")[
        ["rsid", "mpra_lfc", "expression_subscore"]
    ]
    con = sqlite3.connect(ROOT / "tewhey_raw_delta_cache.db")
    raw = pd.read_sql("SELECT rsid, max_signed_raw FROM raw_deltas", con)
    con.close()
    return parq.merge(raw, on="rsid", how="inner")


def _matched_quantile(x: np.ndarray, null_sorted: np.ndarray) -> np.ndarray:
    n_lt = np.searchsorted(null_sorted, x, side="left")
    n_eq = np.searchsorted(null_sorted, x, side="right") - n_lt
    return 2.0 * ((n_lt + 0.5 * n_eq) / len(null_sorted)) - 1.0


# ── shapes ────────────────────────────────────────────────────────────────────

def test_dataframe_shapes(null_df, tewhey):
    assert len(null_df) == 5933                       # clean matched null
    assert pd.read_parquet(ROOT / "tewhey_mpra.parquet").shape[0] == 3301
    con = sqlite3.connect(ROOT / "tewhey_raw_delta_cache.db")
    n_raw = con.execute("SELECT COUNT(*) FROM raw_deltas").fetchone()[0]
    con.close()
    assert n_raw == 3301


# ── matched null saturation (README §Methods smoking-gun) ─────────────────────

def test_matched_null_saturation(null_df):
    x = null_df["raw_max_signed_delta"].dropna().to_numpy(float)
    assert len(x) == 5933
    assert np.mean(np.abs(x) > 0.9) * 100 == pytest.approx(0.42, abs=PCT_TOL)
    assert np.mean(np.abs(x) > 0.5) * 100 == pytest.approx(1.28, abs=PCT_TOL)

    for col in ("raw_mean_signed_delta", "raw_median_signed_delta"):
        m = null_df[col].dropna().to_numpy(float)
        assert np.mean(np.abs(m) > 0.9) == 0.0
        assert np.mean(np.abs(m) > 0.5) == 0.0

    q = null_df["single_track_quantile"].dropna().to_numpy(float)
    assert np.mean(np.abs(q) > 0.9) * 100 == pytest.approx(41.4, abs=PCT_TOL)


# ── Tewhey headline: ρ 0.037 → 0.123, n = 3,246 ───────────────────────────────

def test_tewhey_spearman_headline(null_df, tewhey):
    null_sorted = np.sort(null_df["raw_max_signed_delta"].dropna().to_numpy(float))
    d = tewhey.dropna(subset=["mpra_lfc", "expression_subscore", "max_signed_raw"]).copy()
    assert len(d) == 3246

    d["matched_q"] = _matched_quantile(d["max_signed_raw"].to_numpy(float), null_sorted)

    rho_orig, _ = spearmanr(d["expression_subscore"], d["mpra_lfc"])
    rho_matched, _ = spearmanr(d["matched_q"], d["mpra_lfc"])
    rho_raw, _ = spearmanr(d["max_signed_raw"], d["mpra_lfc"])

    assert rho_orig == pytest.approx(0.0367, abs=RHO_TOL)
    assert rho_matched == pytest.approx(0.1225, abs=RHO_TOL)
    assert rho_raw == pytest.approx(0.1225, abs=RHO_TOL)
    # raw and matched are monotone transforms -> identical rank correlation
    assert rho_matched == pytest.approx(rho_raw, abs=1e-3)


def test_tewhey_saturation(null_df, tewhey):
    null_sorted = np.sort(null_df["raw_max_signed_delta"].dropna().to_numpy(float))
    d = tewhey.dropna(subset=["mpra_lfc", "expression_subscore", "max_signed_raw"]).copy()
    d["matched_q"] = _matched_quantile(d["max_signed_raw"].to_numpy(float), null_sorted)
    assert np.mean(np.abs(d["expression_subscore"]) > 0.9) * 100 == pytest.approx(94.9, abs=PCT_TOL)
    assert np.mean(np.abs(d["matched_q"]) > 0.9) * 100 == pytest.approx(12.9, abs=PCT_TOL)
    assert np.mean(np.abs(d["max_signed_raw"]) > 0.9) * 100 == pytest.approx(1.1, abs=PCT_TOL)


def test_full_panel_subscore_saturation():
    parq = pd.read_parquet(ROOT / "tewhey_mpra.parquet")
    q = parq["expression_subscore"].dropna().to_numpy(float)
    assert len(q) == 3259
    assert np.mean(np.abs(q) > 0.9) * 100 == pytest.approx(94.9, abs=PCT_TOL)


# ── disease-modality saturation gradient (README Findings) ────────────────────

def test_disease_modality_saturation():
    con = sqlite3.connect(ROOT / "scored_variants.db")
    df = pd.read_sql(
        "SELECT modality, signed_max_score FROM modality_scores "
        "WHERE signed_max_score IS NOT NULL",
        con,
    )
    con.close()
    expected = {"expression": (767, 99.6), "chromatin": (767, 69.4), "tf_binding": (767, 31.9)}
    for modality, (n_exp, sat_exp) in expected.items():
        s = df.loc[df["modality"] == modality, "signed_max_score"].to_numpy(float)
        assert len(s) == n_exp
        assert np.mean(np.abs(s) > 0.9) * 100 == pytest.approx(sat_exp, abs=PCT_TOL)


# ── blood-trait replication (README Findings) ─────────────────────────────────

def test_blood_trait_saturation():
    con = sqlite3.connect(ROOT / "phase3_blood_cache.db")
    df = pd.read_sql(
        "SELECT trait, expression_subscore FROM scores WHERE expression_subscore IS NOT NULL",
        con,
    )
    con.close()
    expected = {"hemoglobin": (195, 100.0), "platelet_count": (198, 99.0)}
    for trait, (n_exp, sat_exp) in expected.items():
        s = df.loc[df["trait"] == trait, "expression_subscore"].to_numpy(float)
        assert len(s) == n_exp
        assert np.mean(np.abs(s) > 0.9) * 100 == pytest.approx(sat_exp, abs=PCT_TOL)


# ── committed CSV tables match the data / README ──────────────────────────────

def test_matched_calibration_comparison_csv():
    csv = pd.read_csv(ROOT / "matched_calibration_comparison.csv")
    assert (csv["n"] == 3246).all()
    by = {r.predictor: r for r in csv.itertuples()}
    orig = next(v for k, v in by.items() if k.startswith("Original"))
    matched = next(v for k, v in by.items() if "Matched" in k)
    assert orig.spearman_r == pytest.approx(0.0367, abs=RHO_TOL)
    assert matched.spearman_r == pytest.approx(0.1225, abs=RHO_TOL)


def test_three_recipe_csv():
    csv = pd.read_csv(ROOT / "matched_recipes_comparison.csv").set_index("recipe")
    assert csv.loc["max", "null_n"] == 5933
    assert csv.loc["max", "null_abs_gt_0_9"] * 100 == pytest.approx(0.42, abs=PCT_TOL)
    assert csv.loc["mean", "null_abs_gt_0_9"] == 0.0
    assert csv.loc["max", "spearman_r"] == pytest.approx(0.1233, abs=RHO_TOL)
    assert csv.loc["mean", "spearman_r"] == pytest.approx(0.1176, abs=RHO_TOL)
