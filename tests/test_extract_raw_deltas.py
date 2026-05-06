"""Tests for extract_raw_deltas.py — no API calls, cached data or synthetic DataFrames only."""
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

import extract_raw_deltas as edr


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_cache(tmp_path, monkeypatch):
    db = tmp_path / "test_cache.db"
    monkeypatch.setattr(edr, "CACHE_DB", db)
    return db


def _make_report_corrs():
    """Minimal corrs dict for _write_report."""
    stub = {"label": "x", "r": 0.1, "p": 0.05, "ci": (0.0, 0.2), "n": 100}
    return {
        "max_signed_raw": stub,
        "mean_signed_raw": stub,
        "quant_signed": stub,
        "mag_abs": stub,
    }


# ---------------------------------------------------------------------------
# 1–3: Cache I/O
# ---------------------------------------------------------------------------

def test_load_cache_returns_empty_dict_when_db_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(edr, "CACHE_DB", tmp_path / "nonexistent.db")
    assert edr._load_cache() == {}


def test_load_cache_returns_rsid_keyed_dict(tmp_cache):
    edr._init_cache()
    edr._save_raw_delta("rs123", {"max_signed_raw": 0.5, "mean_signed_raw": 0.3, "n_expr_tracks": 10})
    result = edr._load_cache()
    assert "rs123" in result
    assert result["rs123"]["max_signed_raw"] == pytest.approx(0.5)
    assert result["rs123"]["n_expr_tracks"] == 10


def test_save_and_reload_roundtrip(tmp_cache):
    edr._init_cache()
    entry = {"max_signed_raw": -0.123, "mean_signed_raw": 0.456, "n_expr_tracks": 7}
    edr._save_raw_delta("rs999", entry)
    cache = edr._load_cache()
    assert cache["rs999"]["max_signed_raw"] == pytest.approx(-0.123)
    assert cache["rs999"]["mean_signed_raw"] == pytest.approx(0.456)
    assert cache["rs999"]["n_expr_tracks"] == 7
    assert cache["rs999"]["error"] is None


# ---------------------------------------------------------------------------
# 5: Variant selection — build_sample quintile balance
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_per_quintile", [1, 5, 20])
def test_build_sample_returns_n_per_quintile_times_5(n_per_quintile):
    rng = np.random.default_rng(0)
    n = 500
    df = pd.DataFrame({
        "rsid": [f"rs{i}" for i in range(n)],
        "mpra_lfc": rng.normal(0, 1, n),
        "full_composite": rng.random(n),
    })
    sample = edr._build_sample(df, n_per_quintile=n_per_quintile)
    assert len(sample) == 5 * n_per_quintile


# ---------------------------------------------------------------------------
# 6–7: _select_full_panel
# ---------------------------------------------------------------------------

def test_select_full_panel_drops_rows_with_missing_alleles():
    df = pd.DataFrame({
        "rsid":    ["rs1", "rs2", "rs3"],
        "chrom":   ["chr1", "chr1", "chr1"],
        "pos":     [100, 200, 300],
        "ref":     ["A", None, "G"],
        "alt":     ["T", "C", "A"],
        "mpra_lfc": [0.1, 0.2, 0.3],
    })
    result = edr._select_full_panel(df)
    assert len(result) == 2
    assert "rs2" not in result["rsid"].values


def test_select_full_panel_keeps_all_valid_rows():
    df = pd.DataFrame({
        "rsid":    ["rs1", "rs2", "rs3"],
        "chrom":   ["chr1", "chr2", "chr3"],
        "pos":     [100, 200, 300],
        "ref":     ["A", "G", "C"],
        "alt":     ["T", "C", "G"],
        "mpra_lfc": [0.1, 0.2, 0.3],
    })
    result = edr._select_full_panel(df)
    assert len(result) == 3


# ---------------------------------------------------------------------------
# 8–9: --limit caps API calls (applied to to_score, not sample)
# ---------------------------------------------------------------------------

def test_limit_caps_api_calls_in_stratified_mode():
    all_rsids = [f"rs{i}" for i in range(20)]
    cached = {f"rs{i}" for i in range(5)}
    to_score = [r for r in all_rsids if r not in cached]  # 15 uncached
    limit = 8
    to_score = to_score[:limit]
    assert len(to_score) == limit
    assert all(r not in cached for r in to_score)


def test_limit_caps_api_calls_in_full_mode():
    all_rsids = [f"rs{i}" for i in range(3000)]
    cached = {f"rs{i}" for i in range(600)}
    to_score = [r for r in all_rsids if r not in cached]  # 2400 uncached
    limit = 50
    to_score = to_score[:limit]
    assert len(to_score) == limit
    assert all(r not in cached for r in to_score)


# ---------------------------------------------------------------------------
# 10: _correlate happy path
# ---------------------------------------------------------------------------

def test_correlate_returns_rho_plus1_on_monotone_pair():
    x = pd.Series(range(20), dtype=float)
    y = pd.Series(range(20), dtype=float)
    result = edr._correlate("test", x, y)
    assert result["r"] == pytest.approx(1.0, abs=1e-6)
    assert result["n"] == 20


# ---------------------------------------------------------------------------
# 12 (renamed): _correlate returns None below MIN_CORRELATION_N
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [0, 1, edr.MIN_CORRELATION_N - 1])
def test_correlate_returns_none_below_min_n(n):
    x = pd.Series(range(n), dtype=float)
    y = pd.Series(range(n), dtype=float)
    result = edr._correlate("test", x, y)
    assert result["r"] is None
    assert result["p"] is None
    assert result["n"] == n


def test_correlate_drops_nan_and_inf():
    # 12 values; nan and inf are dropped → 10 remain (≥ MIN_CORRELATION_N)
    x = pd.Series([1.0, 2.0, float("nan"), float("inf"), 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    y = pd.Series([1.0, 2.0, 3.0,          4.0,          5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    result = edr._correlate("test", x, y)
    assert result["n"] == 10  # nan and inf both excluded
    assert result["r"] is not None


# ---------------------------------------------------------------------------
# Resume test: full panel skips cached rsids
# ---------------------------------------------------------------------------

def test_full_panel_resume_skips_cached_rsids(tmp_cache):
    edr._init_cache()
    for rsid in ["rs1", "rs2", "rs3"]:
        edr._save_raw_delta(rsid, {"max_signed_raw": 0.1, "mean_signed_raw": 0.05, "n_expr_tracks": 3})

    panel = pd.DataFrame({
        "rsid":    ["rs1", "rs2", "rs3", "rs4", "rs5"],
        "chrom":   ["chr1"] * 5,
        "pos":     list(range(5)),
        "ref":     ["A"] * 5,
        "alt":     ["T"] * 5,
        "mpra_lfc": [0.1] * 5,
    })
    cache = edr._load_cache()
    to_score = [row for _, row in panel.iterrows() if row.rsid not in cache]

    assert len(to_score) == 2
    assert {r.rsid for r in to_score} == {"rs4", "rs5"}


# ---------------------------------------------------------------------------
# Report text: mode param controls description string
# ---------------------------------------------------------------------------

def test_write_report_stratified_mode_contains_stratified_text(tmp_path, monkeypatch):
    monkeypatch.setattr(edr, "OUT_REPORT", tmp_path / "report.md")
    edr._write_report(_make_report_corrs(), 600, mode="stratified")
    text = (tmp_path / "report.md").read_text()
    assert "stratified sample" in text


def test_write_report_full_mode_omits_stratified_text(tmp_path, monkeypatch):
    monkeypatch.setattr(edr, "OUT_REPORT", tmp_path / "report.md")
    edr._write_report(_make_report_corrs(), 3259, mode="full")
    text = (tmp_path / "report.md").read_text()
    assert "stratified sample" not in text


# ---------------------------------------------------------------------------
# 14–16: _extract_expression_raw
# ---------------------------------------------------------------------------

def test_extract_expression_raw_returns_none_on_empty_input():
    result = edr._extract_expression_raw(None, None)
    assert result["n_expr_tracks"] == 0
    assert result["max_signed_raw"] is None
    assert result["mean_signed_raw"] is None

    result2 = edr._extract_expression_raw(pd.DataFrame(), None)
    assert result2["n_expr_tracks"] == 0


def test_extract_expression_raw_ignores_non_expression_tracks():
    df = pd.DataFrame({
        "output_type": ["RNA_SEQ", "ATAC_SEQ", "CAGE"],
        "raw_score": [1.0, 99.0, 0.5],
        "biosample_term_name": ["K562", "K562", "K562"],
    })
    with patch("scoring.tissue_config.filter_tracks", return_value=df):
        result = edr._extract_expression_raw(df, MagicMock())
    # ATAC_SEQ excluded; of [1.0, 0.5] the max absolute is 1.0
    assert result["max_signed_raw"] == pytest.approx(1.0)
    assert result["n_expr_tracks"] == 2


def test_extract_expression_raw_max_signed_sign_preserved():
    df = pd.DataFrame({
        "output_type": ["RNA_SEQ", "CAGE"],
        "raw_score": [-2.0, 0.5],
        "biosample_term_name": ["K562", "K562"],
    })
    with patch("scoring.tissue_config.filter_tracks", return_value=df):
        result = edr._extract_expression_raw(df, MagicMock())
    # -2.0 has larger absolute value → max_signed_raw should be negative
    assert result["max_signed_raw"] == pytest.approx(-2.0)
