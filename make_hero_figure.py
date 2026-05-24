"""Hero figure for the README — two-panel before/after on the same Tewhey variants.

Panel A: published expression_subscore distribution (saturated at ±1).
Panel B: matched-calibration quantile (full range populated).

Designed to be screenshot-able. Single image, large fonts, two distributions
side-by-side with rho-vs-MPRA-LFC and saturation annotated.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st

DIR = Path(__file__).parent
TEWHEY_PARQUET = DIR / "tewhey_mpra.parquet"
TEWHEY_CACHE   = DIR / "tewhey_raw_delta_cache.db"
NULL_PARQUET   = DIR / "matched_calibration_null.parquet"
OUT_FIG        = DIR / "figures" / "hero_recovery.png"


def matched_quantile_signed(x: np.ndarray, null_sorted: np.ndarray) -> np.ndarray:
    n_null = len(null_sorted)
    n_lt = np.searchsorted(null_sorted, x, side="left")
    n_le = np.searchsorted(null_sorted, x, side="right")
    n_eq = n_le - n_lt
    return 2.0 * (n_lt + 0.5 * n_eq) / n_null - 1.0


def main() -> None:
    parq = pd.read_parquet(TEWHEY_PARQUET)
    conn = sqlite3.connect(TEWHEY_CACHE)
    cache = pd.read_sql_query(
        "SELECT rsid, max_signed_raw, error AS cache_error FROM raw_deltas",
        conn,
    )
    conn.close()
    merged = parq.merge(cache, on="rsid", how="inner")

    null_sorted = np.sort(
        pd.read_parquet(NULL_PARQUET)["raw_max_signed_delta"].dropna().to_numpy(float)
    )

    df = merged.copy()
    df["matched_quantile"] = matched_quantile_signed(
        df["max_signed_raw"].to_numpy(float), null_sorted
    )

    valid = df.dropna(subset=["expression_subscore", "matched_quantile",
                              "max_signed_raw", "mpra_lfc"]).copy()
    n = len(valid)

    pub = valid["expression_subscore"].to_numpy(float)
    matched = valid["matched_quantile"].to_numpy(float)
    lfc = valid["mpra_lfc"].to_numpy(float)

    r_pub, p_pub = st.spearmanr(pub, lfc)
    r_mat, p_mat = st.spearmanr(matched, lfc)

    sat_pub = float((np.abs(pub) > 0.9).mean())
    sat_mat = float((np.abs(matched) > 0.9).mean())

    # ----- figure ----------------------------------------------------------
    plt.rcParams.update({
        "font.family":      "sans-serif",
        "font.size":        12,
        "axes.titlesize":   14,
        "axes.labelsize":   13,
        "xtick.labelsize":  11,
        "ytick.labelsize":  11,
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)

    # ---- Panel A: published -----------------------------------------------
    ax = axes[0]
    bins = np.linspace(-1, 1, 70)
    ax.hist(pub, bins=bins, color="#d6604d", edgecolor="white", alpha=0.92)
    ax.axvline(+0.9, color="#444444", ls="--", lw=1.0)
    ax.axvline(-0.9, color="#444444", ls="--", lw=1.0)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Published quantile-calibrated expression score")
    ax.set_ylabel("Number of Tewhey variants")
    ax.set_title("Published quantile (single-track calibration, max-aggregated)",
                 loc="left", pad=12)
    annot = (f"ρ vs MPRA LFC: {r_pub:+.3f}\n"
             f"p = {p_pub:.1e}\n"
             f"|score|>0.9: {sat_pub:.1%}\n"
             f"n = {n:,}")
    ax.text(0.04, 0.96, annot, transform=ax.transAxes,
            ha="left", va="top", fontsize=12, family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", fc="white",
                      ec="#666666", lw=0.8))
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel B: matched calibration -------------------------------------
    ax = axes[1]
    ax.hist(matched, bins=bins, color="#4575b4", edgecolor="white", alpha=0.92)
    ax.axvline(+0.9, color="#444444", ls="--", lw=1.0)
    ax.axvline(-0.9, color="#444444", ls="--", lw=1.0)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Matched-calibration quantile (max-over-tracks null)")
    ax.set_title("After matched-statistic re-calibration (same variants)",
                 loc="left", pad=12)
    annot = (f"ρ vs MPRA LFC: {r_mat:+.3f}\n"
             f"p = {p_mat:.1e}\n"
             f"|score|>0.9: {sat_mat:.1%}\n"
             f"n = {n:,}")
    ax.text(0.04, 0.96, annot, transform=ax.transAxes,
            ha="left", va="top", fontsize=12, family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", fc="white",
                      ec="#666666", lw=0.8))
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    fig.suptitle(
        "AlphaGenome on Tewhey 2016 MPRA — recovery via matched-statistic calibration",
        fontsize=15, fontweight="bold", y=1.00,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    OUT_FIG.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Hero figure → {OUT_FIG}")
    print(f"  Published:  ρ = {r_pub:+.4f}, saturation |·|>0.9 = {sat_pub:.1%}")
    print(f"  Matched:    ρ = {r_mat:+.4f}, saturation |·|>0.9 = {sat_mat:.1%}")
    print(f"  n = {n:,}")


if __name__ == "__main__":
    main()
