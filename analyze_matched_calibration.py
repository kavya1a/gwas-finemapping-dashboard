"""Component 3: matched-calibration analysis on the Tewhey panel.

For each Tewhey variant in tewhey_raw_delta_cache.db, compute its quantile
within matched_calibration_null.parquet (max-over-tracks K562 expression
delta on common autosomal variants). Then build a four-row Spearman table:

  1. Original quantile (single-track calibration + max-over-tracks applied)
  2. Matched-calibration quantile (max-over-tracks calibration + max applied)
  3. Phred empirical = -10 * log10(1 - matched_quantile + epsilon)
  4. Raw max signed delta (no normalization)

vs Tewhey mpra_lfc, with bootstrap 95% CIs.

Sanity check: rows 2 and 3 must have identical Spearman (phred is a
strict monotone transform of matched_quantile).

Outputs:
  matched_calibration_comparison.csv
  figures/three_way_comparison.png
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
FIG_DIR        = DIR / "figures"
OUT_CSV        = DIR / "matched_calibration_comparison.csv"
OUT_FIG        = FIG_DIR / "three_way_comparison.png"

N_BOOTSTRAP = 1000
SEED        = 42
EPSILON     = 1e-6  # phred = -10 * log10(1 - q + eps); avoids inf at q=+1


# ---------------------------------------------------------------------------
# Load + merge
# ---------------------------------------------------------------------------

def load_tewhey() -> pd.DataFrame:
    parq = pd.read_parquet(TEWHEY_PARQUET)
    conn = sqlite3.connect(TEWHEY_CACHE)
    cache = pd.read_sql_query(
        "SELECT rsid, max_signed_raw, mean_signed_raw, "
        "error AS cache_error FROM raw_deltas",
        conn,
    )
    conn.close()
    merged = parq.merge(cache, on="rsid", how="inner")
    # successful raw delta + valid mpra_lfc
    valid = merged[
        merged["cache_error"].isna()
        & merged["max_signed_raw"].notna()
        & merged["mpra_lfc"].notna()
    ].copy()
    return valid


def load_null() -> np.ndarray:
    nul = pd.read_parquet(NULL_PARQUET)
    return np.sort(nul["raw_max_signed_delta"].dropna().to_numpy(float))


# ---------------------------------------------------------------------------
# Matched quantile (signed) + phred
# ---------------------------------------------------------------------------

def matched_quantile_signed(x: np.ndarray, null_sorted: np.ndarray) -> np.ndarray:
    """ECDF of x within null, mid-rank for ties, mapped linearly to [-1, +1].

    Mirrors AlphaGenome's signed-quantile convention:
      q_unsigned = (n_lt + 0.5 * n_eq) / n_null     in (0, 1)
      q_signed   = 2 * q_unsigned - 1               in (-1, +1)
    """
    n_null = len(null_sorted)
    n_lt = np.searchsorted(null_sorted, x, side="left")
    n_le = np.searchsorted(null_sorted, x, side="right")
    n_eq = n_le - n_lt
    q_unsigned = (n_lt + 0.5 * n_eq) / n_null
    return 2.0 * q_unsigned - 1.0


def phred_empirical(q_signed: np.ndarray) -> np.ndarray:
    """phred = -10 * log10(1 - q + eps).

    Strictly monotonic on q ∈ (-1, +1). Negative phred for q < 0 is
    mathematically meaningful for ranking (preserves Spearman vs signed
    measurements like mpra_lfc) even if non-traditional.
    """
    return -10.0 * np.log10(1.0 - q_signed + EPSILON)


# ---------------------------------------------------------------------------
# Spearman with bootstrap
# ---------------------------------------------------------------------------

def spearman_with_bootstrap(x: np.ndarray, y: np.ndarray,
                            n_boot: int = N_BOOTSTRAP, seed: int = SEED) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    xv, yv = x[mask], y[mask]
    n = int(len(xv))
    r, p = st.spearmanr(xv, yv)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rb, _ = st.spearmanr(xv[idx], yv[idx])
        boots[i] = rb
    ci_lo = float(np.percentile(boots, 2.5))
    ci_hi = float(np.percentile(boots, 97.5))
    return {"r": float(r), "p": float(p), "ci_lo": ci_lo, "ci_hi": ci_hi, "n": n}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

C_ORIG    = "#d6604d"   # red — saturated, problematic
C_MATCHED = "#4575b4"   # blue — matched, hopefully cleaner
C_RAW     = "#1a9850"   # green — raw delta, continuous
C_PHRED   = "#762a83"   # purple — phred (monotone of matched)


def _annotate(ax, label_lines: list[str]) -> None:
    ax.text(0.97, 0.97, "\n".join(label_lines),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#cccccc", alpha=0.9))


def plot_three_way(df: pd.DataFrame, results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(13.5, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.24)

    # ---- Panel A: original quantile (single-track calib) ------------------
    ax = fig.add_subplot(gs[0, 0])
    x = df["expression_subscore"].dropna().to_numpy(float)
    ax.hist(x, bins=np.linspace(-1, 1, 60),
            color=C_ORIG, edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(+0.9, color="black", ls="--", lw=0.5)
    ax.axvline(-0.9, color="black", ls="--", lw=0.5)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Original quantile (single-track calibration, max-over-tracks applied)")
    ax.set_ylabel("# Tewhey variants")
    ax.set_title("(A) Original quantile — published single-track calibration",
                 fontsize=10.5)
    r1 = results["original"]
    _annotate(ax, [
        f"|q|>0.9 saturation: {(np.abs(x)>0.9).mean():.1%}",
        f"ρ vs MPRA LFC: {r1['r']:+.4f}",
        f"95% CI: [{r1['ci_lo']:+.4f}, {r1['ci_hi']:+.4f}]",
        f"n = {r1['n']:,}",
    ])
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel B: matched calibration quantile ---------------------------
    ax = fig.add_subplot(gs[0, 1])
    x = df["matched_quantile"].dropna().to_numpy(float)
    ax.hist(x, bins=np.linspace(-1, 1, 60),
            color=C_MATCHED, edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(+0.9, color="black", ls="--", lw=0.5)
    ax.axvline(-0.9, color="black", ls="--", lw=0.5)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Matched-calibration quantile (max-over-tracks null)")
    ax.set_ylabel("# Tewhey variants")
    ax.set_title("(B) Matched-calibration quantile — null built on max-over-tracks",
                 fontsize=10.5)
    r2 = results["matched"]
    _annotate(ax, [
        f"|q|>0.9 saturation: {(np.abs(x)>0.9).mean():.1%}",
        f"ρ vs MPRA LFC: {r2['r']:+.4f}",
        f"95% CI: [{r2['ci_lo']:+.4f}, {r2['ci_hi']:+.4f}]",
        f"n = {r2['n']:,}",
    ])
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel C: raw max signed delta -----------------------------------
    ax = fig.add_subplot(gs[1, 0])
    x = df["max_signed_raw"].dropna().to_numpy(float)
    lim = float(np.max(np.abs(x))) * 1.05
    ax.hist(x, bins=np.linspace(-lim, lim, 60),
            color=C_RAW, edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("Raw max signed expression delta (no normalization)")
    ax.set_ylabel("# Tewhey variants")
    ax.set_title("(C) Raw max signed Δ — continuous, unnormalized",
                 fontsize=10.5)
    r4 = results["raw"]
    _annotate(ax, [
        f"|Δ|>0.9 saturation: {(np.abs(x)>0.9).mean():.1%}",
        f"ρ vs MPRA LFC: {r4['r']:+.4f}",
        f"95% CI: [{r4['ci_lo']:+.4f}, {r4['ci_hi']:+.4f}]",
        f"n = {r4['n']:,}",
    ])
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel D: phred-vs-matched diagnostic ----------------------------
    ax = fig.add_subplot(gs[1, 1])
    q  = df["matched_quantile"].to_numpy(float)
    ph = df["phred_empirical"].to_numpy(float)
    ax.scatter(q, ph, s=8, c=C_PHRED, alpha=0.4, edgecolors="none",
               rasterized=True)
    qq = np.linspace(q.min(), q.max(), 400)
    ax.plot(qq, -10 * np.log10(1 - qq + EPSILON),
            color="black", ls="--", lw=1.0,
            label="-10·log₁₀(1 − q + ε)")
    r3 = results["phred"]
    _annotate(ax, [
        f"ρ(matched) = {r2['r']:+.6f}",
        f"ρ(phred)   = {r3['r']:+.6f}",
        f"|Δρ|       = {abs(r2['r'] - r3['r']):.2e}",
        f"identical? {'YES' if abs(r2['r'] - r3['r']) < 1e-6 else 'NO'}",
    ])
    ax.set_xlabel("matched_quantile (signed, in [-1, +1])")
    ax.set_ylabel("phred_empirical")
    ax.set_title("(D) Diagnostic — phred is monotone in matched quantile",
                 fontsize=10.5)
    ax.legend(loc="upper left", fontsize=8.5)
    ax.grid(lw=0.3, alpha=0.5)

    fig.suptitle(
        "Tewhey MPRA — three predictors of mpra_lfc, plus phred diagnostic",
        fontsize=12, fontweight="bold", y=0.995,
    )
    fig.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure → {OUT_FIG}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading Tewhey + matched null...")
    tew = load_tewhey()
    null_sorted = load_null()
    print(f"  Tewhey valid (raw + lfc + no error): {len(tew):,}")
    print(f"  Null sorted size:                    {len(null_sorted):,}")

    # Compute matched_quantile + phred for each Tewhey variant
    raw = tew["max_signed_raw"].to_numpy(float)
    tew["matched_quantile"] = matched_quantile_signed(raw, null_sorted)
    tew["phred_empirical"]  = phred_empirical(tew["matched_quantile"].to_numpy(float))

    # Restrict to rows where ALL FOUR predictors are valid for an
    # apples-to-apples comparison (some variants have NaN expression_subscore
    # in the parquet but a valid raw delta, or vice versa).
    cols_all4 = ["expression_subscore", "matched_quantile",
                 "phred_empirical", "max_signed_raw"]
    paired = tew.dropna(subset=cols_all4 + ["mpra_lfc"]).copy()
    print(f"  Variants with all four predictors + lfc: {len(paired):,}")

    y = paired["mpra_lfc"].to_numpy(float)

    print("\nComputing Spearman + bootstrap CIs (1,000 iters each)...")
    results = {
        "original": spearman_with_bootstrap(paired["expression_subscore"].to_numpy(float), y),
        "matched":  spearman_with_bootstrap(paired["matched_quantile"].to_numpy(float), y),
        "phred":    spearman_with_bootstrap(paired["phred_empirical"].to_numpy(float), y),
        "raw":      spearman_with_bootstrap(paired["max_signed_raw"].to_numpy(float), y),
    }

    # CSV
    rows = [
        ("Original quantile (single-track calib, max applied)", "original"),
        ("Matched-calibration quantile (max-track calib)",       "matched"),
        ("Phred empirical (monotone transform of #2)",            "phred"),
        ("Raw max signed delta (no normalization)",               "raw"),
    ]
    df_table = pd.DataFrame([
        {"predictor": label,
         "spearman_r": results[k]["r"],
         "ci_lo": results[k]["ci_lo"],
         "ci_hi": results[k]["ci_hi"],
         "p_value": results[k]["p"],
         "n": results[k]["n"]}
        for label, k in rows
    ])
    df_table.to_csv(OUT_CSV, index=False)
    print(f"Table → {OUT_CSV}")

    # Print
    print("\n=== Four-row Spearman comparison (Tewhey mpra_lfc) ===")
    hdr = f"{'Predictor':<55}  {'ρ':>9}  {'95% CI':>22}  {'p':>10}  {'n':>7}"
    print(hdr)
    print("-" * len(hdr))
    for label, k in rows:
        r = results[k]
        ci = f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]"
        print(f"{label:<55}  {r['r']:+.4f}  {ci:>22}  {r['p']:.2e}  {r['n']:>7,}")

    # Sanity check: matched and phred must have identical Spearman
    diff = abs(results["matched"]["r"] - results["phred"]["r"])
    print(f"\nSanity check (monotone transform preserves rank):")
    print(f"  ρ(matched) = {results['matched']['r']:+.8f}")
    print(f"  ρ(phred)   = {results['phred']['r']:+.8f}")
    print(f"  |Δρ|       = {diff:.3e}")
    if diff < 1e-6:
        print("  PASS: identical to numerical precision.")
    else:
        print("  FAIL: phred and matched_quantile Spearman diverge — investigate.")

    # Saturation snapshot
    print("\nSaturation snapshot (Tewhey panel):")
    for label, col in [("Original quantile",     "expression_subscore"),
                       ("Matched quantile",      "matched_quantile"),
                       ("Raw max signed Δ",      "max_signed_raw")]:
        x = paired[col].to_numpy(float)
        sat = (np.abs(x) > 0.9).mean()
        print(f"  {label:<22} |·|>0.9: {sat:.1%}")

    # Plot
    plot_three_way(paired, results)
    print("\nDone.")


if __name__ == "__main__":
    main()
