"""Three-recipe matched-calibration comparison: max / mean / median.

For each aggregation recipe r ∈ {max, mean, median}:

  1. Build a null distribution from raw_{r}_signed_delta on the matched
     common-variant pool (matched_calibration_cache.db).
  2. Re-quantile each Tewhey variant's raw_{r}_signed_delta against that
     null, signed and mapped to [-1, +1].
  3. Compute Spearman vs Tewhey mpra_lfc + 95% CI by bootstrap.
  4. Report tail saturation on both the null and the Tewhey panel.

If the calibration-statistic mismatch / order-statistic explanation is
right, max-aggregation should show high tail saturation under both null
and Tewhey re-quantiling that mean/median should mostly eliminate, while
all three should give similar Spearman vs MPRA LFC.

Outputs:
  matched_recipes_comparison.csv
  figures/matched_recipes_comparison.png
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
MATCHED_CACHE  = DIR / "matched_calibration_cache.db"
FIG_DIR        = DIR / "figures"
OUT_CSV        = DIR / "matched_recipes_comparison.csv"
OUT_FIG        = FIG_DIR / "matched_recipes_comparison.png"

N_BOOTSTRAP = 1000
SEED        = 42

# Mapping: recipe → (tewhey cache column, matched cache column, display)
RECIPES = [
    ("max",    "max_signed_raw",    "raw_max_signed_delta",    "#1a9850"),
    ("mean",   "mean_signed_raw",   "raw_mean_signed_delta",   "#762a83"),
    ("median", None,                "raw_median_signed_delta", "#d6604d"),  # tewhey median not cached
]


def load_tewhey() -> pd.DataFrame:
    parq = pd.read_parquet(TEWHEY_PARQUET)
    conn = sqlite3.connect(TEWHEY_CACHE)
    cache = pd.read_sql_query(
        "SELECT rsid, max_signed_raw, mean_signed_raw, "
        "error AS cache_error FROM raw_deltas",
        conn,
    )
    conn.close()
    return parq.merge(cache, on="rsid", how="inner")


def load_matched_nulls() -> dict[str, np.ndarray]:
    conn = sqlite3.connect(MATCHED_CACHE)
    df = pd.read_sql_query(
        "SELECT raw_max_signed_delta, raw_mean_signed_delta, "
        "raw_median_signed_delta, error FROM scores",
        conn,
    )
    conn.close()
    df = df[df["error"].isna()]
    return {
        "max":    np.sort(df["raw_max_signed_delta"].dropna().to_numpy(float)),
        "mean":   np.sort(df["raw_mean_signed_delta"].dropna().to_numpy(float)),
        "median": np.sort(df["raw_median_signed_delta"].dropna().to_numpy(float)),
    }


def matched_quantile_signed(x: np.ndarray, null_sorted: np.ndarray) -> np.ndarray:
    """Signed ECDF rank: 2 * (n_lt + 0.5*n_eq)/n_null - 1, mapped to [-1, +1]."""
    n_null = len(null_sorted)
    if n_null == 0:
        return np.full_like(x, np.nan, dtype=float)
    n_lt = np.searchsorted(null_sorted, x, side="left")
    n_le = np.searchsorted(null_sorted, x, side="right")
    n_eq = n_le - n_lt
    q_unsigned = (n_lt + 0.5 * n_eq) / n_null
    return 2.0 * q_unsigned - 1.0


def spearman_with_bootstrap(x: np.ndarray, y: np.ndarray,
                            n_boot: int = N_BOOTSTRAP, seed: int = SEED) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    xv, yv = x[mask], y[mask]
    n = int(len(xv))
    r, p = st.spearmanr(xv, yv) if n >= 3 else (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rb, _ = st.spearmanr(xv[idx], yv[idx])
        boots[i] = rb
    return {
        "r": float(r), "p": float(p), "n": n,
        "ci_lo": float(np.percentile(boots, 2.5)),
        "ci_hi": float(np.percentile(boots, 97.5)),
    }


def sat_stats(x: np.ndarray) -> dict:
    a = np.abs(x[np.isfinite(x)])
    return {
        "n":          int(len(a)),
        "abs_gt_0_5": float((a > 0.5).mean()) if len(a) else float("nan"),
        "abs_gt_0_9": float((a > 0.9).mean()) if len(a) else float("nan"),
        "exact_1":    float(((np.abs(x) == 1.0)).mean()) if len(a) else float("nan"),
    }


def _annotate(ax, lines: list[str]) -> None:
    ax.text(0.97, 0.97, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#cccccc", alpha=0.9))


def plot(tew: pd.DataFrame, nulls: dict, results: dict, sat_null: dict,
        sat_tewhey_q: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.30)

    for col, (recipe, _tew_col, _mat_col, color) in enumerate(RECIPES):
        # Row 1: null distribution
        ax = fig.add_subplot(gs[0, col])
        n = nulls[recipe]
        lim = float(np.percentile(np.abs(n), 99.5)) * 1.1 if len(n) else 1.0
        ax.hist(n, bins=np.linspace(-lim, lim, 60),
                color=color, edgecolor="white", alpha=0.85)
        ax.axvline(0, color="grey", ls=":", lw=0.7)
        ax.set_title(f"({chr(65+col)}1) Null — common variants, {recipe} aggregation",
                     fontsize=10)
        ax.set_xlabel(f"raw_{recipe}_signed_delta")
        ax.set_ylabel("# variants")
        sn = sat_null[recipe]
        _annotate(ax, [f"n = {sn['n']:,}",
                       f"|·|>0.5: {sn['abs_gt_0_5']:.2%}",
                       f"|·|>0.9: {sn['abs_gt_0_9']:.2%}"])
        ax.grid(axis="y", lw=0.3, alpha=0.5)

        # Row 2: Tewhey re-quantiled against null
        ax = fig.add_subplot(gs[1, col])
        q = results[recipe]["quantile"]
        if q is not None:
            ax.hist(q, bins=np.linspace(-1, 1, 60),
                    color=color, edgecolor="white", alpha=0.85)
            ax.axvline(+0.9, color="black", ls="--", lw=0.5)
            ax.axvline(-0.9, color="black", ls="--", lw=0.5)
            ax.set_xlim(-1.05, 1.05)
            sq = sat_tewhey_q[recipe]
            _annotate(ax, [f"n = {sq['n']:,}",
                           f"|q|>0.5: {sq['abs_gt_0_5']:.2%}",
                           f"|q|>0.9: {sq['abs_gt_0_9']:.2%}",
                           f"q=±1: {sq['exact_1']:.2%}"])
        else:
            ax.text(0.5, 0.5, "Tewhey median\nnot cached",
                    ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"({chr(65+col)}2) Tewhey quantile under {recipe}-aggregation null",
                     fontsize=10)
        ax.set_xlabel(f"matched-{recipe} quantile (signed)")
        ax.set_ylabel("# variants")
        ax.grid(axis="y", lw=0.3, alpha=0.5)

        # Row 3: scatter vs LFC
        ax = fig.add_subplot(gs[2, col])
        r = results[recipe]
        if r["quantile"] is not None:
            ax.scatter(r["quantile"], r["lfc"], s=6, c=color, alpha=0.35,
                       edgecolors="none", rasterized=True)
            ax.axhline(0, color="grey", ls=":", lw=0.6)
            ax.axvline(0, color="grey", ls=":", lw=0.6)
            sp = r["spearman"]
            _annotate(ax, [f"ρ = {sp['r']:+.4f}",
                           f"CI: [{sp['ci_lo']:+.4f}, {sp['ci_hi']:+.4f}]",
                           f"p = {sp['p']:.2e}",
                           f"n = {sp['n']:,}"])
        else:
            ax.text(0.5, 0.5, "n/a", ha="center", va="center",
                    transform=ax.transAxes)
        ax.set_title(f"({chr(65+col)}3) MPRA LFC vs matched-{recipe} quantile",
                     fontsize=10)
        ax.set_xlabel(f"matched-{recipe} quantile")
        ax.set_ylabel("MPRA log fold change")
        ax.grid(lw=0.3, alpha=0.5)

    fig.suptitle("Matched-calibration recipes compared — max / mean / median aggregation",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure → {OUT_FIG}")


def main() -> None:
    print("Loading matched-calibration nulls (max / mean / median)...")
    nulls = load_matched_nulls()
    for r, x in nulls.items():
        print(f"  null[{r}]:  n={len(x):,}, "
              f"|·|>0.5={float((np.abs(x)>0.5).mean()):.3%}, "
              f"|·|>0.9={float((np.abs(x)>0.9).mean()):.3%}")

    print("\nLoading Tewhey...")
    tew = load_tewhey()
    tew_valid = tew[tew["cache_error"].isna() & tew["mpra_lfc"].notna()].copy()
    print(f"  Tewhey rows with valid mpra_lfc + no cache error: {len(tew_valid):,}")

    results = {}
    sat_null = {r: sat_stats(nulls[r]) for r, *_ in RECIPES}
    sat_tewhey_q = {}

    for recipe, tew_col, _mat_col, _ in RECIPES:
        if tew_col is None:
            # Tewhey doesn't have median cached; only show null + saturation
            results[recipe] = {"quantile": None, "lfc": None,
                               "spearman": {"r": float("nan"), "p": float("nan"),
                                            "n": 0, "ci_lo": float("nan"),
                                            "ci_hi": float("nan")}}
            sat_tewhey_q[recipe] = {"n": 0, "abs_gt_0_5": float("nan"),
                                     "abs_gt_0_9": float("nan"),
                                     "exact_1": float("nan")}
            continue

        sub = tew_valid[tew_valid[tew_col].notna()].copy()
        raw = sub[tew_col].to_numpy(float)
        lfc = sub["mpra_lfc"].to_numpy(float)
        q = matched_quantile_signed(raw, nulls[recipe])
        sp = spearman_with_bootstrap(q, lfc)
        sat_tewhey_q[recipe] = sat_stats(q)
        results[recipe] = {"quantile": q, "lfc": lfc, "spearman": sp}

    # Console table
    print("\n=== Three-recipe matched-calibration comparison ===")
    hdr = (f"{'Recipe':<8} | {'Null n':>6} | {'Null|·|>0.5':>11} | {'Null|·|>0.9':>11} | "
           f"{'Tew n':>6} | {'Tew q|·|>0.5':>13} | {'Tew q|·|>0.9':>13} | "
           f"{'ρ':>9} | {'p':>10}")
    print(hdr); print("-" * len(hdr))
    rows_for_csv = []
    for recipe, _tc, _mc, _ in RECIPES:
        sn = sat_null[recipe]
        sq = sat_tewhey_q[recipe]
        sp = results[recipe]["spearman"]
        rho_s = f"{sp['r']:+.4f}" if not np.isnan(sp["r"]) else "  n/a  "
        p_s   = f"{sp['p']:.2e}"  if not np.isnan(sp["p"]) else "  n/a    "
        print(f"{recipe:<8} | {sn['n']:>6,} | "
              f"{sn['abs_gt_0_5']:>10.2%} | {sn['abs_gt_0_9']:>10.2%} | "
              f"{sq['n']:>6,} | "
              f"{sq['abs_gt_0_5']:>12.2%} | {sq['abs_gt_0_9']:>12.2%} | "
              f"{rho_s:>9} | {p_s:>10}")
        rows_for_csv.append({
            "recipe":              recipe,
            "null_n":              sn["n"],
            "null_abs_gt_0_5":     sn["abs_gt_0_5"],
            "null_abs_gt_0_9":     sn["abs_gt_0_9"],
            "tewhey_n":            sq["n"],
            "tewhey_q_abs_gt_0_5": sq["abs_gt_0_5"],
            "tewhey_q_abs_gt_0_9": sq["abs_gt_0_9"],
            "tewhey_q_exact_1":    sq["exact_1"],
            "spearman_r":          sp["r"],
            "spearman_p":          sp["p"],
            "spearman_ci_lo":      sp["ci_lo"],
            "spearman_ci_hi":      sp["ci_hi"],
        })
    pd.DataFrame(rows_for_csv).to_csv(OUT_CSV, index=False)
    print(f"\nTable → {OUT_CSV}")

    plot(tew_valid, nulls, results, sat_null, sat_tewhey_q)
    print("\nDone.")


if __name__ == "__main__":
    main()
