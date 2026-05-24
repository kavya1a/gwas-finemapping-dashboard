"""Order-statistic probe: max vs mean aggregation on Tewhey raw deltas.

The README's `What would change the interpretation` section names mean/median
aggregation as a direct mechanism probe: if max-over-tracks is what produces
saturation, swapping to mean should drop saturation on the same panel.

The Tewhey raw delta cache already stores both per-variant statistics
(see extract_raw_deltas.py: max_signed_raw and mean_signed_raw), so the
analytical portion of the test runs without any new API calls.

Outputs:
  mean_aggregation_comparison.csv
  figures/mean_vs_max_aggregation.png

The full version of this probe — re-scoring the matched-calibration null
for mean_signed_raw so the full published-quantile pipeline can be re-run
under mean aggregation — is a separate API run (see context.md, thread #1b).
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
FIG_DIR        = DIR / "figures"
OUT_CSV        = DIR / "mean_aggregation_comparison.csv"
OUT_FIG        = FIG_DIR / "mean_vs_max_aggregation.png"

N_BOOTSTRAP = 1000
SEED        = 42

C_MAX  = "#1a9850"   # green — max (current pipeline)
C_MEAN = "#762a83"   # purple — mean (alternative aggregation)


def load_tewhey() -> pd.DataFrame:
    parq = pd.read_parquet(TEWHEY_PARQUET)
    conn = sqlite3.connect(TEWHEY_CACHE)
    cache = pd.read_sql_query(
        "SELECT rsid, max_signed_raw, mean_signed_raw, n_expr_tracks, "
        "error AS cache_error FROM raw_deltas",
        conn,
    )
    conn.close()
    merged = parq.merge(cache, on="rsid", how="inner")
    valid = merged[
        merged["cache_error"].isna()
        & merged["max_signed_raw"].notna()
        & merged["mean_signed_raw"].notna()
        & merged["mpra_lfc"].notna()
    ].copy()
    return valid


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


def saturation_rates(x: np.ndarray) -> dict:
    a = np.abs(x[np.isfinite(x)])
    return {
        "abs_gt_0_5": float((a > 0.5).mean()),
        "abs_gt_0_9": float((a > 0.9).mean()),
        "iqr":        float(np.percentile(a, 75) - np.percentile(a, 25)),
        "median_abs": float(np.median(a)),
    }


def _annotate(ax, label_lines: list[str]) -> None:
    ax.text(0.97, 0.97, "\n".join(label_lines),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#cccccc", alpha=0.9))


def plot(df: pd.DataFrame, r_max: dict, r_mean: dict,
         s_max: dict, s_mean: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(13.5, 9.5))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.26)

    x_max  = df["max_signed_raw"].to_numpy(float)
    x_mean = df["mean_signed_raw"].to_numpy(float)
    y      = df["mpra_lfc"].to_numpy(float)

    # ---- Panel A: max distribution -----------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    lim = float(np.max(np.abs(x_max))) * 1.05
    ax.hist(x_max, bins=np.linspace(-lim, lim, 60),
            color=C_MAX, edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(+0.5, color="black", ls="--", lw=0.5)
    ax.axvline(-0.5, color="black", ls="--", lw=0.5)
    ax.set_xlabel("Max signed raw expression Δ (current pipeline)")
    ax.set_ylabel("# Tewhey variants")
    ax.set_title("(A) Max-over-tracks aggregation — what AlphaGenome publishes downstream",
                 fontsize=10.5)
    _annotate(ax, [
        f"|·|>0.5: {s_max['abs_gt_0_5']:.2%}",
        f"|·|>0.9: {s_max['abs_gt_0_9']:.2%}",
        f"median |·|: {s_max['median_abs']:.4f}",
        f"IQR(|·|):   {s_max['iqr']:.4f}",
        f"n = {r_max['n']:,}",
    ])
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel B: mean distribution ----------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    lim = float(np.max(np.abs(x_mean))) * 1.05
    ax.hist(x_mean, bins=np.linspace(-lim, lim, 60),
            color=C_MEAN, edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(+0.5, color="black", ls="--", lw=0.5)
    ax.axvline(-0.5, color="black", ls="--", lw=0.5)
    ax.set_xlabel("Mean signed raw expression Δ (alternative aggregation)")
    ax.set_ylabel("# Tewhey variants")
    ax.set_title("(B) Mean-over-tracks aggregation — no order-statistic inflation",
                 fontsize=10.5)
    _annotate(ax, [
        f"|·|>0.5: {s_mean['abs_gt_0_5']:.2%}",
        f"|·|>0.9: {s_mean['abs_gt_0_9']:.2%}",
        f"median |·|: {s_mean['median_abs']:.4f}",
        f"IQR(|·|):   {s_mean['iqr']:.4f}",
        f"n = {r_mean['n']:,}",
    ])
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # ---- Panel C: max vs LFC scatter ---------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(x_max, y, s=6, c=C_MAX, alpha=0.35,
               edgecolors="none", rasterized=True)
    ax.axhline(0, color="grey", ls=":", lw=0.6)
    ax.axvline(0, color="grey", ls=":", lw=0.6)
    ax.set_xlabel("Max signed raw Δ")
    ax.set_ylabel("MPRA log fold change")
    ax.set_title(f"(C) Max vs MPRA LFC — ρ = {r_max['r']:+.4f}", fontsize=10.5)
    _annotate(ax, [
        f"ρ = {r_max['r']:+.4f}",
        f"95% CI: [{r_max['ci_lo']:+.4f}, {r_max['ci_hi']:+.4f}]",
        f"p = {r_max['p']:.2e}",
        f"n = {r_max['n']:,}",
    ])
    ax.grid(lw=0.3, alpha=0.5)

    # ---- Panel D: mean vs LFC scatter --------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(x_mean, y, s=6, c=C_MEAN, alpha=0.35,
               edgecolors="none", rasterized=True)
    ax.axhline(0, color="grey", ls=":", lw=0.6)
    ax.axvline(0, color="grey", ls=":", lw=0.6)
    ax.set_xlabel("Mean signed raw Δ")
    ax.set_ylabel("MPRA log fold change")
    ax.set_title(f"(D) Mean vs MPRA LFC — ρ = {r_mean['r']:+.4f}", fontsize=10.5)
    _annotate(ax, [
        f"ρ = {r_mean['r']:+.4f}",
        f"95% CI: [{r_mean['ci_lo']:+.4f}, {r_mean['ci_hi']:+.4f}]",
        f"p = {r_mean['p']:.2e}",
        f"n = {r_mean['n']:,}",
    ])
    ax.grid(lw=0.3, alpha=0.5)

    fig.suptitle(
        "Tewhey MPRA — max vs mean aggregation of K562 expression raw deltas",
        fontsize=12, fontweight="bold", y=0.995,
    )
    fig.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure → {OUT_FIG}")


def main() -> None:
    print("Loading Tewhey + raw delta cache...")
    tew = load_tewhey()
    print(f"  Tewhey valid (max + mean + lfc, no error): {len(tew):,}")
    print(f"  Track count per variant — median: {int(tew['n_expr_tracks'].median())}, "
          f"min: {int(tew['n_expr_tracks'].min())}, "
          f"max: {int(tew['n_expr_tracks'].max())}")

    y      = tew["mpra_lfc"].to_numpy(float)
    x_max  = tew["max_signed_raw"].to_numpy(float)
    x_mean = tew["mean_signed_raw"].to_numpy(float)

    print("\nComputing Spearman + bootstrap CIs (1,000 iters each)...")
    r_max  = spearman_with_bootstrap(x_max, y)
    r_mean = spearman_with_bootstrap(x_mean, y)

    s_max  = saturation_rates(x_max)
    s_mean = saturation_rates(x_mean)

    df_table = pd.DataFrame([
        {"aggregation":   "max-over-tracks (current)",
         "spearman_r":    r_max["r"],
         "ci_lo":         r_max["ci_lo"],
         "ci_hi":         r_max["ci_hi"],
         "p_value":       r_max["p"],
         "abs_gt_0_5":    s_max["abs_gt_0_5"],
         "abs_gt_0_9":    s_max["abs_gt_0_9"],
         "median_abs":    s_max["median_abs"],
         "n":             r_max["n"]},
        {"aggregation":   "mean-over-tracks",
         "spearman_r":    r_mean["r"],
         "ci_lo":         r_mean["ci_lo"],
         "ci_hi":         r_mean["ci_hi"],
         "p_value":       r_mean["p"],
         "abs_gt_0_5":    s_mean["abs_gt_0_5"],
         "abs_gt_0_9":    s_mean["abs_gt_0_9"],
         "median_abs":    s_mean["median_abs"],
         "n":             r_mean["n"]},
    ])
    df_table.to_csv(OUT_CSV, index=False)
    print(f"Table → {OUT_CSV}")

    print("\n=== Aggregation comparison on Tewhey raw deltas ===")
    hdr = (f"{'Aggregation':<28}  {'ρ':>9}  {'95% CI':>22}  "
           f"{'p':>10}  {'|·|>0.5':>8}  {'|·|>0.9':>8}  {'n':>7}")
    print(hdr)
    print("-" * len(hdr))
    for label, r, s in [("max-over-tracks", r_max,  s_max),
                        ("mean-over-tracks", r_mean, s_mean)]:
        ci = f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]"
        print(f"{label:<28}  {r['r']:+.4f}  {ci:>22}  {r['p']:.2e}  "
              f"{s['abs_gt_0_5']:>7.2%}  {s['abs_gt_0_9']:>7.2%}  {r['n']:>7,}")

    plot(tew, r_max, r_mean, s_max, s_mean)
    print("\nDone.")


if __name__ == "__main__":
    main()
