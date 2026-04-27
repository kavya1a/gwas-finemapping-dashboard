"""
Tewhey 2016 MPRA correlation figures.

Reads tewhey_mpra.parquet (output of tewhey_analysis.py) and produces
a 4-panel publication figure:
  (a) expression_subscore vs measured LFC
  (b) full_composite vs measured LFC
  (c) CADD PHRED vs measured LFC
  (d) LFC distribution histogram

Each scatter panel annotates Spearman ρ + 95% CI and draws the Pearson
regression line. Points are colored by |LFC| quintile (viridis, lighter =
lower effect, darker = stronger effect). Axes are consistent: all scatter
panels share the same x-range (LFC domain).

Output: tewhey_correlation_panel.png
Designed to read cleanly at both 1× (9×9 in) and 0.5× (preprint column) scale.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.stats as st

DIR = Path(__file__).parent
PARQUET_IN = DIR / "tewhey_mpra.parquet"
FIG_OUT = DIR / "tewhey_correlation_panel.png"

N_BOOTSTRAP = 1000
rng = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

FIGSIZE = (9.5, 9.0)
BASE_FONT = 11
ANNOT_FONT = 9
TICK_FONT = 9
CMAP = "viridis"
SCATTER_ALPHA = 0.25
SCATTER_SIZE = 6
REGLINE_KW = dict(color="#c0392b", lw=1.5, ls="--", zorder=3, label="Pearson OLS")
HIST_KW = dict(bins=60, color="#2980b9", edgecolor="none", alpha=0.8)

plt.rcParams.update({
    "font.size": BASE_FONT,
    "axes.labelsize": BASE_FONT,
    "axes.titlesize": BASE_FONT,
    "xtick.labelsize": TICK_FONT,
    "ytick.labelsize": TICK_FONT,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------

def _bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOTSTRAP
) -> tuple[float, float]:
    n = len(x)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = st.spearmanr(x[idx], y[idx])
        boot.append(r)
    arr = np.array(boot)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def _correlate(x: pd.Series, y: pd.Series) -> dict:
    xn = pd.to_numeric(x, errors="coerce")
    yn = pd.to_numeric(y, errors="coerce")
    valid = xn.notna() & yn.notna() & np.isfinite(xn.fillna(0)) & np.isfinite(yn.fillna(0))
    x, y = xn, yn
    xv, yv = x[valid].to_numpy(float), y[valid].to_numpy(float)
    n = len(xv)
    if n < 10:
        return {"r_sp": None, "p_sp": None, "ci": (None, None), "n": n,
                "slope": None, "intercept": None, "r_pe": None}
    r_sp, p_sp = st.spearmanr(xv, yv)
    ci = _bootstrap_spearman_ci(xv, yv)
    slope, intercept, r_pe, *_ = st.linregress(xv, yv)
    return {"r_sp": r_sp, "p_sp": p_sp, "ci": ci, "n": n,
            "slope": slope, "intercept": intercept, "r_pe": r_pe,
            "xv": xv, "yv": yv}


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _lfc_color_bins(lfc: pd.Series, n_bins: int = 5) -> np.ndarray:
    """Map |LFC| to [0,1] per-point values for coloring by quintile."""
    abs_lfc = lfc.abs()
    quantiles = abs_lfc.rank(pct=True)
    return quantiles.to_numpy()


def _scatter_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    score_col: str,
    score_label: str,
    color_vals: np.ndarray,
    lfc_xlim: tuple[float, float],
) -> dict:
    valid = df[score_col].notna() & df["mpra_lfc"].notna()
    xv = df.loc[valid, "mpra_lfc"].to_numpy(float)
    yv = df.loc[valid, score_col].to_numpy(float)
    cv = color_vals[valid.to_numpy()]

    stats = _correlate(df.loc[valid, "mpra_lfc"], df.loc[valid, score_col])

    sc = ax.scatter(xv, yv, c=cv, cmap=CMAP, vmin=0, vmax=1,
                    alpha=SCATTER_ALPHA, s=SCATTER_SIZE, rasterized=True, zorder=2)

    # Pearson regression line
    if stats["slope"] is not None:
        x_line = np.linspace(lfc_xlim[0], lfc_xlim[1], 200)
        ax.plot(x_line, stats["slope"] * x_line + stats["intercept"], **REGLINE_KW)

    # Zero-lines
    ax.axhline(0, color="grey", lw=0.5, ls=":", zorder=1)
    ax.axvline(0, color="grey", lw=0.5, ls=":", zorder=1)

    ax.set_xlim(lfc_xlim)
    ax.set_xlabel("MPRA allelic LFC (B − A)", fontsize=BASE_FONT)
    ax.set_ylabel(score_label, fontsize=BASE_FONT)

    # Annotation
    if stats["r_sp"] is not None:
        ci_lo, ci_hi = stats["ci"]
        p_str = f"{stats['p_sp']:.2e}" if stats["p_sp"] is not None else "—"
        ann = (
            f"Spearman ρ = {stats['r_sp']:+.3f}\n"
            f"95% CI [{ci_lo:+.3f}, {ci_hi:+.3f}]\n"
            f"p = {p_str},  n = {stats['n']}"
        )
        ax.text(
            0.04, 0.97, ann,
            transform=ax.transAxes,
            va="top", ha="left",
            fontsize=ANNOT_FONT,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9),
            zorder=5,
        )

    return {"sc": sc, "stats": stats}


def _hist_panel(ax: plt.Axes, lfc: pd.Series) -> None:
    valid = lfc.notna() & np.isfinite(lfc.fillna(0))
    vals = lfc[valid].to_numpy(float)
    ax.hist(vals, **HIST_KW)
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    mean_lfc = vals.mean()
    ax.axvline(mean_lfc, color="#e74c3c", lw=1.0, ls="--",
               label=f"mean = {mean_lfc:.3f}")
    ax.set_xlabel("MPRA allelic LFC (B − A)", fontsize=BASE_FONT)
    ax.set_ylabel("Variant count", fontsize=BASE_FONT)
    ax.set_title("LFC distribution (all variants)", fontsize=BASE_FONT)
    ax.legend(fontsize=ANNOT_FONT)

    sd = vals.std()
    n = len(vals)
    ax.text(
        0.97, 0.97,
        f"n = {n}\nmean = {mean_lfc:.3f}\nSD = {sd:.3f}",
        transform=ax.transAxes,
        va="top", ha="right",
        fontsize=ANNOT_FONT,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#cccccc", alpha=0.9),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

PANELS = [
    ("expression_subscore", "expression_subscore\n(AlphaGenome, signed)", "a"),
    ("full_composite",      "full_composite\n(AlphaGenome, 0–1)",         "b"),
    ("cadd_phred",          "CADD PHRED",                                  "c"),
]


def make_figure(parquet_path: Path = PARQUET_IN, out_path: Path = FIG_OUT) -> None:
    if not parquet_path.exists():
        raise FileNotFoundError(f"{parquet_path} not found — run tewhey_analysis.py first")

    df = pd.read_parquet(parquet_path)
    print(f"Loaded {len(df)} rows from {parquet_path.name}")

    for col in ["mpra_lfc", "expression_subscore", "full_composite", "cadd_phred"]:
        if col not in df.columns:
            print(f"  WARNING: column '{col}' missing — panel will be empty")

    # Color by |LFC| rank (percentile, 0→1) across all variants with valid LFC
    valid_lfc = df["mpra_lfc"].notna() & np.isfinite(df["mpra_lfc"].fillna(0))
    color_vals = pd.Series(np.nan, index=df.index)
    color_vals[valid_lfc] = df.loc[valid_lfc, "mpra_lfc"].abs().rank(pct=True)

    lfc_vals = df.loc[valid_lfc, "mpra_lfc"]
    lfc_pad = 0.2
    lfc_xlim = (lfc_vals.min() - lfc_pad, lfc_vals.max() + lfc_pad)

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE)
    fig.suptitle(
        "Tewhey 2016 MPRA (LCL, GSE75661) — AlphaGenome vs measured allelic LFC",
        fontsize=BASE_FONT + 1, fontweight="bold", y=0.99,
    )

    scatter_sc = None
    for i, (col, ylabel, panel_id) in enumerate(PANELS):
        ax = axes[i // 2][i % 2]
        ax.set_title(f"({panel_id})", loc="left", fontsize=BASE_FONT, fontweight="bold")
        result = _scatter_panel(ax, df, col, ylabel, color_vals.to_numpy(), lfc_xlim)
        if scatter_sc is None:
            scatter_sc = result["sc"]

    # Panel (d): histogram
    ax_hist = axes[1][1]
    ax_hist.set_title("(d)", loc="left", fontsize=BASE_FONT, fontweight="bold")
    _hist_panel(ax_hist, df["mpra_lfc"])

    # Shared colorbar for the three scatter panels
    if scatter_sc is not None:
        cbar = fig.colorbar(
            scatter_sc,
            ax=axes[0, :].tolist() + [axes[1, 0]],
            orientation="vertical",
            fraction=0.015,
            pad=0.02,
            shrink=0.85,
        )
        cbar.set_label("|LFC| percentile rank", fontsize=ANNOT_FONT)
        cbar.ax.tick_params(labelsize=ANNOT_FONT)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    make_figure()
