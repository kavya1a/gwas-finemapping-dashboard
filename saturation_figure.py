"""Phase 1 Task 2 — Saturation empirical CDF figure.

Produces figures/saturation_cdf.png showing that quantile-normalized
AlphaGenome scores saturate on regulatory-enriched variant sets.

Three distributions plotted for each of three modalities
(expression, chromatin, tf_binding):
  (a) Tewhey 2016 MPRA variants (GWAS regulatory loci, n~3259)
  (b) Disease GWAS variants from scored_variants.db (n=767 expression pairs)
  (c) Uniform reference — the theoretical distribution for random
      common variants by construction of the AlphaGenome quantile
      calibration (~300K gnomAD common variants define the background;
      a random draw from that population has uniform quantile scores)

The headline figure (left panel, expression modality) is the CDF of
|signed_score|. A uniform reference would be a straight diagonal line.
Saturation appears as the curve hugging the bottom then jumping to 1
near |score| = 0.9-1.0.

Also saves figures/raw_vs_normalized_dist.png (Phase 2 Task 1 preview):
raw delta distribution vs quantile score distribution on the same Tewhey
600-sample, showing raw is continuous while normalized is bimodal.

Run:
    python saturation_figure.py
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DIR = Path(__file__).parent
FIG_DIR = DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

MODALITIES = ["expression", "chromatin", "tf_binding"]
MOD_LABELS = {
    "expression": "Expression\n(RNA-seq + CAGE + PRO-cap)",
    "chromatin":  "Chromatin accessibility\n(ATAC + DNase + ChIP-Histone)",
    "tf_binding": "TF binding\n(ChIP-TF)",
}

# Colors
C_TEWHEY = "#2166ac"   # blue
C_GWAS   = "#d6604d"   # red-orange
C_UNIF   = "#999999"   # grey (reference)
C_RAW    = "#1a9641"   # green

LW = 1.8
ALPHA_FILL = 0.08


def _ecdf(vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (x, y) for empirical CDF of |vals|."""
    x = np.sort(np.abs(vals))
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def _load_tewhey_scores() -> dict[str, np.ndarray]:
    """Load per-modality signed scores for Tewhey variants from parquet + cache."""
    df = pd.read_parquet(DIR / "tewhey_mpra.parquet")
    result = {}
    # expression_subscore is already in parquet
    expr = df["expression_subscore"].dropna().to_numpy(float)
    result["expression"] = expr

    # chromatin and tf_binding are NOT in the parquet — they're in tewhey_scores_cache.db,
    # which only stores full_composite and expression_subscore.
    # Use the GWAS scored_variants.db modality_scores as a proxy for the per-modality
    # saturation pattern (both datasets were scored against the same calibration).
    # Flag this in the plot as "GWAS variants (proxy for chromatin/TF saturation)".
    return result


def _load_gwas_scores() -> dict[str, np.ndarray]:
    """Load per-modality signed_max_score for disease GWAS variants."""
    conn = sqlite3.connect(DIR / "scored_variants.db")
    result = {}
    for mod in MODALITIES:
        rows = conn.execute(
            "SELECT signed_max_score FROM modality_scores "
            "WHERE modality=? AND signed_max_score IS NOT NULL",
            (mod,),
        ).fetchall()
        result[mod] = np.array([r[0] for r in rows], dtype=float)
    conn.close()
    return result


def _load_raw_deltas() -> np.ndarray:
    """Load raw max_signed_raw from tewhey_raw_delta_cache for expression."""
    conn = sqlite3.connect(DIR / "tewhey_raw_delta_cache.db")
    rows = conn.execute(
        "SELECT max_signed_raw FROM raw_deltas WHERE max_signed_raw IS NOT NULL"
    ).fetchall()
    conn.close()
    return np.array([r[0] for r in rows], dtype=float)


def _saturation_stats(arr: np.ndarray, label: str) -> None:
    abs_arr = np.abs(arr)
    print(f"  {label}: n={len(arr)}, "
          f"|x|>0.9: {(abs_arr>0.9).mean():.1%}, "
          f"|x|<0.5: {(abs_arr<0.5).mean():.1%}, "
          f"median|x|={np.median(abs_arr):.4f}")


# ---------------------------------------------------------------------------
# Figure 1: Saturation CDF — headline figure
# ---------------------------------------------------------------------------

def make_saturation_cdf(tewhey: dict, gwas: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)
    fig.suptitle(
        "Quantile-normalized AlphaGenome scores saturate on regulatory-enriched variant sets",
        fontsize=11, fontweight="bold", y=1.01,
    )

    # Uniform reference CDF (theoretical): |x| ~ Uniform(0,1) → CDF(t) = t
    t_ref = np.linspace(0, 1, 300)

    for ax, mod in zip(axes, MODALITIES):
        gwas_arr = gwas.get(mod, np.array([]))
        tewhey_arr = tewhey.get(mod)   # may be None for chromatin/tf_binding

        # Uniform reference
        ax.plot(t_ref, t_ref, color=C_UNIF, lw=LW, ls="--",
                label="Uniform reference\n(random common variants,\nby AlphaGenome construction)",
                zorder=1)

        # GWAS disease variants
        if len(gwas_arr) > 0:
            xg, yg = _ecdf(gwas_arr)
            ax.plot(xg, yg, color=C_GWAS, lw=LW,
                    label=f"Disease GWAS variants\n(n={len(gwas_arr)} variant-disease pairs)",
                    zorder=3)
            ax.fill_between(xg, 0, yg, color=C_GWAS, alpha=ALPHA_FILL)

        # Tewhey variants (only available for expression in parquet)
        if tewhey_arr is not None and len(tewhey_arr) > 0:
            xt, yt = _ecdf(tewhey_arr)
            ax.plot(xt, yt, color=C_TEWHEY, lw=LW,
                    label=f"Tewhey MPRA loci\n(n={len(tewhey_arr)} GWAS regulatory variants)",
                    zorder=4)
            ax.fill_between(xt, 0, yt, color=C_TEWHEY, alpha=ALPHA_FILL)
        elif mod != "expression":
            # Annotate that Tewhey chromatin/TF not separately cached
            ax.text(0.5, 0.35, "Tewhey per-modality scores\nnot separately cached\n"
                    "(see GWAS curve for proxy)",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=7.5, color=C_TEWHEY,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor=C_TEWHEY, alpha=0.8))

        # Saturation reference line at |score| = 0.9
        ax.axvline(0.9, color="#888888", lw=0.8, ls=":", zorder=0)
        ax.text(0.905, 0.08, "|score|=0.9", fontsize=7, color="#666666", rotation=90)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("|Signed quantile score|", fontsize=9)
        ax.set_title(MOD_LABELS[mod], fontsize=9)
        ax.grid(axis="both", lw=0.3, alpha=0.5)

        if ax is axes[0]:
            ax.set_ylabel("Cumulative fraction of variants", fontsize=9)
            ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)

    fig.text(
        0.5, -0.04,
        "Saturation: both regulatory-enriched sets shift entirely to |score| > 0.9,\n"
        "collapsing the predictor to a near-binary signal (by construction of the genome-wide quantile calibration).",
        ha="center", fontsize=8.5, style="italic", color="#333333",
    )

    out = FIG_DIR / "saturation_cdf.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Figure 2: Raw delta vs normalized score distribution (Phase 2 Task 1 preview)
# ---------------------------------------------------------------------------

def make_raw_vs_normalized(tewhey_expr: np.ndarray, raw_deltas: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(
        "Distribution of quantile-normalized scores vs raw model outputs (Tewhey variants)",
        fontsize=10, fontweight="bold",
    )

    # Panel A: quantile-normalized expression_subscore
    ax = axes[0]
    ax.hist(tewhey_expr, bins=80, color=C_TEWHEY, alpha=0.75, edgecolor="none")
    ax.set_xlabel("Quantile-normalized expression_subscore", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.set_title(f"Quantile score (bimodal)\nn = {len(tewhey_expr):,}", fontsize=9)
    pct_extreme = (np.abs(tewhey_expr) > 0.9).mean()
    ax.text(0.05, 0.95, f"{pct_extreme:.1%} of variants\nhave |score| > 0.9",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_TEWHEY, alpha=0.9))
    ax.axvline(0.9, color="grey", lw=0.8, ls=":")
    ax.axvline(-0.9, color="grey", lw=0.8, ls=":")

    # Panel B: raw delta
    ax = axes[1]
    # Clip extreme outliers for display (keep 99th percentile range)
    p99 = np.percentile(np.abs(raw_deltas), 99)
    clipped = raw_deltas[np.abs(raw_deltas) <= p99 * 1.1]
    ax.hist(clipped, bins=80, color=C_RAW, alpha=0.75, edgecolor="none")
    ax.set_xlabel("Raw max signed expression delta", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.set_title(f"Raw delta (continuous)\nn = {len(raw_deltas):,}", fontsize=9)
    pct_small = (np.abs(raw_deltas) < 0.1).mean()
    ax.text(0.55, 0.95, f"{pct_small:.1%} of variants\nhave |delta| < 0.1\n(continuous gradient)",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_RAW, alpha=0.9))

    out = FIG_DIR / "raw_vs_normalized_dist.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading scores...")
    tewhey = _load_tewhey_scores()
    gwas   = _load_gwas_scores()
    raw    = _load_raw_deltas()

    print("\n=== Saturation statistics ===")
    for mod in MODALITIES:
        print(f"\n  [{mod}]")
        if mod in tewhey:
            _saturation_stats(tewhey[mod], "Tewhey")
        _saturation_stats(gwas[mod], f"GWAS ({mod})")

    print(f"\n  [raw delta (expression)]")
    _saturation_stats(raw, "Raw delta")

    print("\nGenerating figures...")
    make_saturation_cdf(tewhey, gwas)
    make_raw_vs_normalized(tewhey["expression"], raw)
    print("\nDone. Figures in figures/")


if __name__ == "__main__":
    main()
