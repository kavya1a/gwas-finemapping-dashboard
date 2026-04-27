"""Phase 2 figures — Demonstrating the fix.

Produces:
  figures/phase2_distribution.png   — raw delta vs normalized score distributions
  figures/phase2_lfc_bins.png       — Spearman correlation by |LFC| bin, raw vs normalized
  figures/phase2_combined.png       — 3-panel combined figure for README

Run:
    python phase2_figures.py
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats as st

DIR = Path(__file__).parent
FIG_DIR = DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

SEED = 42
N_BOOT = 1000

C_NORM = "#2166ac"   # blue  — normalized score
C_RAW  = "#1a9641"   # green — raw delta
C_ZERO = "#888888"   # grey  — zero reference


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (full_tewhey_df, raw_delta_sample_df)."""
    df = pd.read_parquet(DIR / "tewhey_mpra.parquet")

    conn = sqlite3.connect(DIR / "tewhey_raw_delta_cache.db")
    raw = pd.read_sql(
        "SELECT rsid, max_signed_raw, mean_signed_raw FROM raw_deltas "
        "WHERE max_signed_raw IS NOT NULL",
        conn,
    )
    conn.close()

    sample = df.merge(raw, on="rsid", how="inner")
    return df, sample


# ---------------------------------------------------------------------------
# Bootstrap Spearman CI
# ---------------------------------------------------------------------------

def _bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                         n_boot: int = N_BOOT) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    n = len(x)
    boot = [st.spearmanr(x[idx := rng.integers(0, n, n)], y[idx])[0]
            for _ in range(n_boot)]
    arr = np.array(boot)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def spearman_by_lfc_bin(
    lfc: np.ndarray, score: np.ndarray, bins: list[float]
) -> list[dict]:
    results = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (np.abs(lfc) >= lo) & (np.abs(lfc) < hi)
        xv, yv = lfc[mask], score[mask]
        n = mask.sum()
        if n < 10:
            results.append({"lo": lo, "hi": hi, "n": n,
                             "r": None, "ci_lo": None, "ci_hi": None})
            continue
        r, _ = st.spearmanr(xv, yv)
        ci_lo, ci_hi = _bootstrap_spearman(xv, yv)
        results.append({"lo": lo, "hi": hi, "n": n, "r": r,
                         "ci_lo": ci_lo, "ci_hi": ci_hi})
    return results


# ---------------------------------------------------------------------------
# Figure A: distributions
# ---------------------------------------------------------------------------

def make_distribution_fig(full_df: pd.DataFrame, sample: pd.DataFrame) -> None:
    norm_scores = full_df["expression_subscore"].dropna().to_numpy(float)
    raw_deltas  = sample["max_signed_raw"].to_numpy(float)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(
        "Quantile-normalized scores collapse to bimodal; raw deltas retain a continuous gradient",
        fontsize=10, fontweight="bold",
    )

    # Panel A — normalized
    ax = axes[0]
    ax.hist(norm_scores, bins=80, color=C_NORM, alpha=0.80, edgecolor="none")
    pct = (np.abs(norm_scores) > 0.9).mean()
    ax.axvline(0.9,  color=C_ZERO, lw=0.9, ls=":")
    ax.axvline(-0.9, color=C_ZERO, lw=0.9, ls=":")
    ax.set_xlabel("Quantile-normalized expression score", fontsize=9)
    ax.set_ylabel("Variant count", fontsize=9)
    ax.set_title(f"Normalized score — Tewhey variants  (n = {len(norm_scores):,})", fontsize=9)
    ax.text(0.05, 0.95,
            f"{pct:.1%} of variants\nhave |score| > 0.9\n→ effectively binary",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_NORM, alpha=0.9))

    # Panel B — raw delta
    ax = axes[1]
    p99 = np.percentile(np.abs(raw_deltas), 99.5)
    clipped = raw_deltas[np.abs(raw_deltas) <= p99 * 1.05]
    ax.hist(clipped, bins=80, color=C_RAW, alpha=0.80, edgecolor="none")
    pct_small = (np.abs(raw_deltas) < 0.1).mean()
    ax.axvline(0, color=C_ZERO, lw=0.8, ls="--")
    ax.set_xlabel("Raw max signed expression delta", fontsize=9)
    ax.set_ylabel("Variant count", fontsize=9)
    ax.set_title(f"Raw delta — stratified sample  (n = {len(raw_deltas):,})", fontsize=9)
    ax.text(0.55, 0.95,
            f"{pct_small:.1%} of variants\nhave |Δ| < 0.1\n→ continuous gradient",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_RAW, alpha=0.9))

    fig.tight_layout()
    out = FIG_DIR / "phase2_distribution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Figure B: correlation by |LFC| bin
# ---------------------------------------------------------------------------

LFC_BINS = [0.0, 0.05, 0.10, 0.20, 0.50, 3.1]
BIN_LABELS = ["0–0.05", "0.05–0.10", "0.10–0.20", "0.20–0.50", "0.50+"]


def make_lfc_bin_fig(full_df: pd.DataFrame, sample: pd.DataFrame) -> None:
    # Normalized: full 3259-variant set
    lfc_norm  = full_df["mpra_lfc"].dropna().to_numpy(float)
    expr_norm = full_df.loc[full_df["mpra_lfc"].notna(), "expression_subscore"]
    valid_n = full_df["mpra_lfc"].notna() & full_df["expression_subscore"].notna()
    lfc_n  = full_df.loc[valid_n, "mpra_lfc"].to_numpy(float)
    scr_n  = full_df.loc[valid_n, "expression_subscore"].to_numpy(float)

    # Raw: 600-variant stratified sample
    lfc_r  = sample["mpra_lfc"].to_numpy(float)
    scr_r  = sample["max_signed_raw"].to_numpy(float)

    bins_norm = spearman_by_lfc_bin(lfc_n, scr_n, LFC_BINS)
    bins_raw  = spearman_by_lfc_bin(lfc_r, scr_r, LFC_BINS)

    # Print table
    print("\n=== Correlation by |LFC| bin ===")
    print(f"{'|LFC| bin':<12} {'norm r':>8} {'norm n':>8} {'raw r':>8} {'raw n':>8}")
    for bn, br in zip(bins_norm, bins_raw):
        rn = f"{bn['r']:+.3f}" if bn['r'] is not None else "  —  "
        rr = f"{br['r']:+.3f}" if br['r'] is not None else "  —  "
        print(f"{BIN_LABELS[bins_norm.index(bn)]:<12} {rn:>8} {bn['n']:>8} {rr:>8} {br['n']:>8}")

    x = np.arange(len(BIN_LABELS))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_title(
        "Spearman ρ (score vs MPRA LFC) by |LFC| bin\n"
        "Raw deltas correlate across the full range; normalized scores fail below |LFC| = 0.2",
        fontsize=10, fontweight="bold",
    )

    for j, (bins, color, label, hatch) in enumerate([
        (bins_norm, C_NORM, f"Normalized score (n={len(lfc_n):,})", ""),
        (bins_raw,  C_RAW,  f"Raw delta (n={len(lfc_r):,})",        "//"),
    ]):
        rs    = [b["r"]    if b["r"]    is not None else 0.0 for b in bins]
        ci_lo = [b["ci_lo"] if b["ci_lo"] is not None else 0.0 for b in bins]
        ci_hi = [b["ci_hi"] if b["ci_hi"] is not None else 0.0 for b in bins]
        ns    = [b["n"] for b in bins]
        err_lo = [r - lo for r, lo in zip(rs, ci_lo)]
        err_hi = [hi - r for r, hi in zip(rs, ci_hi)]

        bars = ax.bar(
            x + j * width - width / 2, rs, width,
            color=color, alpha=0.82, hatch=hatch,
            label=label, edgecolor="white", lw=0.5,
        )
        ax.errorbar(
            x + j * width - width / 2, rs,
            yerr=[err_lo, err_hi],
            fmt="none", color="black", capsize=3, lw=1.0,
        )
        # Annotate n
        for xi, (r, n) in enumerate(zip(rs, ns)):
            if n >= 10:
                ax.text(xi + j * width - width / 2, max(r, 0) + 0.008,
                        f"n={n}", ha="center", va="bottom", fontsize=6.5, color="#333333")

    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"|LFC|\n{b}" for b in BIN_LABELS], fontsize=9)
    ax.set_ylabel("Spearman ρ (vs MPRA allelic LFC)", fontsize=9)
    ax.set_xlabel("|MPRA allelic LFC| bin", fontsize=9)
    ax.set_ylim(-0.15, 0.55)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    # Annotation box
    ax.text(0.98, 0.97,
            "Error bars: 95% bootstrap CI (n=1000)\n"
            "Normalized scores use full Tewhey set\n"
            "Raw deltas use stratified sample (120/quintile)",
            transform=ax.transAxes, va="top", ha="right", fontsize=7.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#aaaaaa", alpha=0.9))

    fig.tight_layout()
    out = FIG_DIR / "phase2_lfc_bins.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Figure C: combined 3-panel for README
# ---------------------------------------------------------------------------

def make_combined_fig(full_df: pd.DataFrame, sample: pd.DataFrame) -> None:
    norm_scores = full_df["expression_subscore"].dropna().to_numpy(float)
    raw_deltas  = sample["max_signed_raw"].to_numpy(float)

    valid_n = full_df["mpra_lfc"].notna() & full_df["expression_subscore"].notna()
    lfc_n  = full_df.loc[valid_n, "mpra_lfc"].to_numpy(float)
    scr_n  = full_df.loc[valid_n, "expression_subscore"].to_numpy(float)
    lfc_r  = sample["mpra_lfc"].to_numpy(float)
    scr_r  = sample["max_signed_raw"].to_numpy(float)

    bins_norm = spearman_by_lfc_bin(lfc_n, scr_n, LFC_BINS)
    bins_raw  = spearman_by_lfc_bin(lfc_r, scr_r, LFC_BINS)

    fig = plt.figure(figsize=(14, 4.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1.4], wspace=0.35)

    # Panel 1 — normalized histogram
    ax1 = fig.add_subplot(gs[0])
    ax1.hist(norm_scores, bins=80, color=C_NORM, alpha=0.8, edgecolor="none")
    ax1.axvline(0.9,  color=C_ZERO, lw=0.8, ls=":")
    ax1.axvline(-0.9, color=C_ZERO, lw=0.8, ls=":")
    pct = (np.abs(norm_scores) > 0.9).mean()
    ax1.text(0.05, 0.95, f"{pct:.1%} have\n|score| > 0.9",
             transform=ax1.transAxes, va="top", fontsize=8,
             bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_NORM, alpha=0.9))
    ax1.set_xlabel("Quantile-normalized expression score", fontsize=8.5)
    ax1.set_ylabel("Count", fontsize=8.5)
    ax1.set_title(f"(A) Normalized score\nn = {len(norm_scores):,} Tewhey variants", fontsize=9)

    # Panel 2 — raw delta histogram
    ax2 = fig.add_subplot(gs[1])
    p99 = np.percentile(np.abs(raw_deltas), 99.5)
    clipped = raw_deltas[np.abs(raw_deltas) <= p99 * 1.05]
    ax2.hist(clipped, bins=80, color=C_RAW, alpha=0.8, edgecolor="none")
    ax2.axvline(0, color=C_ZERO, lw=0.8, ls="--")
    pct_small = (np.abs(raw_deltas) < 0.1).mean()
    ax2.text(0.52, 0.95, f"{pct_small:.1%} have\n|Δ| < 0.1",
             transform=ax2.transAxes, va="top", fontsize=8,
             bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_RAW, alpha=0.9))
    ax2.set_xlabel("Raw max signed expression delta", fontsize=8.5)
    ax2.set_ylabel("Count", fontsize=8.5)
    ax2.set_title(f"(B) Raw delta\nn = {len(raw_deltas):,} stratified sample", fontsize=9)

    # Panel 3 — correlation by bin
    ax3 = fig.add_subplot(gs[2])
    x = np.arange(len(BIN_LABELS))
    width = 0.35
    for j, (bins, color, label, hatch) in enumerate([
        (bins_norm, C_NORM, "Normalized score", ""),
        (bins_raw,  C_RAW,  "Raw delta", "//"),
    ]):
        rs    = [b["r"]    if b["r"]    is not None else 0.0 for b in bins]
        ci_lo = [b["ci_lo"] if b["ci_lo"] is not None else 0.0 for b in bins]
        ci_hi = [b["ci_hi"] if b["ci_hi"] is not None else 0.0 for b in bins]
        err_lo = [r - lo for r, lo in zip(rs, ci_lo)]
        err_hi = [hi - r for r, hi in zip(rs, ci_hi)]
        ax3.bar(x + j*width - width/2, rs, width,
                color=color, alpha=0.82, hatch=hatch,
                label=label, edgecolor="white", lw=0.5)
        ax3.errorbar(x + j*width - width/2, rs,
                     yerr=[err_lo, err_hi],
                     fmt="none", color="black", capsize=3, lw=1.0)

    ax3.axhline(0, color="black", lw=0.8)
    ax3.set_xticks(x)
    ax3.set_xticklabels(BIN_LABELS, fontsize=8)
    ax3.set_xlabel("|MPRA LFC| bin", fontsize=8.5)
    ax3.set_ylabel("Spearman ρ (vs MPRA LFC)", fontsize=8.5)
    ax3.set_title("(C) Spearman ρ by |LFC| bin\n(error bars: 95% bootstrap CI)", fontsize=9)
    ax3.set_ylim(-0.15, 0.55)
    ax3.legend(fontsize=8, loc="upper left")
    ax3.grid(axis="y", lw=0.3, alpha=0.5)

    fig.suptitle(
        "Raw expression deltas restore directional signal lost by quantile normalization",
        fontsize=11, fontweight="bold", y=1.02,
    )

    out = FIG_DIR / "phase2_combined.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading data...")
    full_df, sample = load_data()
    print(f"  Full Tewhey: {len(full_df)} rows, {full_df['expression_subscore'].notna().sum()} scored")
    print(f"  Raw delta sample: {len(sample)} matched rows")

    make_distribution_fig(full_df, sample)
    make_lfc_bin_fig(full_df, sample)
    make_combined_fig(full_df, sample)
    print("\nPhase 2 figures complete.")


if __name__ == "__main__":
    main()
