"""Which transform do AlphaGenome's shipped phred-scale quantiles correspond to?

STATUS: NOT YET RUN — author-suggested next experiment, scaffolded only.

The Nature-paper SDK used here (AlphaGenome v0.6.1) exposes no phred-scaled
output — verified across the SDK source, all ten tagged releases, the issue
tracker, the docs, and the supplement. If a newer SDK ships phred-scale quantile
values, this script settles a concrete question: are those shipped values the
**single-track** phred transform, or the **matched-recalibration** phred?

  (a) single-track phred  = -10·log10(1 - q_single + eps)
        where q_single is the published per-track quantile (expression_subscore)
        mapped to [0, 1). This is what you'd get by phred-transforming the
        existing saturated single-track quantile — it inherits the saturation.

  (b) matched-recal phred = -10·log10(1 - q_matched + eps)
        where q_matched is this work's matched-statistic quantile (max-over-tracks
        null). This is `phred_empirical` in docs/matched_calibration.md.

Both are computed here from committed data. The unknown is the SDK's own phred
column. The script ranks the SDK values against (a) and (b) by Spearman and by
absolute agreement, and reports which the shipped values track — i.e., whether
the SDK's phred inherits the single-track saturation or the matched recovery.

Run:
    # Needs a source of SDK phred values for the Tewhey rsIDs, one of:
    #   --sdk-phred-csv PATH   CSV with columns: rsid, sdk_phred
    # (A newer SDK that returns phred directly can be wired into _load_sdk_phred.)
    python phred_scale_check.py --sdk-phred-csv sdk_phred_tewhey.csv

It FAILS LOUDLY if no SDK phred source is available. (a) and (b) come from
committed caches; the SDK column must come from a real SDK/run — this script
never invents it.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DIR = Path(__file__).parent
TEWHEY_PARQUET = DIR / "tewhey_mpra.parquet"
TEWHEY_RAW_CACHE = DIR / "tewhey_raw_delta_cache.db"
NULL_PARQUET = DIR / "matched_calibration_null.parquet"

EPS = 1e-6


def _fail(msg: str, code: int = 2) -> None:
    print("=" * 78, file=sys.stderr)
    print("phred_scale_check.py — cannot run", file=sys.stderr)
    print("=" * 78, file=sys.stderr)
    print(msg.rstrip(), file=sys.stderr)
    sys.exit(code)


def _phred(signed_quantile: np.ndarray) -> np.ndarray:
    """-10·log10(1 - q + eps); matches docs/matched_calibration.md phred_empirical."""
    return -10.0 * np.log10(1.0 - signed_quantile + EPS)


def _matched_quantile(x: np.ndarray, null_sorted: np.ndarray) -> np.ndarray:
    n_lt = np.searchsorted(null_sorted, x, side="left")
    n_eq = np.searchsorted(null_sorted, x, side="right") - n_lt
    return 2.0 * ((n_lt + 0.5 * n_eq) / len(null_sorted)) - 1.0


def _load_sdk_phred(args) -> pd.DataFrame:
    """Return a DataFrame [rsid, sdk_phred] from a real SDK source, or fail loud.

    There is deliberately no fallback: if no phred source is supplied this raises,
    because AlphaGenome v0.6.1 ships no phred output and fabricating one would
    defeat the entire purpose of the check.
    """
    if args.sdk_phred_csv:
        path = Path(args.sdk_phred_csv)
        if not path.exists():
            _fail(f"--sdk-phred-csv {path} does not exist.")
        df = pd.read_csv(path)
        missing = {"rsid", "sdk_phred"} - set(df.columns)
        if missing:
            _fail(f"{path} is missing required column(s): {sorted(missing)}.")
        return df[["rsid", "sdk_phred"]].dropna()

    _fail(
        "No source of SDK phred-scale values was provided.\n\n"
        "AlphaGenome v0.6.1 (pinned in requirements.txt) exposes no phred-scaled\n"
        "output, so there is nothing to compare against out of the box. To run this\n"
        "check, supply one of:\n"
        "  --sdk-phred-csv PATH   a CSV (columns: rsid, sdk_phred) exported from an\n"
        "                         SDK version that ships phred quantiles.\n\n"
        "The single-track and matched-recalibration phred transforms (a) and (b) are\n"
        "computed from committed data, but the SDK's own phred column must come from a\n"
        "real SDK — this script will not invent it."
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sdk-phred-csv", default=None,
                   help="CSV with columns rsid, sdk_phred (from an SDK that ships phred)")
    args = p.parse_args()

    sdk = _load_sdk_phred(args)  # fails loud if absent

    # (a) and (b) from committed caches.
    parq = pd.read_parquet(TEWHEY_PARQUET)[["rsid", "expression_subscore"]]
    con = sqlite3.connect(TEWHEY_RAW_CACHE)
    raw = pd.read_sql("SELECT rsid, max_signed_raw FROM raw_deltas", con)
    con.close()
    null_sorted = np.sort(
        pd.read_parquet(NULL_PARQUET)["raw_max_signed_delta"].dropna().to_numpy(float)
    )

    df = parq.merge(raw, on="rsid", how="inner").merge(sdk, on="rsid", how="inner")
    df = df.dropna(subset=["expression_subscore", "max_signed_raw", "sdk_phred"]).copy()
    if df.empty:
        _fail("No Tewhey rsIDs overlap between the SDK phred CSV and the committed caches.")

    df["phred_single_track"] = _phred(df["expression_subscore"].to_numpy(float))
    df["phred_matched"] = _phred(
        _matched_quantile(df["max_signed_raw"].to_numpy(float), null_sorted)
    )

    r_single, _ = spearmanr(df["sdk_phred"], df["phred_single_track"])
    r_matched, _ = spearmanr(df["sdk_phred"], df["phred_matched"])

    print(f"n overlapping variants: {len(df)}")
    print(f"Spearman(sdk_phred, single-track phred) = {r_single:+.4f}")
    print(f"Spearman(sdk_phred, matched-recal phred) = {r_matched:+.4f}")
    verdict = "single-track" if r_single > r_matched else "matched-recalibration"
    print(f"\nShipped SDK phred values track the {verdict} transform more closely.")


if __name__ == "__main__":
    main()
