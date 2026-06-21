"""Matched re-calibration under DNase / ATAC (accessibility) output types.

STATUS: NOT YET RUN — author-suggested next experiment, scaffolded only.
No accessibility cache is committed to this repo, and this script has not been
validated end-to-end against the live SDK. It mirrors the committed expression
pipeline (`build_matched_calibration.py` + `analyze_matched_calibration.py`) but
scores under DNase/ATAC output types instead of RNA_SEQ + CAGE + PROCAP.

Why it matters
--------------
The headline result is that AlphaGenome's published per-track expression quantile,
composed with a max-over-tracks aggregation, saturates — and that a matched null
(same max-over-tracks statistic, re-derived on common variants) removes the
artifact. Accessibility (DNase/ATAC) is the natural second test:

  * AlphaGenome exposes far FEWER accessibility tracks than expression tracks.
    The order-statistic inflation is 1 - (1 - p)^n in the number of tracks n, so
    with fewer accessibility tracks the *unmatched* max should already saturate
    LESS than expression does (expression hit 41.4% on the common-variant null;
    accessibility should sit lower). This script measures that directly.
  * If the matched-vs-unmatched gap shrinks the way the track-count argument
    predicts, that is independent confirmation that the mechanism is order
    statistics, not anything specific to the expression tracks.

What it does (identical structure to the expression analysis)
-------------------------------------------------------------
  1. Sample the SAME common-variant windows (seed = 2026, reusing
     `build_matched_calibration.collect_variants`) so the null is comparable.
  2. Score each variant through AlphaGenome, but keep the max-over-tracks signed
     `raw_score` across ACCESSIBILITY output types {ATAC, DNASE} and the peak
     accessibility track's published `quantile_score`.
  3. Build the matched accessibility null and report saturation |·|>0.9 / >0.5 for
     (a) matched max-over-tracks raw delta and (b) the published single-track
     accessibility quantile — the same smoking-gun table as expression.
  4. Re-quantile the Tewhey panel's accessibility deltas against that null and
     report Spearman vs MPRA LFC, exactly as `analyze_matched_calibration.py`
     does for expression. Compare matched vs unmatched side by side.

Run (requires an AlphaGenome API key and hours of API time):
    export ALPHAGENOME_API_KEY=...
    python dnase_atac_recalibration.py            # build accessibility null + Tewhey re-cal

It FAILS LOUDLY if the API key (or, for --post, a previously-built accessibility
cache) is absent. It never emits placeholder or simulated numbers — every value
comes from a real scoring run.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from pathlib import Path

DIR = Path(__file__).parent
ACCESS_CACHE = DIR / "dnase_atac_calibration_cache.db"          # created on first real run
ACCESS_NULL_PARQUET = DIR / "dnase_atac_calibration_null.parquet"
TEWHEY_ACCESS_CACHE = DIR / "tewhey_dnase_atac_delta_cache.db"  # created on first real run
OUT_CSV = DIR / "dnase_atac_comparison.csv"

# Mirror of build_matched_calibration.EXPRESSION_OUTPUT_TYPES, for accessibility.
ACCESSIBILITY_OUTPUT_TYPES = {"ATAC", "DNASE"}


def _fail(msg: str, code: int = 2) -> None:
    """Loud, explanatory exit. Never returns; never produces numbers."""
    print("=" * 78, file=sys.stderr)
    print("dnase_atac_recalibration.py — cannot run", file=sys.stderr)
    print("=" * 78, file=sys.stderr)
    print(msg.rstrip(), file=sys.stderr)
    sys.exit(code)


def extract_accessibility_max_signed(tidy_df, profile) -> dict:
    """Per-variant accessibility stats from one tidy_df.

    Mirrors `build_matched_calibration.extract_max_signed_with_quantile`, but
    filters to ACCESSIBILITY_OUTPUT_TYPES instead of expression. Returns
    {raw_max_signed_delta, single_track_quantile, n_access_tracks}.
    """
    from scoring.tissue_config import filter_tracks  # lazy: only needed on a real run

    empty = {"raw_max_signed_delta": None, "single_track_quantile": None,
             "n_access_tracks": 0}
    if tidy_df is None or getattr(tidy_df, "empty", True):
        return empty

    filtered = filter_tracks(tidy_df, profile)
    if filtered is None or filtered.empty:
        return empty
    ot = filtered["output_type"].astype(str).str.replace("OutputType.", "", regex=False)
    acc = filtered[ot.isin(ACCESSIBILITY_OUTPUT_TYPES)]
    raw = acc["raw_score"].dropna()
    if raw.empty:
        return empty

    idx_max = raw.abs().idxmax()
    max_signed = float(raw.loc[idx_max])
    single_track_q = (
        float(acc.loc[idx_max, "quantile_score"])
        if "quantile_score" in acc.columns and acc.loc[idx_max, "quantile_score"] is not None
        else None
    )
    return {
        "raw_max_signed_delta": max_signed,
        "single_track_quantile": single_track_q,
        "n_access_tracks": int(len(raw)),
    }


def run(seed: int, cap: int, post: bool) -> None:
    # ── prerequisite checks: fail loud before any heavy import or compute ──────
    if post:
        if not ACCESS_CACHE.exists():
            _fail(
                f"--post rebuilds the null + figure from a committed accessibility cache,\n"
                f"but {ACCESS_CACHE.name} does not exist. This experiment has not been run,\n"
                f"so there is no cache to summarize. Run without --post (and with an API key)\n"
                f"to score the accessibility tracks first."
            )
    else:
        if not os.environ.get("ALPHAGENOME_API_KEY", ""):
            _fail(
                "ALPHAGENOME_API_KEY is not set.\n\n"
                "This script scores DNase/ATAC accessibility tracks for ~6,000 common\n"
                "variants plus the Tewhey panel through the AlphaGenome API (hours of\n"
                "runtime). It has no committed cache to fall back on, because the\n"
                "experiment has not been run yet. Set ALPHAGENOME_API_KEY and re-run,\n"
                "or see the README 'Next experiments' section for what this produces."
            )

    # ── real run path (reuses the committed expression machinery) ─────────────
    # Imported lazily so the fail-loud paths above work without the SDK installed.
    from build_matched_calibration import (  # noqa: F401  (reused, output-type-agnostic)
        collect_variants,
        load_k562_profile,
        score_one,
    )
    from scoring.score_worker import ScoreWorker

    profile = load_k562_profile()
    variants = collect_variants(seed=seed)[:cap]
    api_key = os.environ["ALPHAGENOME_API_KEY"]
    worker = ScoreWorker(api_key)

    _init_cache()
    print(f"Scoring {len(variants)} common variants under {sorted(ACCESSIBILITY_OUTPUT_TYPES)} ...")
    try:
        for i, v in enumerate(variants):
            # score_one returns the tidy_df; we re-extract the ACCESSIBILITY summary
            # rather than the expression one. (See module docstring step 2.)
            tidy_df = score_one(worker, v, profile)
            entry = extract_accessibility_max_signed(tidy_df, profile)
            _save(v["rsid"], v, entry)
            if (i + 1) % 100 == 0:
                print(f"  [{i + 1}/{len(variants)}]")
    finally:
        worker.close()

    # NOTE: null construction, the matched-vs-unmatched saturation table, and the
    # Tewhey re-quantiling step mirror build_null_and_figure() +
    # analyze_matched_calibration.py exactly, against ACCESS_NULL_PARQUET /
    # TEWHEY_ACCESS_CACHE. They are intentionally left for the first real run so
    # that no number is written without a validated end-to-end scoring pass.
    raise SystemExit(
        "Accessibility scoring cache written. Null construction + Tewhey re-quantiling "
        "follow the expression pipeline 1:1 and should be wired against the freshly "
        "scored cache before reporting any numbers — do not infer results from a partial run."
    )


def _init_cache() -> None:
    conn = sqlite3.connect(ACCESS_CACHE)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS scores (
               rsid TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
               maf REAL, raw_max_signed_delta REAL, single_track_quantile REAL,
               n_access_tracks INTEGER, error TEXT, scored_at INTEGER)"""
    )
    conn.commit()
    conn.close()


def _save(rsid: str, v: dict, entry: dict) -> None:
    conn = sqlite3.connect(ACCESS_CACHE)
    conn.execute(
        "INSERT OR REPLACE INTO scores (rsid,chrom,pos,ref,alt,maf,"
        "raw_max_signed_delta,single_track_quantile,n_access_tracks,error,scored_at) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (rsid, v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"), v.get("maf"),
         entry.get("raw_max_signed_delta"), entry.get("single_track_quantile"),
         entry.get("n_access_tracks"), entry.get("error"), int(time.time())),
    )
    conn.commit()
    conn.close()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=2026, help="Window-sampling seed (default 2026)")
    p.add_argument("--cap", type=int, default=6000, help="Max common variants to score")
    p.add_argument("--post", action="store_true",
                   help="Rebuild null + figure from an existing accessibility cache (none committed yet)")
    args = p.parse_args()
    run(seed=args.seed, cap=args.cap, post=args.post)


if __name__ == "__main__":
    main()
