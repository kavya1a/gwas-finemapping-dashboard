"""Build a matched-statistic calibration distribution for the K562 raw delta.

AlphaGenome's published quantile_score is calibrated against a single-track
common-variant null. Our pipeline applies a max-over-tracks summary on top of
that, which inflates saturation. This script builds the matched null: per-variant
K562/blood-lineage RNA_SEQ + CAGE + PROCAP expression raw deltas, aggregated
three ways (max / mean / median), on ~5,000 random common autosomal variants.
Component 3 then quantile-ranks each Tewhey variant against each of the three
recipes.

Resume-aware: rows with a successful raw_max but no raw_mean (e.g. cached by
an earlier version of this script) will be re-scored to backfill mean+median.
Prior errors are not retried.

Sampling: 66 windows of 50 kb anchored at uniform-random offsets (seed 2026)
within autosome interiors (>=1 Mb from telomeres), allocated proportional
to chromosome length. gnomAD v3 AF>0.01 biallelic SNVs only. Tewhey rsids
excluded.

Scoring: same K562 tissue profile, same expression-modality filter,
same 60 s per-variant timeout as extract_raw_deltas.py. 4-way thread-pool
parallelism on top, with a write lock around SQLite. Resume-able via
matched_calibration_cache.db.

Run:
    python build_matched_calibration.py --pilot   # 5 variants serial (10-15 min)
    python build_matched_calibration.py            # full ~5,000-variant parallel run
    python build_matched_calibration.py --post     # post-process cache only
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import re
import sqlite3
import sys
import threading
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import yaml
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

DIR          = Path(__file__).parent
CACHE_DB     = DIR / "matched_calibration_cache.db"
NULL_PARQUET = DIR / "matched_calibration_null.parquet"
FIG_DIR      = DIR / "figures"
HIST_FIG     = FIG_DIR / "matched_calibration_histogram.png"
TEWHEY_PARQUET = DIR / "tewhey_mpra.parquet"
CONFIG_PATH    = DIR / "config.yaml"

# Sampling config (reproducible from seed alone)
SEED          = 2026
N_WINDOWS     = 66
WINDOW_BP     = 50_000
MIN_TELO_BP   = 1_000_000
MAF_THRESHOLD = 0.01

# Scoring config (mirrors batch_score.py / extract_raw_deltas.py)
VARIANT_TIMEOUT_SECS = 60
PARALLEL_WORKERS     = 4
PILOT_N              = 5

# Diagnostic stop: if common variants already saturate at >30% after 200, halt.
SAT_CHECK_AFTER       = 200
SAT_HALT_FRACTION     = 0.30

EXPRESSION_OUTPUT_TYPES = {"RNA_SEQ", "CAGE", "PROCAP"}
VALID_BASE = re.compile(r"^[ACGTacgt]$")

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
_GNOMAD_QUERY = """
query RegionVariants($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r3) {
      variant_id
      rsid
      genome { af }
    }
  }
}
"""

# GRCh38 autosome lengths (chr1..chr22; sex chromosomes excluded per spec).
GRCH38_AUTOSOMES: dict[str, int] = {
    "1":  248_956_422, "2":  242_193_529, "3":  198_295_559, "4":  190_214_555,
    "5":  181_538_259, "6":  170_805_979, "7":  159_345_973, "8":  145_138_636,
    "9":  138_394_717, "10": 133_797_422, "11": 135_086_622, "12": 133_275_309,
    "13": 114_364_328, "14": 107_043_718, "15": 101_991_189, "16":  90_338_345,
    "17":  83_257_441, "18":  80_373_285, "19":  58_617_616, "20":  64_444_167,
    "21":  46_709_983, "22":  50_818_468,
}
TOTAL_AUTOSOME_BP = sum(GRCH38_AUTOSOMES.values())


# ---------------------------------------------------------------------------
# Window sampling — reproducible from SEED alone
# ---------------------------------------------------------------------------

def sample_windows(seed: int = SEED, n_windows: int = N_WINDOWS,
                   window_bp: int = WINDOW_BP) -> list[tuple[str, int, int]]:
    """Allocate windows proportional to chromosome length, place uniformly
    within autosome interiors (>= MIN_TELO_BP from each end)."""
    rng = np.random.default_rng(seed)

    # Proportional allocation, with a floor of 1 per chromosome.
    raw = [(c, n_windows * length / TOTAL_AUTOSOME_BP)
           for c, length in GRCH38_AUTOSOMES.items()]
    alloc = [(c, max(1, int(round(x)))) for c, x in raw]
    # Adjust to exact total by tweaking the largest chromosomes
    diff = n_windows - sum(n for _, n in alloc)
    alloc.sort(key=lambda kv: GRCH38_AUTOSOMES[kv[0]], reverse=True)
    i = 0
    while diff != 0 and i < len(alloc):
        c, n = alloc[i]
        if diff > 0:
            alloc[i] = (c, n + 1); diff -= 1
        elif n > 1:
            alloc[i] = (c, n - 1); diff += 1
        i = (i + 1) % len(alloc)
        if all(n == 1 for _, n in alloc) and diff < 0:
            break

    windows: list[tuple[str, int, int]] = []
    for chrom, n in alloc:
        length = GRCH38_AUTOSOMES[chrom]
        lo = MIN_TELO_BP
        hi = length - MIN_TELO_BP - window_bp
        if hi <= lo:
            continue
        starts = rng.integers(lo, hi, size=n)
        for s in sorted(int(s) for s in starts):
            windows.append((chrom, s, s + window_bp))
    return windows


# ---------------------------------------------------------------------------
# gnomAD variant fetch (reuse pattern from random_variant_control)
# ---------------------------------------------------------------------------

def fetch_region(chrom: str, start: int, end: int) -> list[dict]:
    """All variants in the region from gnomAD v3, no caller-side filtering yet."""
    last_exc = None
    for attempt in range(3):
        try:
            resp = requests.post(
                GNOMAD_API,
                json={"query": _GNOMAD_QUERY,
                      "variables": {"chrom": chrom, "start": start, "stop": end}},
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
            if resp.status_code == 429:
                wait = 30 * (attempt + 1)
                print(f"    rate-limited, sleeping {wait}s")
                time.sleep(wait)
                continue
            resp.raise_for_status()
            return (resp.json().get("data") or {}).get("region", {}).get("variants", []) or []
        except Exception as exc:
            last_exc = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"gnomAD fetch failed after 3 attempts: {last_exc}")


def common_snvs_from_region(raw_variants: list[dict],
                            tewhey_rsids: set[str]) -> tuple[list[dict], int]:
    """Filter to MAF>0.01 biallelic SNVs, drop Tewhey overlap. Return (variants,
    n_dropped_for_tewhey)."""
    out: list[dict] = []
    n_tewhey_drop = 0
    for v in raw_variants:
        af = (v.get("genome") or {}).get("af") or 0.0
        if af < MAF_THRESHOLD:
            continue
        vid = v.get("variant_id", "")
        parts = vid.split("-")
        if len(parts) != 4:
            continue
        chrom_str, pos_str, ref, alt = parts
        if not (VALID_BASE.match(ref) and VALID_BASE.match(alt)):
            continue  # SNV only
        rsid = v.get("rsid") or vid
        if rsid in tewhey_rsids:
            n_tewhey_drop += 1
            continue
        out.append({
            "rsid": rsid,
            "chrom": f"chr{chrom_str}",
            "pos": int(pos_str),
            "ref": ref,
            "alt": alt,
            "maf": float(af),
        })
    return out, n_tewhey_drop


def collect_variants(seed: int = SEED) -> list[dict]:
    """Run window sampling + gnomAD fetch + filtering. Reproducible from seed."""
    tewhey_rsids = set(pd.read_parquet(TEWHEY_PARQUET)["rsid"].dropna().astype(str))
    print(f"Tewhey rsids loaded: {len(tewhey_rsids):,}")

    windows = sample_windows(seed=seed)
    print(f"Sampled {len(windows)} windows ({WINDOW_BP/1000:.0f} kb each, "
          f"{sum(GRCH38_AUTOSOMES.values())/1e9:.2f} Gb genome, "
          f"{len(windows)*WINDOW_BP/1e6:.1f} Mb total)")

    pool: list[dict] = []
    seen_rsids: set[str] = set()
    n_tewhey_drops = 0
    for i, (chrom, start, end) in enumerate(windows, 1):
        if i > 1:
            time.sleep(2)  # gnomAD politeness
        try:
            raw = fetch_region(chrom, start, end)
        except Exception as exc:
            print(f"  WARN window {i}/{len(windows)} chr{chrom}:{start:,}-{end:,} "
                  f"failed: {exc}")
            continue
        kept, n_drop = common_snvs_from_region(raw, tewhey_rsids)
        n_tewhey_drops += n_drop
        # Dedupe across windows (rare, but possible if windows overlap)
        new = [v for v in kept if v["rsid"] not in seen_rsids]
        seen_rsids.update(v["rsid"] for v in new)
        pool.extend(new)
        if i % 10 == 0 or i == len(windows):
            print(f"  [{i}/{len(windows)}] chr{chrom} window → "
                  f"{len(kept)} kept ({n_drop} tewhey-dropped). "
                  f"Pool now {len(pool)}")

    print(f"\nVariant pool: {len(pool):,} unique common SNVs after MAF + Tewhey filter")
    print(f"  Tewhey-overlap drops across all windows: {n_tewhey_drops}")
    return pool


# ---------------------------------------------------------------------------
# Cache (per-row resume; thread-safe writes)
# ---------------------------------------------------------------------------

_DB_LOCK = threading.Lock()


def init_cache() -> None:
    conn = sqlite3.connect(CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scores (
            rsid                  TEXT PRIMARY KEY,
            chrom                 TEXT,
            pos                   INTEGER,
            ref                   TEXT,
            alt                   TEXT,
            maf                   REAL,
            raw_max_signed_delta  REAL,
            single_track_quantile REAL,
            raw_mean_signed_delta   REAL,
            raw_median_signed_delta REAL,
            n_expr_tracks         INTEGER,
            error                 TEXT,
            scored_at             INTEGER
        )
    """)
    # Idempotent backfill of mean/median columns on pre-existing caches.
    existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(scores)").fetchall()}
    for col, sqltype in (("raw_mean_signed_delta", "REAL"),
                         ("raw_median_signed_delta", "REAL")):
        if col not in existing_cols:
            conn.execute(f"ALTER TABLE scores ADD COLUMN {col} {sqltype}")
    conn.commit()
    conn.close()


def save_row(v: dict, entry: dict) -> None:
    with _DB_LOCK:
        conn = sqlite3.connect(CACHE_DB)
        conn.execute(
            """INSERT OR REPLACE INTO scores
               (rsid,chrom,pos,ref,alt,maf,
                raw_max_signed_delta,single_track_quantile,
                raw_mean_signed_delta,raw_median_signed_delta,
                n_expr_tracks,error,scored_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (v["rsid"], v.get("chrom"), v.get("pos"), v.get("ref"),
             v.get("alt"), v.get("maf"),
             entry.get("raw_max_signed_delta"),
             entry.get("single_track_quantile"),
             entry.get("raw_mean_signed_delta"),
             entry.get("raw_median_signed_delta"),
             entry.get("n_expr_tracks"),
             entry.get("error"),
             int(time.time())),
        )
        conn.commit()
        conn.close()


def load_cache() -> dict[str, dict]:
    if not CACHE_DB.exists():
        return {}
    conn = sqlite3.connect(CACHE_DB)
    rows = conn.execute(
        "SELECT rsid, raw_max_signed_delta, single_track_quantile, "
        "raw_mean_signed_delta, raw_median_signed_delta, "
        "n_expr_tracks, error FROM scores"
    ).fetchall()
    conn.close()
    return {r[0]: {"raw_max_signed_delta": r[1], "single_track_quantile": r[2],
                   "raw_mean_signed_delta": r[3], "raw_median_signed_delta": r[4],
                   "n_expr_tracks": r[5], "error": r[6]} for r in rows}


def needs_rescoring(entry: dict | None) -> bool:
    """Re-score if the row is missing (not cached), or has a successful
    raw_max_signed_delta but no raw_mean_signed_delta (the new column).
    Skip rows that already have an error — those failed for a reason and
    will fail again."""
    if entry is None:
        return True
    if entry.get("error") is not None:
        return False  # already tried, failed; leave alone
    if entry.get("raw_mean_signed_delta") is None:
        return True
    return False


# ---------------------------------------------------------------------------
# K562 profile + raw delta + matched single-track quantile
# ---------------------------------------------------------------------------

def load_k562_profile():
    from scoring.tissue_config import TissueProfile
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    k562 = cfg.get("tissue_profiles", {}).get("tewhey_k562", {})
    return TissueProfile(
        display_name="K562 / Blood-lineage (matched calibration)",
        biosample_keywords=k562.get("biosample_keywords", []),
        gtex_keywords=k562.get("gtex_keywords", []),
    )


def extract_max_signed_with_quantile(tidy_df: pd.DataFrame, profile) -> dict:
    """Per-variant K562 expression stats from one tidy_df.

    Returns:
      raw_max_signed_delta    — sign(argmax|raw|) * max|raw|  (peak-track)
      single_track_quantile   — quantile_score of the same peak-track row
      raw_mean_signed_delta   — arithmetic mean of signed raw across tracks
      raw_median_signed_delta — median of signed raw across tracks
      n_expr_tracks           — number of K562 expression tracks contributing
    """
    from scoring.tissue_config import filter_tracks

    empty = {"raw_max_signed_delta": None, "single_track_quantile": None,
             "raw_mean_signed_delta": None, "raw_median_signed_delta": None,
             "n_expr_tracks": 0}

    if tidy_df is None or tidy_df.empty:
        return empty

    filtered = filter_tracks(tidy_df, profile)
    ot = filtered["output_type"].astype(str).str.replace("OutputType.", "", regex=False)
    expr = filtered[ot.isin(EXPRESSION_OUTPUT_TYPES)]
    if expr.empty:
        return empty

    raw = expr["raw_score"].dropna()
    if raw.empty:
        return empty

    idx_max = raw.abs().idxmax()
    max_signed = float(raw.loc[idx_max])
    sq = expr.loc[idx_max, "quantile_score"] if "quantile_score" in expr.columns else None
    single_track_q = float(sq) if sq is not None and pd.notna(sq) else None

    return {
        "raw_max_signed_delta":    max_signed,
        "single_track_quantile":   single_track_q,
        "raw_mean_signed_delta":   float(raw.mean()),
        "raw_median_signed_delta": float(raw.median()),
        "n_expr_tracks":           int(len(raw)),
    }


# ---------------------------------------------------------------------------
# Per-variant scoring with timeout + rate-limit guard
# ---------------------------------------------------------------------------

class RateLimitError(RuntimeError):
    """Raised to bail the entire run on suspected gating."""


def score_one(worker, vi, profile) -> dict:
    """60-second timeout via subprocess worker. Raises RateLimitError on
    suspected rate-limiting so the main thread can stop the run cleanly.

    The worker's .score() is uniform across ScoreWorker (serial) and
    ScoreWorkerPool (parallel) — both expose the same call signature.
    """
    result = worker.score(vi, profile, timeout=VARIANT_TIMEOUT_SECS)
    err = result.get("error") if isinstance(result, dict) else None

    if err == "api_timeout":
        return {"error": "api_timeout"}

    if err and "composite_score" not in result:
        # Worker-level exception path (rare — score_single_variant catches
        # most things and returns a full-shape error dict instead).
        msg = err.lower()
        if any(k in msg for k in ("rate", "quota", "limit", "429")):
            raise RateLimitError(err)
        if "sequence length" in msg and "not supported" in msg:
            return {"error": "insufficient_flanking_sequence"}
        return {"error": err}

    if err:
        return {"error": err}

    return extract_max_signed_with_quantile(result.get("tidy_df"), profile)


def score_variant_dispatch(worker, profile, v: dict) -> tuple[dict, dict, float]:
    """Worker entry. Returns (input_v, entry_dict, elapsed_seconds).

    `worker` is a ScoreWorker (serial) or ScoreWorkerPool (parallel) — both
    expose a thread-safe .score(vi, profile, timeout=...) method.
    """
    from scoring.composite import VariantInput

    rsid = v["rsid"]
    if not (VALID_BASE.match(v["ref"]) and VALID_BASE.match(v["alt"])):
        return v, {"error": "invalid_bases"}, 0.0

    try:
        vi = VariantInput(rsid=rsid, chrom=v["chrom"], pos=int(v["pos"]),
                          ref=v["ref"], alt=v["alt"], maf=v.get("maf"))
    except Exception as exc:
        return v, {"error": f"input_error:{exc}"}, 0.0

    t0 = time.monotonic()
    entry = score_one(worker, vi, profile)
    return v, entry, time.monotonic() - t0


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def _line(i: int, total: int, v: dict, entry: dict, elapsed: float,
          counters: dict) -> None:
    err = entry.get("error")
    if err:
        flag = f" ERROR={err}"
    else:
        mean = entry.get("raw_mean_signed_delta")
        med  = entry.get("raw_median_signed_delta")
        mean_s = f" mean={mean:+.4f}" if mean is not None else ""
        med_s  = f" med={med:+.4f}"   if med  is not None else ""
        flag = (f" raw_max={entry['raw_max_signed_delta']:+.4f}"
                f" q={entry['single_track_quantile']:+.4f}"
                f"{mean_s}{med_s}"
                f" n_tracks={entry['n_expr_tracks']}")
    print(f"  [{i+1}/{total}] {v['rsid']}{flag}  ({elapsed:.1f}s)  "
          f"ok={counters['ok']} err={counters['err']} timeout={counters['timeout']}")


def run_serial(variants: list[dict], worker, profile,
               max_n: int | None = None) -> dict:
    """Pilot path: serial loop, prints per-variant lines, no DB lock contention."""
    counters = {"ok": 0, "err": 0, "timeout": 0}
    total = len(variants) if max_n is None else min(max_n, len(variants))
    for i, v in enumerate(variants[:total]):
        try:
            v, entry, elapsed = score_variant_dispatch(worker, profile, v)
        except RateLimitError as exc:
            print(f"!! Rate limit hit on {v['rsid']}: {exc}\n   Stopping.")
            sys.exit(1)
        save_row(v, entry)
        if entry.get("error") == "api_timeout":
            counters["timeout"] += 1
        elif entry.get("error"):
            counters["err"] += 1
        else:
            counters["ok"] += 1
        _line(i, total, v, entry, elapsed, counters)
    return counters


def run_parallel(variants: list[dict], pool, profile) -> dict:
    """Full path: 4-way ThreadPool over variants, dispatching each call to the
    subprocess `pool`. Each task uses pool.score() which checks out a free
    subprocess worker; on timeout that worker is killed + respawned. Aborts on
    any RateLimitError."""
    counters = {"ok": 0, "err": 0, "timeout": 0}
    total = len(variants)
    saturated = []  # raw_max_signed_delta values for the 200-variant stop check

    print(f"Launching {PARALLEL_WORKERS}-way parallel scoring over {total} variants...")
    started = time.monotonic()
    rate_limit_hit = threading.Event()
    rate_limit_msg: list[str] = []

    def _task(v):
        if rate_limit_hit.is_set():
            return v, {"error": "aborted_rate_limit"}, 0.0
        try:
            return score_variant_dispatch(pool, profile, v)
        except RateLimitError as exc:
            rate_limit_msg.append(str(exc))
            rate_limit_hit.set()
            return v, {"error": "rate_limit_hit"}, 0.0

    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=PARALLEL_WORKERS) as ex:
        futures = {ex.submit(_task, v): v for v in variants}
        for fut in concurrent.futures.as_completed(futures):
            v, entry, elapsed = fut.result()
            save_row(v, entry)
            err = entry.get("error")
            if err == "api_timeout":
                counters["timeout"] += 1
            elif err:
                counters["err"] += 1
            else:
                counters["ok"] += 1
                if entry.get("raw_max_signed_delta") is not None:
                    saturated.append(abs(entry["raw_max_signed_delta"]))

            completed += 1
            if completed % 25 == 0 or completed == total:
                rate = completed / max(time.monotonic() - started, 1e-9)
                eta_min = (total - completed) / max(rate, 1e-9) / 60
                print(f"  [{completed}/{total}] "
                      f"ok={counters['ok']} err={counters['err']} "
                      f"timeout={counters['timeout']}  "
                      f"rate={rate:.2f}/s  eta={eta_min:.1f}min")

            # Diagnostic stop: common variants saturating defeats the experiment.
            if completed == SAT_CHECK_AFTER and len(saturated) >= SAT_CHECK_AFTER // 2:
                sat_frac = float(np.mean(np.array(saturated) > 0.9))
                print(f"  Sanity check @ {completed} variants: "
                      f"|raw_max| > 0.9 fraction = {sat_frac:.1%}")
                if sat_frac > SAT_HALT_FRACTION:
                    print(f"!! Common-variant saturation {sat_frac:.1%} > "
                          f"{SAT_HALT_FRACTION:.0%} threshold — halting. "
                          f"Diagnose upstream before continuing.")
                    rate_limit_hit.set()  # reuse abort flag

            if rate_limit_hit.is_set() and rate_limit_msg:
                # Cancel pending work
                for f in futures:
                    if not f.done():
                        f.cancel()
                break

    if rate_limit_msg:
        print(f"!! Rate-limit / abort: {rate_limit_msg[0]}")

    return counters


# ---------------------------------------------------------------------------
# Post-processing: null parquet + histogram
# ---------------------------------------------------------------------------

def build_null_and_figure() -> None:
    if not CACHE_DB.exists():
        print("No cache yet — nothing to post-process.")
        return

    conn = sqlite3.connect(CACHE_DB)
    df = pd.read_sql_query(
        "SELECT rsid, chrom, pos, ref, alt, maf, "
        "raw_max_signed_delta, single_track_quantile, "
        "raw_mean_signed_delta, raw_median_signed_delta, "
        "n_expr_tracks, error "
        "FROM scores",
        conn,
    )
    conn.close()

    n_total = len(df)
    n_err   = df["error"].notna().sum()
    clean   = df[df["error"].isna() & df["raw_max_signed_delta"].notna()].copy()
    clean   = clean.sort_values("raw_max_signed_delta").reset_index(drop=True)
    clean.to_parquet(NULL_PARQUET, index=False)
    print(f"\nNull written → {NULL_PARQUET}  ({len(clean):,} clean rows / "
          f"{n_total:,} total / {n_err:,} errors)")

    # Summary stats — three aggregation recipes plus the published peak-track quantile
    x_max  = clean["raw_max_signed_delta"].to_numpy(float)
    x_mean = clean["raw_mean_signed_delta"].dropna().to_numpy(float)
    x_med  = clean["raw_median_signed_delta"].dropna().to_numpy(float)
    q      = clean["single_track_quantile"].dropna().to_numpy(float)

    def _print_block(name: str, x: np.ndarray) -> None:
        if len(x) == 0:
            print(f"\n=== {name}: no data ==="); return
        print(f"\n=== {name}  (n = {len(x):,}) ===")
        print(f"  mean       = {x.mean():+.4f}")
        print(f"  median     = {np.median(x):+.4f}")
        print(f"  std        = {x.std():.4f}")
        print(f"  IQR        = [{np.percentile(x,25):+.4f}, {np.percentile(x,75):+.4f}]")
        print(f"  range      = [{x.min():+.4f}, {x.max():+.4f}]")
        print(f"  |·|>0.5 fraction = {(np.abs(x)>0.5).mean():.3%}")
        print(f"  |·|>0.9 fraction = {(np.abs(x)>0.9).mean():.3%}")
        print(f"  exactly ±1       = {((x==1.0) | (x==-1.0)).mean():.3%}")

    _print_block("Matched null — raw_max_signed_delta (current pipeline)", x_max)
    _print_block("Matched null — raw_mean_signed_delta (alternative)",    x_mean)
    _print_block("Matched null — raw_median_signed_delta (alternative)",  x_med)
    if len(q) > 0:
        print(f"\n=== Published single_track_quantile (peak-track, for reference) ===")
        print(f"  n               = {len(q):,}")
        print(f"  |q|>0.9 fraction = {(np.abs(q)>0.9).mean():.3%}")

    # Histogram (uses the max-aggregation null — primary headline distribution)
    FIG_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = np.linspace(-max(abs(x_max.min()), abs(x_max.max()), 1.0),
                        max(abs(x_max.min()), abs(x_max.max()), 1.0), 80)

    ax = axes[0]
    ax.hist(x_max, bins=bins, color="#4575b4", edgecolor="white", alpha=0.85)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(+0.9, color="#d6604d", ls="--", lw=0.7, label="|Δ|=0.9")
    ax.axvline(-0.9, color="#d6604d", ls="--", lw=0.7)
    ax.set_xlabel("raw_max_signed_delta (max signed K562/blood expression Δ)")
    ax.set_ylabel("# common variants")
    ax.set_title(f"Matched calibration null  (n={len(x_max):,})")
    ann = (f"mean={x_max.mean():+.3f}\nstd={x_max.std():.3f}\n"
           f"|Δ|>0.9: {(np.abs(x_max)>0.9).mean():.1%}")
    ax.text(0.97, 0.97, ann, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.9))
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", lw=0.3, alpha=0.5)

    if len(q) > 0:
        ax = axes[1]
        ax.hist(q, bins=80, color="#91bfdb", edgecolor="white", alpha=0.85)
        ax.axvline(0, color="grey", ls=":", lw=0.8)
        ax.axvline(+0.9, color="#d6604d", ls="--", lw=0.7)
        ax.axvline(-0.9, color="#d6604d", ls="--", lw=0.7)
        ax.set_xlabel("single_track_quantile (published peak-track quantile)")
        ax.set_ylabel("# common variants")
        ax.set_title("Published single-track quantile (same variants, peak track)")
        ann2 = (f"|q|>0.9: {(np.abs(q)>0.9).mean():.1%}\n"
                f"q≈±1 (saturated): {((np.abs(q)>=0.999)).mean():.1%}")
        ax.text(0.97, 0.97, ann2, transform=ax.transAxes, ha="right", va="top",
                fontsize=9, family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.9))
        ax.grid(axis="y", lw=0.3, alpha=0.5)

    fig.suptitle(
        "Matched-statistic null on common autosomal variants  "
        "(MAF>0.01, K562 expression, peak track)",
        fontsize=10.5, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(HIST_FIG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Histogram → {HIST_FIG}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true",
                        help="Score 5 variants serial as a smoke test")
    parser.add_argument("--post", action="store_true",
                        help="Skip scoring; just rebuild null parquet + histogram")
    parser.add_argument("--seed", type=int, default=SEED,
                        help="Window-sampling RNG seed (default: 2026)")
    parser.add_argument("--cap", type=int, default=6000,
                        help="Hard cap on candidates to score (default: 6000)")
    args = parser.parse_args()

    if args.post:
        build_null_and_figure()
        return

    api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set"); sys.exit(1)

    init_cache()
    cache = load_cache()
    print(f"Existing cache: {len(cache):,} rows")

    print("\n=== Step 1: collect variants ===")
    pool = collect_variants(seed=args.seed)

    # Deterministic cap by seed
    if len(pool) > args.cap:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(len(pool), size=args.cap, replace=False)
        pool = [pool[i] for i in sorted(int(i) for i in idx)]
        print(f"Capped pool to {len(pool):,} variants (seed-deterministic subset)")

    to_score = [v for v in pool if needs_rescoring(cache.get(v["rsid"]))]
    n_full_cached  = sum(1 for v in pool
                        if v["rsid"] in cache
                        and not needs_rescoring(cache[v["rsid"]])
                        and cache[v["rsid"]].get("error") is None)
    n_prior_errors = sum(1 for v in pool
                        if v["rsid"] in cache
                        and cache[v["rsid"]].get("error") is not None)
    n_needs_mean   = sum(1 for v in pool
                        if v["rsid"] in cache
                        and cache[v["rsid"]].get("error") is None
                        and cache[v["rsid"]].get("raw_mean_signed_delta") is None)
    print(f"  cached complete (max+mean+median): {n_full_cached}")
    print(f"  cached but missing mean/median:    {n_needs_mean}  (will re-score)")
    print(f"  prior errors (skipped):            {n_prior_errors}")
    print(f"  uncached:                          {len(to_score) - n_needs_mean}")
    print(f"  → {len(to_score)} total to score")

    if not to_score and not args.pilot:
        print("Nothing to score; running post-processing.")
        build_null_and_figure()
        return

    print("\n=== Step 2: load K562 profile + spawn scoring workers ===")
    from scoring.score_worker import ScoreWorker, ScoreWorkerPool
    profile = load_k562_profile()
    print("K562 profile loaded.")

    if args.pilot:
        print(f"\n=== PILOT: {PILOT_N} variants serial ===")
        # Pilot picks the first PILOT_N variants from the new (uncached) list
        pilot_set = to_score[:PILOT_N] if to_score else pool[:PILOT_N]
        worker = ScoreWorker(api_key)
        try:
            counters = run_serial(pilot_set, worker, profile, max_n=PILOT_N)
        finally:
            worker.close()
        print(f"\nPilot result: ok={counters['ok']} err={counters['err']} "
              f"timeout={counters['timeout']}")
        if counters['ok'] == 0:
            print("!! PILOT FAILED — no successful scores. Diagnose before parallel run.")
            sys.exit(2)
        if counters['err'] + counters['timeout'] > PILOT_N // 2:
            print(f"!! PILOT degraded — {counters['err']} err + "
                  f"{counters['timeout']} timeout out of {PILOT_N}. "
                  f"Investigate before parallel.")
            sys.exit(2)
        print("\nPilot looks healthy. Re-run without --pilot to launch the parallel pass.")
        return

    print(f"\n=== Step 3: parallel scoring (workers={PARALLEL_WORKERS}) ===")
    worker_pool = ScoreWorkerPool(PARALLEL_WORKERS, api_key)
    try:
        counters = run_parallel(to_score, worker_pool, profile)
    finally:
        worker_pool.close()
    print(f"\nParallel run complete: ok={counters['ok']} err={counters['err']} "
          f"timeout={counters['timeout']}")

    print("\n=== Step 4: post-process ===")
    build_null_and_figure()


if __name__ == "__main__":
    main()
