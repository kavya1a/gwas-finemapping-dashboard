"""Batch-score preloaded variants through AlphaGenome; cache to scored_variants.db.

Design principles:
  - Incremental writes: each variant is committed to DB immediately after scoring.
    If the process dies mid-run, completed work survives and is skipped on re-run.
  - Per-variant timeout: each AlphaGenome call is wrapped in a 60-second
    ThreadPoolExecutor timeout. On timeout the variant is marked api_timeout in
    the preloaded DB, logged to scoring_timeouts.log, and the run continues.
  - Sequential per-variant scoring (not batched) to support both of the above.
"""

from __future__ import annotations

import concurrent.futures
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

import yaml
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

DIR = Path(__file__).parent
PRELOADED_DB  = DIR / "preloaded_variants.db"
SCORED_DB     = DIR / "scored_variants.db"
BLOCKERS_FILE = DIR / "OVERNIGHT_BLOCKERS.md"
TIMEOUTS_LOG  = DIR / "scoring_timeouts.log"
CONFIG_PATH   = DIR / "config.yaml"

# Load tunable parameters from config.yaml (fallback to hard-coded defaults)
_cfg = yaml.safe_load(CONFIG_PATH.read_text()) if CONFIG_PATH.exists() else {}
CONFIG_VERSION       = _cfg.get("config_version", "unknown")
VARIANT_TIMEOUT_SECS = _cfg.get("pipeline", {}).get("variant_timeout_secs", 60)
MAX_TIMEOUT_FRACTION = _cfg.get("pipeline", {}).get("max_timeout_fraction", 0.10)
MIN_SUCCESS_FRACTION = _cfg.get("pipeline", {}).get("min_success_fraction", 0.50)

DISEASES = ["alzheimers", "t2d", "schizophrenia", "parkinsons"]

_VALID_BASES = re.compile(r"^[ACGTNacgtn]+$")


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def _init_scored_db(conn: sqlite3.Connection) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scores (
            rsid                TEXT NOT NULL,
            disease             TEXT NOT NULL,
            chrom               TEXT,
            pos                 INTEGER,
            ref                 TEXT,
            alt                 TEXT,
            maf                 REAL,
            composite_score     REAL,
            pip_weighted_score  REAL,
            pip                 REAL,
            rare_variant_caution INTEGER DEFAULT 0,
            error               TEXT,
            rank                INTEGER,
            scored_at           INTEGER DEFAULT 0,
            PRIMARY KEY (rsid, disease)
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS modality_scores (
            rsid                  TEXT NOT NULL,
            disease               TEXT NOT NULL,
            modality              TEXT NOT NULL,
            max_abs_score         REAL,
            frac_above_threshold  REAL,
            signed_max_score      REAL,
            PRIMARY KEY (rsid, disease, modality)
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS run_metadata (
            run_id          TEXT NOT NULL,
            config_version  TEXT,
            disease         TEXT NOT NULL,
            started_at      INTEGER,
            completed_at    INTEGER,
            n_scored        INTEGER,
            n_ok            INTEGER,
            n_error         INTEGER,
            n_timeout       INTEGER,
            PRIMARY KEY (run_id, disease)
        )
    """)
    # Migrate existing DBs that predate signed_max_score / run_metadata
    cols = [r[1] for r in conn.execute("PRAGMA table_info(modality_scores)").fetchall()]
    if "signed_max_score" not in cols:
        conn.execute("ALTER TABLE modality_scores ADD COLUMN signed_max_score REAL")
    conn.commit()


def _write_blocker(msg: str) -> None:
    with open(BLOCKERS_FILE, "a") as f:
        f.write(f"\n## BLOCKER [{time.strftime('%Y-%m-%d %H:%M:%S')}] batch_score\n{msg}\n")
    print(f"!! BLOCKER: {msg}")


def _log_timeout(rsid: str, disease: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(TIMEOUTS_LOG, "a") as f:
        f.write(f"{ts}\t{disease}\t{rsid}\n")


def _mark_preloaded(rsid: str, disease: str, reason: str) -> None:
    try:
        conn = sqlite3.connect(PRELOADED_DB)
        conn.execute(
            "UPDATE variants SET scoring_skipped_reason = ? "
            "WHERE rsid = ? AND disease = ? AND scoring_skipped_reason IS NULL",
            (reason, rsid, disease),
        )
        conn.commit()
        conn.close()
    except Exception:
        pass


def _save_single_result(result: dict, disease: str) -> None:
    """Write one variant result to scored_variants.db immediately (incremental)."""
    conn = sqlite3.connect(SCORED_DB)
    _init_scored_db(conn)
    try:
        rsid = result["rsid"]
        if result.get("error"):
            _mark_preloaded(rsid, disease, "api_error")

        conn.execute(
            """INSERT OR REPLACE INTO scores
               (rsid,disease,chrom,pos,ref,alt,maf,
                composite_score,pip_weighted_score,pip,
                rare_variant_caution,error,scored_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                rsid, disease,
                result.get("chrom"), result.get("pos"),
                result.get("ref"), result.get("alt"),
                result.get("maf"),
                result.get("composite_score"),
                result.get("pip_weighted_score"),
                result.get("pip"),
                int(bool(result.get("rare_variant_caution", False))),
                result.get("error"),
                int(time.time()),
            ),
        )

        mod_df = result.get("modality_breakdown")
        if mod_df is not None and not mod_df.empty:
            for _, mrow in mod_df.iterrows():
                conn.execute(
                    """INSERT OR REPLACE INTO modality_scores
                       (rsid,disease,modality,max_abs_score,frac_above_threshold,signed_max_score)
                       VALUES (?,?,?,?,?,?)""",
                    (rsid, disease,
                     mrow.get("modality"),
                     mrow.get("max_abs_score"),
                     mrow.get("frac_above_threshold"),
                     mrow.get("signed_max_score")),
                )
        conn.commit()
    finally:
        conn.close()


def _load_preloaded(disease: str) -> list[dict]:
    if not PRELOADED_DB.exists():
        return []
    conn = sqlite3.connect(PRELOADED_DB)
    try:
        rows = conn.execute(
            "SELECT rsid, chrom, pos, ref, alt, maf, p_value "
            "FROM variants "
            "WHERE disease = ? AND scoring_skipped_reason IS NULL",
            (disease,),
        ).fetchall()
    finally:
        conn.close()
    return [
        {"rsid": r[0], "chrom": r[1], "pos": r[2], "ref": r[3],
         "alt": r[4], "maf": r[5], "p_value": r[6]}
        for r in rows
    ]


def _already_scored(disease: str) -> set[str]:
    if not SCORED_DB.exists():
        return set()
    conn = sqlite3.connect(SCORED_DB)
    _init_scored_db(conn)
    try:
        rows = conn.execute(
            "SELECT rsid FROM scores WHERE disease = ?", (disease,)
        ).fetchall()
    finally:
        conn.close()
    return {r[0] for r in rows}


# ---------------------------------------------------------------------------
# Per-variant scoring with timeout
# ---------------------------------------------------------------------------

def _score_one_with_timeout(
    model,
    variant_input,
    tissue_profile,
    rsid: str,
    disease: str,
) -> dict | None:
    """Score a single variant with a 60-second wall-clock timeout.

    Returns result dict on success, None on timeout (already logged + marked).
    Creates a fresh ThreadPoolExecutor per call so a timed-out thread never
    blocks a subsequent submission.
    """
    from scoring.composite import score_single_variant

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    future = executor.submit(score_single_variant, model, variant_input, tissue_profile)

    try:
        return future.result(timeout=VARIANT_TIMEOUT_SECS)
    except concurrent.futures.TimeoutError:
        _log_timeout(rsid, disease)
        _mark_preloaded(rsid, disease, "api_timeout")
        print(f"  TIMEOUT ({VARIANT_TIMEOUT_SECS}s): {rsid} — logged, skipping")
        return None
    except Exception as exc:
        err_str = str(exc).lower()
        if any(k in err_str for k in ("rate", "quota", "limit", "429")):
            _write_blocker(
                f"AlphaGenome API rate limit hit scoring {rsid} ({disease}): {exc}\n"
                "Stopping to avoid further charges."
            )
            executor.shutdown(wait=False)
            sys.exit(1)
        # Variants within ~524kb of a chromosome end cannot be scored: the 1Mb
        # context model requires at least 524,288bp of flanking sequence on each
        # side, and positions closer than that to a telomere produce an
        # unsupported sequence-length error.  Tag distinctly so the yield report
        # can separate these from true API failures.
        if "sequence length" in err_str and "not supported" in err_str:
            _mark_preloaded(rsid, disease, "insufficient_flanking_sequence")
            print(f"  Flanking seq too short {rsid}: {exc}")
            return {
                "rsid": rsid, "chrom": None, "pos": None, "ref": None, "alt": None,
                "maf": None, "composite_score": None, "pip_weighted_score": None,
                "pip": None, "rare_variant_caution": False,
                "error": str(exc), "modality_breakdown": None, "tidy_df": None,
            }
        # Other API errors: record and continue
        _mark_preloaded(rsid, disease, "api_error")
        print(f"  API error {rsid}: {exc}")
        return {
            "rsid": rsid, "chrom": None, "pos": None, "ref": None, "alt": None,
            "maf": None, "composite_score": None, "pip_weighted_score": None,
            "pip": None, "rare_variant_caution": False,
            "error": str(exc), "modality_breakdown": None, "tidy_df": None,
        }
    finally:
        executor.shutdown(wait=False)


# ---------------------------------------------------------------------------
# Per-disease scoring loop
# ---------------------------------------------------------------------------

def _write_run_metadata(
    run_id: str, disease: str, started_at: int, counters: dict
) -> None:
    conn = sqlite3.connect(SCORED_DB)
    _init_scored_db(conn)
    conn.execute(
        """INSERT OR REPLACE INTO run_metadata
           (run_id,config_version,disease,started_at,completed_at,
            n_scored,n_ok,n_error,n_timeout)
           VALUES (?,?,?,?,?,?,?,?,?)""",
        (
            run_id, CONFIG_VERSION, disease, started_at, int(time.time()),
            counters["ok"] + counters["error"],
            counters["ok"], counters["error"], counters["timeout"],
        ),
    )
    conn.commit()
    conn.close()


def _score_disease(disease: str, to_score: list[dict], api_key: str, run_id: str) -> dict:
    """Score all variants for one disease, writing incrementally. Returns counters."""
    from alphagenome.models import dna_client
    from scoring.composite import VariantInput, _build_interval
    from scoring.tissue_config import get_profile

    profile = get_profile(disease)
    model = dna_client.create(api_key)

    counters = {"ok": 0, "error": 0, "timeout": 0, "skipped": 0}
    started_at = int(time.time())

    for i, v in enumerate(to_score):
        rsid = v["rsid"]
        chrom = v["chrom"] if v["chrom"].startswith("chr") else f"chr{v['chrom']}"

        # Belt-and-suspenders allele check
        if not (v.get("ref") and v.get("alt") and
                _VALID_BASES.match(v["ref"]) and _VALID_BASES.match(v["alt"])):
            _mark_preloaded(rsid, disease, "indel_not_supported")
            counters["skipped"] += 1
            continue

        variant_input = VariantInput(
            rsid=rsid,
            chrom=chrom,
            pos=int(v["pos"]),
            ref=v["ref"],
            alt=v["alt"],
            maf=v.get("maf"),
            p_value=v.get("p_value"),
        )

        t0 = time.monotonic()
        result = _score_one_with_timeout(model, variant_input, profile, rsid, disease)
        elapsed = time.monotonic() - t0

        if result is None:
            counters["timeout"] += 1
            continue

        _save_single_result(result, disease)

        if result.get("error"):
            counters["error"] += 1
            flag = " ERROR"
        else:
            counters["ok"] += 1
            flag = f" composite={result['composite_score']:.4f}"

        print(f"  [{i+1}/{len(to_score)}] {rsid}{flag}  ({elapsed:.1f}s)")

    _write_run_metadata(run_id, disease, started_at, counters)
    return counters


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def batch_score_all(target_diseases: list[str] | None = None) -> None:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not api_key:
        _write_blocker("ALPHAGENOME_API_KEY not set.")
        sys.exit(1)

    conn = sqlite3.connect(SCORED_DB)
    _init_scored_db(conn)
    conn.close()

    run_id = f"run_{int(time.time())}"
    diseases = target_diseases or DISEASES
    total_timeouts = 0
    total_variants = 0

    for disease in diseases:
        print(f"\n{'='*60}")
        print(f"Disease: {disease}")

        preloaded = _load_preloaded(disease)
        already = _already_scored(disease)
        to_score = [v for v in preloaded if v["rsid"] not in already]

        print(f"  {len(preloaded)} scoreable, {len(already)} already in DB, {len(to_score)} to score")

        if not to_score:
            print("  All scored — skipping")
            continue

        counters = _score_disease(disease, to_score, api_key, run_id)

        ok, err, timeout, skipped = (
            counters["ok"], counters["error"], counters["timeout"], counters["skipped"]
        )
        total_timeouts += timeout
        total_variants += len(to_score)

        print(f"\n  {disease}: {ok} ok / {err} api_error / {timeout} timeout / {skipped} skipped")

        # Per-disease stop conditions
        attempted = ok + err + timeout
        if attempted > 0:
            success_rate = ok / attempted
            if success_rate < MIN_SUCCESS_FRACTION:
                _write_blocker(
                    f"{disease}: success rate {success_rate:.1%} ({ok}/{attempted}) "
                    f"is below the 50% threshold."
                )
                sys.exit(1)

        if timeout > 0:
            timeout_rate = timeout / max(attempted, 1)
            if total_timeouts / max(total_variants, 1) > MAX_TIMEOUT_FRACTION:
                _write_blocker(
                    f"Cumulative timeout rate {total_timeouts}/{total_variants} "
                    f"({100*total_timeouts/total_variants:.1f}%) exceeds 10% threshold — "
                    "systemic API issue."
                )
                sys.exit(1)

    print(f"\n{'='*60}")
    print(f"All diseases done. Total timeouts: {total_timeouts}/{total_variants}")


if __name__ == "__main__":
    targets = sys.argv[1:] if len(sys.argv) > 1 else None
    batch_score_all(targets)
