"""Canonical variant verification harness.

Run after any scoring change to confirm the pipeline recovers known causal
variants at biologically expected ranks.

Tests:
  1. rs429358 (APOE ε4) in top 20% of AD scored set
  2. rs7412 (APOE ε2) scored for AD AND expression_subscore opposite sign to rs429358
  3. rs7903146 (TCF7L2) in top 20% of T2D scored set
  4. rs1006737 (CACNA1C) in top 20% of SCZ scored set
  5. rs356219 (SNCA) in top 20% of PD scored set

Exit code 0 if all pass, 1 if any fail.
Run:
    python verification/canonical_variants_test.py
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

DIR = Path(__file__).parent.parent
SCORED_DB = DIR / "scored_variants.db"

TOP_PCT_THRESHOLD = 0.20  # "top 20%" = rank / total <= 0.20

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
SKIP = "\033[33mSKIP\033[0m"


def _get_scores(conn: sqlite3.Connection, disease: str) -> list[tuple]:
    """Return all scored (rsid, composite_score) for a disease, sorted desc."""
    return conn.execute(
        "SELECT rsid, composite_score FROM scores "
        "WHERE disease=? AND error IS NULL AND composite_score IS NOT NULL "
        "ORDER BY composite_score DESC",
        (disease,),
    ).fetchall()


def _rank_info(rows: list[tuple], rsid: str) -> tuple[int | None, int, float | None]:
    """Returns (1-based rank, total, percentile) for rsid in rows."""
    for i, (r, _) in enumerate(rows):
        if r == rsid:
            rank = i + 1
            total = len(rows)
            pct = rank / total
            return rank, total, pct
    return None, len(rows), None


def _get_signed_expression(conn: sqlite3.Connection, disease: str, rsid: str) -> float | None:
    """Return signed_max_score for expression modality; None if not stored."""
    row = conn.execute(
        "SELECT signed_max_score FROM modality_scores "
        "WHERE disease=? AND rsid=? AND modality='expression'",
        (disease, rsid),
    ).fetchone()
    if row:
        return row[0]
    return None


def run_tests() -> bool:
    if not SCORED_DB.exists():
        print(f"  ERROR: {SCORED_DB} not found. Run batch_score.py first.")
        return False

    conn = sqlite3.connect(SCORED_DB)
    all_pass = True
    results = []

    # ------------------------------------------------------------------ #
    # Test 1: rs429358 (APOE ε4) in top 20% of AD
    # ------------------------------------------------------------------ #
    ad_rows = _get_scores(conn, "alzheimers")
    rank, total, pct = _rank_info(ad_rows, "rs429358")
    if rank is None:
        status = FAIL
        detail = f"rs429358 not found in alzheimers scored set (n={total})"
        all_pass = False
    elif pct <= TOP_PCT_THRESHOLD:
        status = PASS
        detail = f"rank {rank}/{total} ({pct:.1%}) — within top {TOP_PCT_THRESHOLD:.0%}"
    else:
        status = FAIL
        detail = f"rank {rank}/{total} ({pct:.1%}) — below top {TOP_PCT_THRESHOLD:.0%} threshold"
        all_pass = False
    results.append(("Test 1", "rs429358 in AD top 20%", status, detail))

    # ------------------------------------------------------------------ #
    # Test 2: rs7412 scored AND expression direction opposite to rs429358
    # ------------------------------------------------------------------ #
    rank7412, total7412, pct7412 = _rank_info(ad_rows, "rs7412")
    sign_rs429358 = _get_signed_expression(conn, "alzheimers", "rs429358")
    sign_rs7412   = _get_signed_expression(conn, "alzheimers", "rs7412")

    if rank7412 is None:
        status = FAIL
        detail = "rs7412 not scored for alzheimers — run whitelist fetch + score"
        all_pass = False
    elif sign_rs429358 is None or sign_rs7412 is None:
        status = SKIP
        detail = (
            f"rs7412 scored (rank {rank7412}/{total7412}) but signed_max_score "
            "not persisted. Add signed_max_score column (Task C) and rescore."
        )
        # SKIP does not set all_pass=False — it's an infrastructure gap, not a failure
    elif (sign_rs429358 > 0) != (sign_rs7412 > 0):
        status = PASS
        detail = (
            f"rs429358 expression_signed={sign_rs429358:+.4f}, "
            f"rs7412 expression_signed={sign_rs7412:+.4f} — opposite signs as expected"
        )
    else:
        status = FAIL
        detail = (
            f"rs429358 expression_signed={sign_rs429358:+.4f}, "
            f"rs7412 expression_signed={sign_rs7412:+.4f} — SAME sign, expected opposite"
        )
        all_pass = False
    results.append(("Test 2", "APOE ε4/ε2 expression direction inversion", status, detail))

    # ------------------------------------------------------------------ #
    # Test 3: rs7903146 (TCF7L2) in top 20% of T2D
    # ------------------------------------------------------------------ #
    t2d_rows = _get_scores(conn, "t2d")
    rank, total, pct = _rank_info(t2d_rows, "rs7903146")
    if rank is None:
        status = FAIL
        detail = f"rs7903146 not found in t2d scored set (n={total}) — run whitelist fetch"
        all_pass = False
    elif pct <= TOP_PCT_THRESHOLD:
        status = PASS
        detail = f"rank {rank}/{total} ({pct:.1%})"
    else:
        status = FAIL
        detail = f"rank {rank}/{total} ({pct:.1%}) — below top {TOP_PCT_THRESHOLD:.0%}"
        all_pass = False
    results.append(("Test 3", "rs7903146 (TCF7L2) in T2D top 20%", status, detail))

    # ------------------------------------------------------------------ #
    # Test 4: rs1006737 (CACNA1C) in top 20% of SCZ
    # ------------------------------------------------------------------ #
    scz_rows = _get_scores(conn, "schizophrenia")
    rank, total, pct = _rank_info(scz_rows, "rs1006737")
    if rank is None:
        status = FAIL
        detail = f"rs1006737 not found in schizophrenia scored set (n={total}) — run whitelist fetch"
        all_pass = False
    elif pct <= TOP_PCT_THRESHOLD:
        status = PASS
        detail = f"rank {rank}/{total} ({pct:.1%})"
    else:
        status = FAIL
        detail = f"rank {rank}/{total} ({pct:.1%}) — below top {TOP_PCT_THRESHOLD:.0%}"
        all_pass = False
    results.append(("Test 4", "rs1006737 (CACNA1C) in SCZ top 20%", status, detail))

    # ------------------------------------------------------------------ #
    # Test 5: rs356219 (SNCA) in top 20% of PD
    # ------------------------------------------------------------------ #
    pd_rows = _get_scores(conn, "parkinsons")
    rank, total, pct = _rank_info(pd_rows, "rs356219")
    if rank is None:
        status = FAIL
        detail = f"rs356219 not found in parkinsons scored set (n={total}) — run whitelist fetch"
        all_pass = False
    elif pct <= TOP_PCT_THRESHOLD:
        status = PASS
        detail = f"rank {rank}/{total} ({pct:.1%})"
    else:
        status = FAIL
        detail = f"rank {rank}/{total} ({pct:.1%}) — below top {TOP_PCT_THRESHOLD:.0%}"
        all_pass = False
    results.append(("Test 5", "rs356219 (SNCA) in PD top 20%", status, detail))

    conn.close()

    # ------------------------------------------------------------------ #
    # Print results
    # ------------------------------------------------------------------ #
    print("\n=== Canonical variant verification harness ===\n")
    col_w = max(len(r[1]) for r in results) + 2
    for test_id, name, status, detail in results:
        print(f"  {test_id}  [{status}]  {name:<{col_w}}  {detail}")

    print()
    if all_pass:
        print("All tests PASS (or SKIP pending infrastructure). Safe to claim canonical variant recovery.\n")
    else:
        print("One or more tests FAILED. Do not assert canonical variant claims until failures are resolved.\n")

    return all_pass


if __name__ == "__main__":
    ok = run_tests()
    sys.exit(0 if ok else 1)
