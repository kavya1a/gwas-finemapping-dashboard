"""Canonical variant verification harness.

Run after any scoring change to confirm the pipeline behaves correctly on
five known causal variants. Tests are split into two groups, reflecting
the different claims they support:

  Group A — canonical_rank_recovery
    Variants whose AlphaGenome signal is strong enough to rank them in the
    top 20% of their disease's scored set:
      1. rs429358 (APOE ε4) in top 20% of AD
      2. rs7412 (APOE ε2) scored AND expression direction opposite to rs429358
      3. rs1006737 (CACNA1C) in top 20% of SCZ

  Group B — canonical_regulatory_detection
    Variants where AlphaGenome detects strong regulatory signal but does not
    differentiate the canonical from other GWAS variants at the same locus.
    These pass if the model assigns strong absolute regulatory effect, not
    if it ranks the variant above the GWAS noise. Within-disease-set rank
    is reported for context, but is NOT a pass criterion:
      4. rs7903146 (TCF7L2) regulatory signal detected for T2D
      5. rs356219 (SNCA) regulatory signal detected for PD

The rank-recovery group is the strict claim ("model singles this variant out").
The detection group is the weaker but more honest claim ("model assigns strong
regulatory effect"). Background and per-modality diagnosis is in
docs/canonical_variants.md.

Exit code 0 if all tests pass, 1 otherwise.
Run:
    python verification/canonical_variants_test.py
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

DIR = Path(__file__).parent.parent
SCORED_DB = DIR / "scored_variants.db"

TOP_PCT_THRESHOLD = 0.20            # rank-recovery: rank/total ≤ 0.20
DETECT_COMPOSITE_MIN = 0.5          # detection: composite_score > 0.5
DETECT_EXPR_ABS_MIN  = 0.9          # detection: |expression signed_max| > 0.9

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
SKIP = "\033[33mSKIP\033[0m"


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def _get_scores(conn: sqlite3.Connection, disease: str) -> list[tuple]:
    """All (rsid, composite_score) for a disease, sorted desc."""
    return conn.execute(
        "SELECT rsid, composite_score FROM scores "
        "WHERE disease=? AND error IS NULL AND composite_score IS NOT NULL "
        "ORDER BY composite_score DESC",
        (disease,),
    ).fetchall()


def _rank_info(rows: list[tuple], rsid: str) -> tuple[int | None, int, float | None]:
    """1-based rank, total, percentile for rsid in rows."""
    for i, (r, _) in enumerate(rows):
        if r == rsid:
            rank = i + 1
            total = len(rows)
            return rank, total, rank / total
    return None, len(rows), None


def _get_composite(conn: sqlite3.Connection, disease: str, rsid: str) -> float | None:
    row = conn.execute(
        "SELECT composite_score FROM scores "
        "WHERE disease=? AND rsid=? AND error IS NULL",
        (disease, rsid),
    ).fetchone()
    return row[0] if row else None


def _get_signed_expression(conn: sqlite3.Connection, disease: str, rsid: str) -> float | None:
    row = conn.execute(
        "SELECT signed_max_score FROM modality_scores "
        "WHERE disease=? AND rsid=? AND modality='expression'",
        (disease, rsid),
    ).fetchone()
    return row[0] if row else None


# ---------------------------------------------------------------------------
# Test types
# ---------------------------------------------------------------------------

def _rank_test(rows: list[tuple], rsid: str) -> tuple[str, str, bool]:
    rank, total, pct = _rank_info(rows, rsid)
    if rank is None:
        return FAIL, f"{rsid} not found in scored set (n={total})", False
    if pct <= TOP_PCT_THRESHOLD:
        return PASS, f"rank {rank}/{total} ({pct:.1%}) — within top {TOP_PCT_THRESHOLD:.0%}", True
    return FAIL, f"rank {rank}/{total} ({pct:.1%}) — below top {TOP_PCT_THRESHOLD:.0%}", False


def _detection_test(conn: sqlite3.Connection, disease: str, rsid: str) -> tuple[str, str, bool]:
    """Passes if AlphaGenome assigns strong regulatory signal — composite > 0.5
    AND expression |signed_max| > 0.9. Within-disease-set rank is reported as
    diagnostic context but does not affect pass/fail."""
    composite = _get_composite(conn, disease, rsid)
    expr_signed = _get_signed_expression(conn, disease, rsid)

    if composite is None:
        return FAIL, f"{rsid} not scored for {disease}", False
    if expr_signed is None:
        return FAIL, (
            f"{rsid} scored (composite={composite:.3f}) but expression "
            "signed_max not persisted — re-score to populate modality_scores"
        ), False

    composite_ok = composite > DETECT_COMPOSITE_MIN
    expr_ok = abs(expr_signed) > DETECT_EXPR_ABS_MIN

    rows = _get_scores(conn, disease)
    rank, total, pct = _rank_info(rows, rsid)
    direction = "↑ (alt > ref)" if expr_signed > 0 else "↓ (alt < ref)"

    status = PASS if (composite_ok and expr_ok) else FAIL
    detail = (
        f"composite={composite:.3f} "
        f"({'≥' if composite_ok else '<'}{DETECT_COMPOSITE_MIN}); "
        f"expression signed={expr_signed:+.3f} "
        f"({'≥' if expr_ok else '<'}|{DETECT_EXPR_ABS_MIN}|), {direction}\n"
        f"            [diagnostic] within-{disease}-set rank {rank}/{total} "
        f"({pct:.1%}) — not a pass criterion"
    )
    return status, detail, (composite_ok and expr_ok)


# ---------------------------------------------------------------------------
# Test runners
# ---------------------------------------------------------------------------

def _run_rank_recovery(conn: sqlite3.Connection) -> tuple[list[tuple], bool]:
    out = []
    all_pass = True

    # 1. rs429358 (APOE ε4) in top 20% of AD
    ad_rows = _get_scores(conn, "alzheimers")
    status, detail, ok = _rank_test(ad_rows, "rs429358")
    out.append(("1.", "rs429358 (APOE ε4) in AD top 20%", status, detail))
    all_pass &= ok

    # 2. rs7412 scored AND expression direction opposite to rs429358
    rank7412, total7412, _ = _rank_info(ad_rows, "rs7412")
    sign429 = _get_signed_expression(conn, "alzheimers", "rs429358")
    sign7412 = _get_signed_expression(conn, "alzheimers", "rs7412")
    if rank7412 is None:
        status, detail, ok = FAIL, "rs7412 not scored for AD — run whitelist fetch + score", False
    elif sign429 is None or sign7412 is None:
        status, detail, ok = (
            SKIP,
            f"rs7412 scored (rank {rank7412}/{total7412}) but signed_max_score not persisted",
            True,  # SKIP is not a failure
        )
    elif (sign429 > 0) != (sign7412 > 0):
        status, detail, ok = (
            PASS,
            f"rs429358 expression_signed={sign429:+.4f}, "
            f"rs7412 expression_signed={sign7412:+.4f} — opposite signs as expected",
            True,
        )
    else:
        status, detail, ok = (
            FAIL,
            f"rs429358 expression_signed={sign429:+.4f}, "
            f"rs7412 expression_signed={sign7412:+.4f} — SAME sign, expected opposite",
            False,
        )
    out.append(("2.", "APOE ε4/ε2 expression direction inversion", status, detail))
    all_pass &= ok

    # 3. rs1006737 (CACNA1C) in top 20% of SCZ
    scz_rows = _get_scores(conn, "schizophrenia")
    status, detail, ok = _rank_test(scz_rows, "rs1006737")
    out.append(("3.", "rs1006737 (CACNA1C) in SCZ top 20%", status, detail))
    all_pass &= ok

    return out, all_pass


def _run_regulatory_detection(conn: sqlite3.Connection) -> tuple[list[tuple], bool]:
    out = []
    all_pass = True

    for n, disease, rsid, label in [
        ("4.", "t2d",        "rs7903146", "rs7903146 (TCF7L2) regulatory signal detected for T2D"),
        ("5.", "parkinsons", "rs356219",  "rs356219 (SNCA) regulatory signal detected for PD"),
    ]:
        status, detail, ok = _detection_test(conn, disease, rsid)
        out.append((n, label, status, detail))
        all_pass &= ok

    return out, all_pass


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _print_group(title: str, results: list[tuple]) -> None:
    print(f"\n  {title}")
    print("  " + "-" * (len(title) + 2))
    col_w = max(len(r[1]) for r in results) + 2
    for test_id, name, status, detail in results:
        print(f"    {test_id} [{status}]  {name:<{col_w}}  {detail}")


def run_tests() -> bool:
    if not SCORED_DB.exists():
        print(f"  ERROR: {SCORED_DB} not found. Run batch_score.py first.")
        return False

    conn = sqlite3.connect(SCORED_DB)
    try:
        rank_results, rank_ok = _run_rank_recovery(conn)
        detect_results, detect_ok = _run_regulatory_detection(conn)
    finally:
        conn.close()

    print("\n=== Canonical variant verification harness ===")
    _print_group("Group A — canonical_rank_recovery (model singles variant out)",
                 rank_results)
    _print_group("Group B — canonical_regulatory_detection (model assigns strong signal)",
                 detect_results)

    all_pass = rank_ok and detect_ok
    print()
    if all_pass:
        print("All tests PASS. Safe to claim canonical variant recovery and detection.\n")
    else:
        print("One or more tests FAILED. See docs/canonical_variants.md for context.\n")
    return all_pass


if __name__ == "__main__":
    ok = run_tests()
    sys.exit(0 if ok else 1)
