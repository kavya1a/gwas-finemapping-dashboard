"""Generate pipeline yield table across all 4 diseases.

Reads preloaded_variants.db (fetch + allele + skip stages) and
scored_variants.db (scoring outcome). Prints a markdown table.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

PRELOADED_DB = Path(__file__).parent / "preloaded_variants.db"
SCORED_DB = Path(__file__).parent / "scored_variants.db"

DISEASES = ["alzheimers", "t2d", "schizophrenia", "parkinsons"]
DISEASE_LABELS = {
    "alzheimers":    "Alzheimer's",
    "t2d":           "T2D",
    "schizophrenia": "Schizophrenia",
    "parkinsons":    "Parkinson's",
}


def main() -> None:
    pconn = sqlite3.connect(PRELOADED_DB) if PRELOADED_DB.exists() else None
    sconn = sqlite3.connect(SCORED_DB) if SCORED_DB.exists() else None

    rows = []
    for disease in DISEASES:
        label = DISEASE_LABELS[disease]

        # --- Prefetch stage ---
        fetched = 0
        indel_skip = 0
        allele_fail = 0
        if pconn:
            fetched = pconn.execute(
                "SELECT COUNT(*) FROM variants WHERE disease=?", (disease,)
            ).fetchone()[0]
            indel_skip = pconn.execute(
                "SELECT COUNT(*) FROM variants WHERE disease=? AND scoring_skipped_reason='indel_not_supported'",
                (disease,)
            ).fetchone()[0]
            allele_fail = pconn.execute(
                "SELECT COUNT(*) FROM variants WHERE disease=? AND scoring_skipped_reason='allele_resolution_failed'",
                (disease,)
            ).fetchone()[0]

        # --- Flanking skip (new reason added in batch_score.py) ---
        flanking_skip = 0
        if pconn:
            flanking_skip = pconn.execute(
                "SELECT COUNT(*) FROM variants "
                "WHERE disease=? AND scoring_skipped_reason='insufficient_flanking_sequence'",
                (disease,)
            ).fetchone()[0]

        # --- Scoring stage ---
        # Collect api_error rsIDs from both DBs into a set to deduplicate.
        # The same event is marked in preloaded_variants.db (scoring_skipped_reason)
        # AND in scored_variants.db (error column) — union prevents double-counting.
        error_rsids: set[str] = set()
        if pconn:
            for (rsid,) in pconn.execute(
                "SELECT rsid FROM variants WHERE disease=? AND scoring_skipped_reason='api_error'",
                (disease,),
            ).fetchall():
                error_rsids.add(rsid)
        scored_ok = 0
        if sconn:
            for (rsid,) in sconn.execute(
                "SELECT rsid FROM scores WHERE disease=? AND (error IS NOT NULL OR composite_score IS NULL)",
                (disease,),
            ).fetchall():
                error_rsids.add(rsid)
            scored_ok = sconn.execute(
                "SELECT COUNT(*) FROM scores "
                "WHERE disease=? AND error IS NULL AND composite_score IS NOT NULL",
                (disease,),
            ).fetchone()[0]

        rows.append({
            "disease": label,
            "fetched": fetched,
            "indel_skip": indel_skip,
            "allele_fail": allele_fail,
            "flanking_skip": flanking_skip,
            "api_err": len(error_rsids),
            "scored_ok": scored_ok,
        })

    if pconn:
        pconn.close()
    if sconn:
        sconn.close()

    # --- Print table ---
    print("\n## Pipeline Yield Table\n")
    header = (
        f"{'Disease':<15} {'Fetched':>8} {'Indel skip':>11} {'Allele fail':>12} "
        f"{'Flanking':>9} {'API error':>10} {'Scored OK':>10} {'Yield':>7}"
    )
    print(header)
    print("-" * len(header))
    total_fetched = total_indel = total_afail = total_flanking = total_apierr = total_ok = 0
    for r in rows:
        fetched     = r["fetched"]
        indel       = r["indel_skip"]
        afail       = r["allele_fail"]
        flanking    = r["flanking_skip"]
        apierr      = r["api_err"]
        ok          = r["scored_ok"]
        yield_pct   = f"{100*ok/fetched:.0f}%" if fetched else "—"
        print(
            f"{r['disease']:<15} {fetched:>8} {indel:>11} {afail:>12} "
            f"{flanking:>9} {apierr:>10} {ok:>10} {yield_pct:>7}"
        )
        total_fetched += fetched;  total_indel   += indel
        total_afail   += afail;    total_flanking += flanking
        total_apierr  += apierr;   total_ok       += ok
    print("-" * len(header))
    total_yield = f"{100*total_ok/total_fetched:.0f}%" if total_fetched else "—"
    print(
        f"{'TOTAL':<15} {total_fetched:>8} {total_indel:>11} {total_afail:>12} "
        f"{total_flanking:>9} {total_apierr:>10} {total_ok:>10} {total_yield:>7}"
    )
    print()


if __name__ == "__main__":
    main()
