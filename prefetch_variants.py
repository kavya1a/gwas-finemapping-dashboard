"""Pre-fetch top GWAS variants per disease and resolve alleles.

Saves to preloaded_variants.db keyed by (disease, rsid).
Marks scoring_skipped_reason at fetch time:
  - allele_resolution_failed: Ensembl could not resolve ref/alt
  - indel_not_supported: ref or alt is '-', empty, or non-ACGTN
Safe to re-run — already-fetched variants are skipped.
"""

from __future__ import annotations

import re
import sqlite3
import time
from pathlib import Path

import yaml

from gwas_catalog import fetch_gwas_variants

DIR = Path(__file__).parent
DB_PATH = DIR / "preloaded_variants.db"
BLOCKERS_FILE = DIR / "OVERNIGHT_BLOCKERS.md"

_cfg = yaml.safe_load((DIR / "config.yaml").read_text()) if (DIR / "config.yaml").exists() else {}
VARIANTS_PER_DISEASE: int = _cfg.get("pipeline", {}).get("variants_per_disease", 300)

DISEASES = list({
    "alzheimers": None,
    "t2d": None,
    "schizophrenia": None,
    "parkinsons": None,
}.keys())

_VALID_BASES = re.compile(r'^[ACGTNacgtn]+$')


def _skip_reason(ref, alt) -> str | None:
    if ref is None or alt is None:
        return "allele_resolution_failed"
    if not ref or ref == '-' or not _VALID_BASES.match(ref):
        return "indel_not_supported"
    if not alt or alt == '-' or not _VALID_BASES.match(alt):
        return "indel_not_supported"
    return None


def _init_db(conn: sqlite3.Connection) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            disease               TEXT NOT NULL,
            rsid                  TEXT NOT NULL,
            chrom                 TEXT,
            pos                   INTEGER,
            ref                   TEXT,
            alt                   TEXT,
            maf                   REAL,
            multi_allelic         INTEGER DEFAULT 0,
            p_value               REAL,
            fetched_at            INTEGER DEFAULT 0,
            scoring_skipped_reason TEXT,
            whitelist             INTEGER DEFAULT 0,
            PRIMARY KEY (disease, rsid)
        )
    """)
    # Migrate existing DBs
    cols = [r[1] for r in conn.execute("PRAGMA table_info(variants)").fetchall()]
    if "scoring_skipped_reason" not in cols:
        conn.execute("ALTER TABLE variants ADD COLUMN scoring_skipped_reason TEXT")
    if "whitelist" not in cols:
        conn.execute("ALTER TABLE variants ADD COLUMN whitelist INTEGER DEFAULT 0")
    conn.commit()


def _already_fetched(conn: sqlite3.Connection, disease: str) -> int:
    return conn.execute(
        "SELECT COUNT(*) FROM variants WHERE disease = ?", (disease,)
    ).fetchone()[0]


def _write_blocker(msg: str) -> None:
    with open(BLOCKERS_FILE, "a") as f:
        f.write(f"\n## BLOCKER [{time.strftime('%Y-%m-%d %H:%M:%S')}] prefetch_variants\n{msg}\n")
    print(f"!! BLOCKER: {msg}")


def prefetch_all() -> None:
    conn = sqlite3.connect(DB_PATH)
    _init_db(conn)
    conn.close()

    for slug in DISEASES:
        conn = sqlite3.connect(DB_PATH)
        already = _already_fetched(conn, slug)
        conn.close()

        if already >= VARIANTS_PER_DISEASE:
            print(f"[{slug}] already have {already} variants (>= {VARIANTS_PER_DISEASE}) — skipping fetch")
            continue

        print(f"\n[{slug}] Fetching up to {VARIANTS_PER_DISEASE} variants...")
        try:
            variants = fetch_gwas_variants(slug, max_variants=VARIANTS_PER_DISEASE)
        except Exception as exc:
            _write_blocker(f"GWAS Catalog fetch failed for {slug}: {exc}")
            print(f"  ERROR — see {BLOCKERS_FILE}. Continuing with next disease.")
            continue

        if not variants:
            _write_blocker(f"GWAS Catalog returned 0 variants for {slug}")
            continue

        n_scoreable = sum(1 for v in variants if _skip_reason(v.get("ref"), v.get("alt")) is None)
        n_indel = sum(1 for v in variants if _skip_reason(v.get("ref"), v.get("alt")) == "indel_not_supported")
        n_afail = sum(1 for v in variants if _skip_reason(v.get("ref"), v.get("alt")) == "allele_resolution_failed")
        print(f"  Got {len(variants)}: {n_scoreable} scoreable, {n_indel} indel_skip, {n_afail} allele_fail")

        if n_scoreable == 0:
            _write_blocker(
                f"0 scoreable variants for {slug} — all filtered at allele stage. "
                "Check allele_resolver.py and variants.db."
            )

        conn = sqlite3.connect(DB_PATH)
        _init_db(conn)
        try:
            for v in variants:
                reason = _skip_reason(v.get("ref"), v.get("alt"))
                conn.execute(
                    """INSERT OR REPLACE INTO variants
                       (disease,rsid,chrom,pos,ref,alt,maf,multi_allelic,
                        p_value,fetched_at,scoring_skipped_reason,whitelist)
                       VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (
                        slug, v["rsid"], v.get("chrom"), v.get("pos"),
                        v.get("ref"), v.get("alt"), v.get("maf"),
                        int(v.get("multi_allelic", False)),
                        v.get("p_value"), int(time.time()), reason,
                        int(bool(v.get("whitelist", False))),
                    ),
                )
            conn.commit()
            print(f"  Saved {len(variants)} variants to {DB_PATH}")
        finally:
            conn.close()

    # Summary
    conn = sqlite3.connect(DB_PATH)
    total = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    scoreable = conn.execute(
        "SELECT COUNT(*) FROM variants WHERE scoring_skipped_reason IS NULL"
    ).fetchone()[0]
    conn.close()
    print(f"\nDone. {total} total variants ({scoreable} scoreable) in preloaded_variants.db.")


if __name__ == "__main__":
    prefetch_all()
