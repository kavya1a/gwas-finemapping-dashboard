"""
Local CADD lookup via remote tabix (CADD GRCh38 v1.7, whole-genome SNV file).

Replaces the REST API to avoid rate limits. Uses tabix HTTP byte-range requests
against the CADD server (Accept-Ranges: bytes), so the 87 GB file is never fully
downloaded — only the relevant genomic windows are fetched on demand. Results are
cached in cadd_cache.db (same schema used by tewhey_analysis.py).

Requirements: htslib tabix on PATH (install via: brew install htslib).

Usage:
    from cadd_local import lookup_cadd_single, populate_cadd_from_parquet
    phred, raw = lookup_cadd_single("chr19", 44908684, "T", "C")

    # After tewhey_analysis.py finishes:
    populate_cadd_from_parquet()  # fills cadd_phred column in tewhey_mpra.parquet
"""

from __future__ import annotations

import shutil
import sqlite3
import subprocess
import time
from pathlib import Path

import pandas as pd

DIR = Path(__file__).parent
CADD_CACHE_DB = DIR / "cadd_cache.db"
PARQUET_PATH = DIR / "tewhey_mpra.parquet"

CADD_TABIX_URL = (
    "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
    "whole_genome_SNVs.tsv.gz"
)

# CADD tabix index uses bare chromosome numbers (no "chr" prefix).
def _bare(chrom: str) -> str:
    return chrom.lstrip("chr")


# ---------------------------------------------------------------------------
# Cache init
# ---------------------------------------------------------------------------

def _init_cache() -> None:
    conn = sqlite3.connect(CADD_CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cadd (
            chrom       TEXT,
            pos         INTEGER,
            ref         TEXT,
            alt         TEXT,
            phred       REAL,
            raw_score   REAL,
            not_found   INTEGER DEFAULT 0,
            fetched_at  INTEGER,
            PRIMARY KEY (chrom, pos, ref, alt)
        )
    """)
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Core tabix query
# ---------------------------------------------------------------------------

def _tabix_region(chrom: str, pos: int) -> list[dict]:
    """Query all CADD alleles at a single position via remote tabix."""
    if not shutil.which("tabix"):
        raise RuntimeError("tabix not found on PATH — install htslib: brew install htslib")
    region = f"{_bare(chrom)}:{pos}-{pos}"
    try:
        proc = subprocess.run(
            ["tabix", CADD_TABIX_URL, region],
            capture_output=True, text=True, timeout=60,
        )
        if proc.returncode != 0 or not proc.stdout.strip():
            return []
        rows = []
        for line in proc.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            rows.append({
                "chrom": parts[0], "pos": int(parts[1]),
                "ref": parts[2], "alt": parts[3],
                "raw_score": float(parts[4]), "phred": float(parts[5]),
            })
        return rows
    except subprocess.TimeoutExpired:
        return []
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def lookup_cadd_single(
    chrom: str, pos: int, ref: str, alt: str
) -> tuple[float | None, float | None]:
    """
    Return (PHRED, RawScore) for one variant; (None, None) if absent from CADD.

    Checks cadd_cache.db first. On miss, queries the remote file via tabix
    and caches the result (including not_found entries) for future calls.
    """
    _init_cache()

    # Cache check
    conn = sqlite3.connect(CADD_CACHE_DB)
    row = conn.execute(
        "SELECT phred, raw_score, not_found FROM cadd "
        "WHERE chrom=? AND pos=? AND ref=? AND alt=?",
        (chrom, pos, ref, alt),
    ).fetchone()
    conn.close()
    if row is not None:
        return (None, None) if row[2] else (row[0], row[1])

    # Remote lookup
    alleles = _tabix_region(chrom, pos)
    phred = raw = None
    for entry in alleles:
        if entry["ref"].upper() == ref.upper() and entry["alt"].upper() == alt.upper():
            phred, raw = entry["phred"], entry["raw_score"]
            break

    # Cache result
    conn = sqlite3.connect(CADD_CACHE_DB)
    conn.execute(
        "INSERT OR REPLACE INTO cadd "
        "(chrom,pos,ref,alt,phred,raw_score,not_found,fetched_at) "
        "VALUES (?,?,?,?,?,?,?,?)",
        (chrom, pos, ref, alt, phred, raw, int(phred is None), int(time.time())),
    )
    conn.commit()
    conn.close()
    return phred, raw


def populate_cadd_from_parquet(
    parquet_path: Path = PARQUET_PATH,
    force: bool = False,
) -> None:
    """
    Fill cadd_phred column in tewhey_mpra.parquet using local tabix lookup.

    Run this after tewhey_analysis.py completes. On re-run, skips variants
    already in cadd_cache.db unless force=True.

    Args:
        parquet_path: Path to the Tewhey parquet file.
        force: If True, re-fetch all variants even if cached.
    """
    if not parquet_path.exists():
        raise FileNotFoundError(f"{parquet_path} not found — run tewhey_analysis.py first")

    df = pd.read_parquet(parquet_path)
    scoreable = df.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
    print(f"Looking up CADD scores for {len(scoreable)} variants via remote tabix …")

    if force:
        conn = sqlite3.connect(CADD_CACHE_DB)
        conn.execute("DELETE FROM cadd")
        conn.commit()
        conn.close()
        print("  Cache cleared (force=True).")

    _init_cache()
    conn = sqlite3.connect(CADD_CACHE_DB)
    cached_keys = {
        (r[0], r[1], r[2], r[3])
        for r in conn.execute("SELECT chrom, pos, ref, alt FROM cadd").fetchall()
    }
    conn.close()

    to_fetch = [
        row for _, row in scoreable.iterrows()
        if (row.chrom, int(row.pos), row.ref, row.alt) not in cached_keys
    ]
    print(f"  {len(cached_keys)} cached, {len(to_fetch)} to fetch via tabix")

    for i, row in enumerate(to_fetch):
        phred, raw = lookup_cadd_single(row.chrom, int(row.pos), row.ref, row.alt)
        if (i + 1) % 200 == 0 or i == len(to_fetch) - 1:
            print(f"  [{i+1}/{len(to_fetch)}] {row.rsid if hasattr(row, 'rsid') else ''} "
                  f"phred={phred}")

    # Re-read cache and merge into parquet
    conn = sqlite3.connect(CADD_CACHE_DB)
    cadd_rows = conn.execute(
        "SELECT chrom, pos, ref, alt, phred FROM cadd WHERE not_found=0"
    ).fetchall()
    conn.close()

    cadd_df = pd.DataFrame(cadd_rows, columns=["chrom", "pos", "ref", "alt", "cadd_phred"])
    cadd_df["pos"] = cadd_df["pos"].astype(int)

    df["pos"] = df["pos"].astype("Int64")
    cadd_df["pos"] = cadd_df["pos"].astype("Int64")

    df = df.drop(columns=["cadd_phred"], errors="ignore")
    df = df.merge(cadd_df[["chrom", "pos", "ref", "alt", "cadd_phred"]],
                  on=["chrom", "pos", "ref", "alt"], how="left")

    df.to_parquet(parquet_path, index=False)
    n_found = df["cadd_phred"].notna().sum()
    print(f"\nCADD scores added: {n_found}/{len(df)} variants have phred scores.")
    print(f"Updated parquet saved to {parquet_path}")


if __name__ == "__main__":
    # Quick sanity check: look up APOE rs429358 (chr19:44908684 T>C)
    phred, raw = lookup_cadd_single("chr19", 44908684, "T", "C")
    print(f"rs429358 (APOE ε4): PHRED={phred}, RawScore={raw}")
    if phred is not None:
        print("Remote tabix lookup working correctly.")
    else:
        print("Variant not found in CADD (or tabix failed).")
