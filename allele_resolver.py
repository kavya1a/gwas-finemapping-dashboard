"""Resolve ref/alt alleles for rsIDs via Ensembl REST API with SQLite caching.

Rate-limited to 15 req/sec (POST batch, up to 200 IDs at a time).
Multi-allelic variants: picks the alt allele with highest gnomAD global frequency.
When population data is unavailable, falls back to first listed alt allele.
"""

from __future__ import annotations

import json
import sqlite3
import time
import threading
from pathlib import Path
from typing import Optional

import requests

DB_PATH = Path(__file__).parent / "variants.db"
ENSEMBL_POST_URL = "https://rest.ensembl.org/variation/human"
BATCH_SIZE = 200
_MIN_INTERVAL = 1.0 / 15  # 15 req/sec
_last_call_time = 0.0
_rate_lock = threading.Lock()

BLOCKERS_FILE = Path(__file__).parent / "OVERNIGHT_BLOCKERS.md"


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def _init_db(conn: sqlite3.Connection) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS allele_cache (
            rsid           TEXT PRIMARY KEY,
            chrom          TEXT,
            pos            INTEGER,
            ref            TEXT,
            alt            TEXT,
            maf            REAL,
            multi_allelic  INTEGER DEFAULT 0,
            not_found      INTEGER DEFAULT 0,
            fetched_at     INTEGER NOT NULL
        )
    """)
    conn.commit()


# ---------------------------------------------------------------------------
# Rate limiting + HTTP
# ---------------------------------------------------------------------------

def _rate_limit() -> None:
    global _last_call_time
    with _rate_lock:
        now = time.monotonic()
        wait = _MIN_INTERVAL - (now - _last_call_time)
        if wait > 0:
            time.sleep(wait)
        _last_call_time = time.monotonic()


def _post_batch(rsids: list[str], include_pops: bool = True) -> dict:
    """POST up to BATCH_SIZE rsIDs to Ensembl.

    include_pops=True adds ?pops=1 for allele frequency data (needed for
    multi-allelic alt selection). False omits it — lighter response, faster
    for large sets where MAF is not required.
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = json.dumps({"ids": rsids})
    url = ENSEMBL_POST_URL + ("?pops=1" if include_pops else "")

    for attempt in range(6):
        _rate_limit()
        try:
            resp = requests.post(url, data=payload, headers=headers, timeout=90)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 429:
                retry_after = float(resp.headers.get("Retry-After", 2 ** (attempt + 1)))
                print(f"  429 — sleeping {retry_after:.1f}s (attempt {attempt + 1})")
                time.sleep(retry_after)
            elif resp.status_code in (400, 404):
                return {}
            else:
                resp.raise_for_status()
        except requests.RequestException as exc:
            if attempt >= 5:
                raise
            time.sleep(2 ** attempt)
    return {}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _pick_best_alt(
    ref: str,
    alt_alleles: list[str],
    populations: list[dict],
) -> tuple[str, float | None]:
    """Pick the alt allele with highest global frequency; fall back to first alt."""
    if not alt_alleles:
        return "", None

    if len(alt_alleles) == 1:
        # Single alt — still try to find its frequency
        maf = None
        for entry in populations:
            if entry.get("allele") == alt_alleles[0]:
                f = entry.get("frequency", 0.0)
                maf = max(maf or 0.0, f)
        return alt_alleles[0], maf

    # Multi-allelic: build allele → max frequency across populations
    freq_map: dict[str, float] = {}
    for entry in populations:
        a = entry.get("allele", "")
        f = entry.get("frequency", 0.0)
        if a and a != ref and a in alt_alleles:
            freq_map[a] = max(freq_map.get(a, 0.0), f)

    if freq_map:
        best = max(freq_map, key=lambda a: freq_map[a])
        return best, freq_map[best]
    # No population data — take first alt, return None MAF
    return alt_alleles[0], None


def _parse_variant(rsid: str, data: dict) -> Optional[dict]:
    """Parse one Ensembl variant dict; returns None if schema is unparseable."""
    if not isinstance(data, dict):
        return None

    mappings = data.get("mappings", [])
    if not isinstance(mappings, list):
        return None

    grch38 = [m for m in mappings if isinstance(m, dict) and m.get("assembly_name") == "GRCh38"]
    mapping = grch38[0] if grch38 else (mappings[0] if mappings else None)
    if not mapping:
        return None

    chrom_raw = str(mapping.get("seq_region_name", ""))
    if not chrom_raw:
        return None
    chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"

    pos = mapping.get("start")
    if pos is None:
        return None

    allele_string = mapping.get("allele_string", "")
    if not allele_string or "/" not in allele_string:
        return None

    parts = allele_string.split("/")
    ref = parts[0]
    alt_candidates = [a for a in parts[1:] if a and a != ref]
    if not alt_candidates:
        return None

    populations = data.get("populations", [])
    multi_allelic = len(alt_candidates) > 1
    alt, maf = _pick_best_alt(ref, alt_candidates, populations)

    return {
        "rsid": rsid,
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "maf": maf,
        "multi_allelic": int(multi_allelic),
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _write_blocker(message: str) -> None:
    with open(BLOCKERS_FILE, "a") as f:
        f.write(f"\n## BLOCKER [{time.strftime('%Y-%m-%d %H:%M:%S')}]\n{message}\n")
    print(f"!! BLOCKER written to {BLOCKERS_FILE}")


def resolve_alleles_batch(
    rsids: list[str],
    db_path: Path = DB_PATH,
    include_pops: bool = True,
) -> dict[str, Optional[dict]]:
    """
    Resolve ref/alt alleles for a list of rsIDs (cache-first).

    Returns dict: rsid → {rsid, chrom, pos, ref, alt, maf, multi_allelic}
    or None if unresolvable.
    Raises RuntimeError (and writes OVERNIGHT_BLOCKERS.md) on unexpected schema.
    """
    results: dict[str, Optional[dict]] = {}
    to_fetch: list[str] = []

    conn = sqlite3.connect(db_path)
    _init_db(conn)
    try:
        for rsid in rsids:
            row = conn.execute(
                "SELECT chrom, pos, ref, alt, maf, multi_allelic, not_found "
                "FROM allele_cache WHERE rsid = ?",
                (rsid,),
            ).fetchone()
            if row is not None:
                chrom, pos, ref, alt, maf, multi_allelic, not_found = row
                results[rsid] = (
                    None
                    if not_found
                    else {
                        "rsid": rsid, "chrom": chrom, "pos": pos,
                        "ref": ref, "alt": alt, "maf": maf,
                        "multi_allelic": bool(multi_allelic),
                    }
                )
            else:
                to_fetch.append(rsid)
    finally:
        conn.close()

    if not to_fetch:
        return results

    print(f"  Fetching {len(to_fetch)} rsIDs from Ensembl ({BATCH_SIZE}/batch)...")
    n_batches = (len(to_fetch) + BATCH_SIZE - 1) // BATCH_SIZE

    for batch_idx in range(n_batches):
        batch = to_fetch[batch_idx * BATCH_SIZE : (batch_idx + 1) * BATCH_SIZE]
        print(f"  Batch {batch_idx + 1}/{n_batches}: {len(batch)} rsIDs")

        raw = _post_batch(batch, include_pops=include_pops)

        if not isinstance(raw, dict):
            msg = (
                f"Ensembl POST returned unexpected type {type(raw).__name__} "
                f"(expected dict). First rsID in batch: {batch[0]}"
            )
            _write_blocker(msg)
            raise RuntimeError(msg)

        conn = sqlite3.connect(db_path)
        _init_db(conn)
        try:
            for rsid in batch:
                variant_data = raw.get(rsid)
                if variant_data is None or not isinstance(variant_data, dict):
                    conn.execute(
                        "INSERT OR REPLACE INTO allele_cache "
                        "(rsid,chrom,pos,ref,alt,maf,multi_allelic,not_found,fetched_at) "
                        "VALUES (?,NULL,NULL,NULL,NULL,NULL,0,1,?)",
                        (rsid, int(time.time())),
                    )
                    results[rsid] = None
                    continue

                parsed = _parse_variant(rsid, variant_data)
                if parsed is None:
                    # Log but don't raise — some rsIDs have sparse Ensembl records
                    keys_seen = list(variant_data.keys())[:8]
                    print(
                        f"  WARNING: could not parse {rsid} "
                        f"(keys={keys_seen}, mappings={len(variant_data.get('mappings', []))})"
                    )
                    conn.execute(
                        "INSERT OR REPLACE INTO allele_cache "
                        "(rsid,chrom,pos,ref,alt,maf,multi_allelic,not_found,fetched_at) "
                        "VALUES (?,NULL,NULL,NULL,NULL,NULL,0,1,?)",
                        (rsid, int(time.time())),
                    )
                    results[rsid] = None
                else:
                    conn.execute(
                        "INSERT OR REPLACE INTO allele_cache "
                        "(rsid,chrom,pos,ref,alt,maf,multi_allelic,not_found,fetched_at) "
                        "VALUES (?,?,?,?,?,?,?,0,?)",
                        (
                            parsed["rsid"], parsed["chrom"], parsed["pos"],
                            parsed["ref"], parsed["alt"], parsed["maf"],
                            parsed["multi_allelic"], int(time.time()),
                        ),
                    )
                    results[rsid] = parsed
            conn.commit()
        finally:
            conn.close()

    return results
