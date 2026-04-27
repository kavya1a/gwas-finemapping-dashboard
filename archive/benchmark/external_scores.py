"""Fetch CADD scores and GWAS p-values for benchmark variants.

CADD REST API: https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chr}:{pos}
  Returns list of dicts with keys: Chrom, Pos, Ref, Alt, RawScore, PHRED.
  Rate limit: experimental API, not for mass use — we add 0.3s delay between calls.
  For missing variants (404), PHRED defaults to 0.0 (least deletedious).

GWAS Catalog v1 API: /rest/api/singleNucleotidePolymorphisms/{rsid}/associations
  Returns associations for the rsid. We take the smallest p-value found.
  If rsid not in GWAS Catalog, p-value defaults to 1.0 (not associated).

Both lookups are cached in benchmark/cache/external_scores.json.
"""

from __future__ import annotations

import json
import time
import urllib.request
import urllib.error
from pathlib import Path

CACHE_FILE = Path(__file__).parent / "cache" / "external_scores.json"
CADD_URL = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chrom}:{pos}"
GWAS_URL = (
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms"
    "/{rsid}/associations?projection=associationBySnp&page=0&size=10"
)

CADD_DELAY_S = 0.35   # polite delay between CADD calls
GWAS_DELAY_S = 0.12   # polite delay between GWAS catalog calls


def _fetch_json(url: str, timeout: int = 20) -> dict | list | None:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        if e.code in (404, 400):
            return None
        raise
    except Exception:
        return None


def get_cadd_phred(chrom: str, pos: int, ref: str, alt: str) -> float | None:
    """Return CADD PHRED score for a variant, or None if unavailable."""
    raw_chrom = chrom.replace("chr", "")
    url = CADD_URL.format(chrom=raw_chrom, pos=pos)
    data = _fetch_json(url)
    if not data:
        return None
    # data is a list of dicts; match by ref and alt
    for row in data:
        if row.get("Ref") == ref and row.get("Alt") == alt:
            try:
                return float(row["PHRED"])
            except (KeyError, ValueError):
                pass
    # If exact ref/alt match not found, return max PHRED at position
    phreds = []
    for row in data:
        try:
            phreds.append(float(row["PHRED"]))
        except (KeyError, ValueError):
            pass
    return max(phreds) if phreds else None


def get_gwas_pvalue(rsid: str) -> float | None:
    """Return minimum GWAS p-value for rsid across all GWAS Catalog associations."""
    if not rsid.startswith("rs"):
        return None
    url = GWAS_URL.format(rsid=rsid)
    data = _fetch_json(url)
    if not data:
        return None
    assocs = data.get("_embedded", {}).get("associations", [])
    pvals = []
    for a in assocs:
        try:
            mantissa = float(a.get("pvalueMantissa", 1))
            exponent = int(a.get("pvalueExponent", 0))
            pvals.append(mantissa * 10 ** exponent)
        except (TypeError, ValueError):
            pass
    return min(pvals) if pvals else None


def fetch_all_external_scores(
    variants: list[dict],
    use_cache: bool = True,
) -> dict[str, dict]:
    """Fetch CADD + GWAS p-values for all variants; returns dict keyed by rsid.

    Each value: {"cadd_phred": float|None, "gwas_pvalue": float|None}.
    """
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    cache: dict = {}
    if use_cache and CACHE_FILE.exists():
        with open(CACHE_FILE) as f:
            cache = json.load(f)

    changed = False
    for i, v in enumerate(variants):
        key = v["rsid"]
        if key in cache:
            continue

        print(f"  [{i+1}/{len(variants)}] {key} — CADD...", end="", flush=True)
        cadd = get_cadd_phred(v["chrom"], v["pos"], v["ref"], v["alt"])
        time.sleep(CADD_DELAY_S)

        print(f" {cadd or 'N/A'}  GWAS...", end="", flush=True)
        gwas = get_gwas_pvalue(key)
        time.sleep(GWAS_DELAY_S)

        print(f" {gwas or 'N/A'}")
        cache[key] = {"cadd_phred": cadd, "gwas_pvalue": gwas}
        changed = True

    if changed:
        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f)

    return cache
