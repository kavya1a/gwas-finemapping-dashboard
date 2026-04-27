"""Phase 3: Generalize saturation finding to blood trait GWAS variants.

Fetches platelet count (GWAS Catalog EFO_0004615) and hemoglobin
(EFO_0004611) GWAS variants, scores them via AlphaGenome, and shows
the saturation CDF replicates on a second independent regulatory-enriched
dataset distinct from Tewhey MPRA.

Gate 3 threshold: ≥75% of blood trait variants have |expression_subscore| > 0.9.
If met → saturation is a general property of regulatory-enriched selection.
If not → Tewhey-specific, paper framed accordingly.

Run:
    python phase3_blood_traits.py
"""

from __future__ import annotations

import concurrent.futures
import os
import re
import sqlite3
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import requests
import yaml
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

DIR = Path(__file__).parent
FIG_DIR = DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

CACHE_DB = DIR / "phase3_blood_cache.db"
CONFIG_PATH = DIR / "config.yaml"

VARIANT_TIMEOUT_SECS = 60
TARGET_N = 200          # variants per trait to fetch
GATE3_THRESHOLD = 0.75  # fraction with |score| > 0.9 to pass

# Colors consistent with saturation_figure.py
C_TEWHEY = "#2166ac"
C_GWAS   = "#d6604d"
C_BLOOD  = "#762a83"   # purple for blood traits
C_UNIF   = "#999999"
LW = 1.8
ALPHA_FILL = 0.08

# Blood traits to fetch
BLOOD_TRAITS = {
    "platelet_count": "EFO_0004615",
    "hemoglobin":     "EFO_0004611",
}

VALID_BASES = re.compile(r"^[ACGTNacgtn]+$")


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

def _init_cache() -> None:
    conn = sqlite3.connect(CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scores (
            rsid                TEXT NOT NULL,
            trait               TEXT NOT NULL,
            chrom               TEXT,
            pos                 INTEGER,
            ref                 TEXT,
            alt                 TEXT,
            expression_subscore REAL,
            error               TEXT,
            scored_at           INTEGER,
            PRIMARY KEY (rsid, trait)
        )
    """)
    conn.commit()
    conn.close()


def _load_cache(trait: str) -> dict[str, dict]:
    if not CACHE_DB.exists():
        return {}
    conn = sqlite3.connect(CACHE_DB)
    rows = conn.execute(
        "SELECT rsid, expression_subscore, error FROM scores WHERE trait = ?",
        (trait,),
    ).fetchall()
    conn.close()
    return {r[0]: {"expression_subscore": r[1], "error": r[2]} for r in rows}


def _save_score(rsid: str, trait: str, variant: dict, score: float | None, error: str | None) -> None:
    conn = sqlite3.connect(CACHE_DB)
    conn.execute(
        "INSERT OR REPLACE INTO scores "
        "(rsid,trait,chrom,pos,ref,alt,expression_subscore,error,scored_at) "
        "VALUES (?,?,?,?,?,?,?,?,?)",
        (rsid, trait,
         variant.get("chrom"), variant.get("pos"),
         variant.get("ref"), variant.get("alt"),
         score, error, int(time.time())),
    )
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# GWAS Catalog fetch (EFO IDs)
# ---------------------------------------------------------------------------

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
_RATE_LIMIT = 1 / 15


def _fetch_trait_variants(trait: str, efo_id: str, max_n: int = TARGET_N) -> list[dict]:
    """Fetch top-p-value variants for a trait from GWAS Catalog v2."""
    variants = []
    seen_rsids: set[str] = set()
    params = {"efo_id": efo_id, "page": 0, "size": 50}

    while len(variants) < max_n:
        time.sleep(_RATE_LIMIT)
        try:
            resp = requests.get(GWAS_API, params=params, timeout=30)
            resp.raise_for_status()
        except Exception as e:
            print(f"  GWAS API error page {params['page']}: {e}")
            break

        data = resp.json()
        items = data.get("_embedded", {}).get("associations", [])
        if not items:
            break

        for item in items:
            snp_alleles = item.get("snp_allele") or []
            if not snp_alleles:
                continue
            rsid = snp_alleles[0].get("rs_id", "")
            if not rsid or not rsid.startswith("rs"):
                continue
            if rsid in seen_rsids:
                continue

            locations = item.get("locations") or []
            if not locations:
                continue
            loc_parts = str(locations[0]).split(":")
            if len(loc_parts) < 2:
                continue
            try:
                pos = int(loc_parts[1])
            except ValueError:
                continue

            mantissa = item.get("pvalue_mantissa")
            exponent = item.get("pvalue_exponent")
            try:
                p_value = float(mantissa) * 10 ** float(exponent) if mantissa is not None else None
            except (TypeError, ValueError):
                p_value = None

            seen_rsids.add(rsid)
            variants.append({
                "rsid": rsid,
                "chrom": loc_parts[0],
                "pos": pos,
                "ref": None,
                "alt": None,
                "p_value": p_value,
            })

        if not data.get("_links", {}).get("next"):
            break
        params["page"] += 1

    return variants[:max_n]


# ---------------------------------------------------------------------------
# Allele resolution via Ensembl (reuse allele_resolver logic)
# ---------------------------------------------------------------------------

def _resolve_alleles(variants: list[dict]) -> list[dict]:
    from allele_resolver import resolve_alleles_batch, DB_PATH
    rsids = [v["rsid"] for v in variants]
    resolved = resolve_alleles_batch(rsids, db_path=DB_PATH, include_pops=False)
    for v in variants:
        info = resolved.get(v["rsid"])
        if info:
            v["ref"] = info["ref"]
            v["alt"] = info["alt"]
            v["chrom"] = info["chrom"]
            v["pos"] = info["pos"]
    return variants


# ---------------------------------------------------------------------------
# Scoring (expression subscore only)
# ---------------------------------------------------------------------------

def _load_k562_profile():
    from scoring.tissue_config import TissueProfile
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    k562 = cfg.get("tissue_profiles", {}).get("tewhey_k562", {})
    return TissueProfile(
        display_name="K562 / Blood-lineage",
        biosample_keywords=k562.get("biosample_keywords", []),
        gtex_keywords=k562.get("gtex_keywords", []),
    )


def _get_expression_subscore(result: dict) -> float | None:
    """Extract expression_subscore from score_single_variant result."""
    mb = result.get("modality_breakdown")
    if mb is None or mb.empty:
        return None
    expr = mb[mb["modality"] == "expression"]
    if expr.empty:
        return None
    row = expr.iloc[0]
    return float(row.get("signed_max_score")) if row.get("signed_max_score") is not None else None


def _score_one(model, variant_input, profile, rsid: str) -> dict:
    from scoring.composite import score_single_variant
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    future = executor.submit(score_single_variant, model, variant_input, profile)
    try:
        result = future.result(timeout=VARIANT_TIMEOUT_SECS)
        if result.get("error"):
            return {"error": result["error"], "expression_subscore": None}
        return {"expression_subscore": _get_expression_subscore(result), "error": None}
    except concurrent.futures.TimeoutError:
        return {"error": "api_timeout", "expression_subscore": None}
    except Exception as exc:
        return {"error": str(exc), "expression_subscore": None}
    finally:
        executor.shutdown(wait=False)


def _score_trait(trait: str, variants: list[dict], model, profile) -> list[float]:
    from scoring.composite import VariantInput

    cache = _load_cache(trait)
    scoreable = [
        v for v in variants
        if v.get("ref") and v.get("alt")
        and VALID_BASES.match(str(v["ref"])) and VALID_BASES.match(str(v["alt"]))
    ]

    n_total = len(scoreable)
    n_cached = sum(1 for v in scoreable if v["rsid"] in cache and cache[v["rsid"]]["error"] is None)
    print(f"  {trait}: {n_total} scoreable, {n_cached} cached")

    counters = {"ok": 0, "error": 0, "timeout": 0}
    scores = []

    for i, v in enumerate(scoreable):
        rsid = v["rsid"]
        if rsid in cache:
            entry = cache[rsid]
            if entry["expression_subscore"] is not None:
                scores.append(entry["expression_subscore"])
            continue

        chrom = v["chrom"] if v["chrom"].startswith("chr") else f"chr{v['chrom']}"
        try:
            vi = VariantInput(
                rsid=rsid,
                chrom=chrom,
                pos=int(v["pos"]),
                ref=v["ref"],
                alt=v["alt"],
                p_value=v.get("p_value"),
            )
        except Exception as e:
            _save_score(rsid, trait, v, None, f"input_error:{e}")
            counters["error"] += 1
            continue

        result = _score_one(model, vi, profile, rsid)
        error = result.get("error")
        score = result.get("expression_subscore")

        _save_score(rsid, trait, v, score, error)

        if error == "api_timeout":
            counters["timeout"] += 1
            print(f"    [{i+1}/{n_total}] {rsid}: TIMEOUT")
        elif error:
            counters["error"] += 1
            print(f"    [{i+1}/{n_total}] {rsid}: ERROR {error[:60]}")
        else:
            counters["ok"] += 1
            scores.append(score)
            if (i + 1) % 20 == 0 or i + 1 == n_total:
                pct_sat = (np.abs(scores) > 0.9).mean() if scores else 0.0
                print(f"    [{i+1}/{n_total}] {rsid}: {score:+.4f}  "
                      f"(running saturation: {pct_sat:.1%} |score|>0.9)")

    print(f"  {trait}: ok={counters['ok']} error={counters['error']} timeout={counters['timeout']}")
    return [s for s in scores if s is not None]


# ---------------------------------------------------------------------------
# Figure: saturation CDF (blood traits + Tewhey + disease GWAS)
# ---------------------------------------------------------------------------

def _ecdf(vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.abs(vals))
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def _load_existing_curves() -> tuple[np.ndarray, np.ndarray]:
    """Load Tewhey expression_subscore and disease GWAS expression signed_max_score."""
    import pandas as pd
    tewhey_df = pd.read_parquet(DIR / "tewhey_mpra.parquet")
    tewhey_arr = tewhey_df["expression_subscore"].dropna().to_numpy(float)

    conn = sqlite3.connect(DIR / "scored_variants.db")
    rows = conn.execute(
        "SELECT signed_max_score FROM modality_scores "
        "WHERE modality='expression' AND signed_max_score IS NOT NULL"
    ).fetchall()
    conn.close()
    gwas_arr = np.array([r[0] for r in rows], dtype=float)

    return tewhey_arr, gwas_arr


def make_phase3_cdf(blood_scores: dict[str, np.ndarray]) -> None:
    tewhey_arr, gwas_arr = _load_existing_curves()
    t_ref = np.linspace(0, 1, 300)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    fig.suptitle(
        "Saturation replication: blood cell trait GWAS variants\n"
        "(independent of Tewhey MPRA design)",
        fontsize=11, fontweight="bold",
    )

    # Uniform reference
    ax.plot(t_ref, t_ref, color=C_UNIF, lw=LW, ls="--",
            label="Uniform reference\n(random common variants, theoretical)",
            zorder=1)

    # Disease GWAS (existing)
    if len(gwas_arr) > 0:
        xg, yg = _ecdf(gwas_arr)
        ax.plot(xg, yg, color=C_GWAS, lw=LW,
                label=f"Disease GWAS (AD/T2D/SCZ/PD, n={len(gwas_arr)} pairs)",
                zorder=3)
        ax.fill_between(xg, 0, yg, color=C_GWAS, alpha=ALPHA_FILL)

    # Tewhey MPRA
    if len(tewhey_arr) > 0:
        xt, yt = _ecdf(tewhey_arr)
        ax.plot(xt, yt, color=C_TEWHEY, lw=LW,
                label=f"Tewhey MPRA loci (n={len(tewhey_arr)} GWAS regulatory variants)",
                zorder=4)
        ax.fill_between(xt, 0, yt, color=C_TEWHEY, alpha=ALPHA_FILL)

    # Blood trait GWAS (Phase 3)
    blood_colors = ["#762a83", "#9970ab", "#5e4fa2"]
    all_blood = []
    for (trait, arr), col in zip(blood_scores.items(), blood_colors):
        if len(arr) == 0:
            continue
        xb, yb = _ecdf(arr)
        label_name = trait.replace("_", " ").title()
        ax.plot(xb, yb, color=col, lw=LW,
                label=f"{label_name} GWAS (n={len(arr)} variants) [Phase 3]",
                zorder=5)
        ax.fill_between(xb, 0, yb, color=col, alpha=ALPHA_FILL)
        all_blood.extend(arr)

    ax.axvline(0.9, color="#888888", lw=0.8, ls=":", zorder=0)
    ax.text(0.905, 0.08, "|score|=0.9", fontsize=7.5, color="#666666", rotation=90)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("|Expression subscore| (quantile-normalized)", fontsize=10)
    ax.set_ylabel("Cumulative fraction of variants", fontsize=10)
    ax.grid(axis="both", lw=0.3, alpha=0.5)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

    # Gate 3 annotation
    if all_blood:
        pct = (np.abs(all_blood) > 0.9).mean()
        gate = "PASS" if pct >= GATE3_THRESHOLD else "FAIL"
        ax.text(0.98, 0.15,
                f"Blood trait variants:\n{pct:.1%} have |score|>0.9\nGate 3: {gate}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor=C_BLOOD, alpha=0.9))

    out = FIG_DIR / "phase3_saturation_cdf.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    _init_cache()

    print("=== Phase 3: Blood trait GWAS variant scoring ===\n")

    # Step 1: Fetch variants
    all_variants: dict[str, list[dict]] = {}
    for trait, efo_id in BLOOD_TRAITS.items():
        print(f"Fetching {trait} ({efo_id})...")
        variants = _fetch_trait_variants(trait, efo_id, max_n=TARGET_N)
        print(f"  Fetched {len(variants)} variants")
        all_variants[trait] = variants

    # Step 2: Resolve alleles (batch across both traits)
    print("\nResolving alleles via Ensembl...")
    for trait, variants in all_variants.items():
        print(f"  Resolving {len(variants)} rsIDs for {trait}...")
        all_variants[trait] = _resolve_alleles(variants)
        n_resolved = sum(1 for v in all_variants[trait] if v.get("ref") and v.get("alt"))
        print(f"  Resolved: {n_resolved}/{len(variants)}")

    # Step 3: Load AlphaGenome model
    print("\nLoading AlphaGenome model...")
    from alphagenome.models import dna_client
    model = dna_client.create(api_key=os.environ["ALPHAGENOME_API_KEY"])
    profile = _load_k562_profile()
    print("  Model loaded.")

    # Step 4: Score
    print("\nScoring variants (this will take several hours — resumable)...")
    blood_scores: dict[str, np.ndarray] = {}
    for trait, variants in all_variants.items():
        print(f"\n[{trait}]")
        scores = _score_trait(trait, variants, model, profile)
        arr = np.array([s for s in scores if s is not None], dtype=float)
        blood_scores[trait] = arr

        if len(arr) > 0:
            pct_sat = (np.abs(arr) > 0.9).mean()
            print(f"  {trait}: n={len(arr)}, |score|>0.9: {pct_sat:.1%}, "
                  f"median|score|={np.median(np.abs(arr)):.4f}")

    # Step 5: Gate 3 check
    print("\n=== Gate 3 Check ===")
    all_blood = np.concatenate(list(blood_scores.values())) if blood_scores else np.array([])
    if len(all_blood) > 0:
        pct = (np.abs(all_blood) > 0.9).mean()
        print(f"Blood trait variants: n={len(all_blood)}, |score|>0.9: {pct:.1%}")
        if pct >= GATE3_THRESHOLD:
            print(f"Gate 3: PASS (≥{GATE3_THRESHOLD:.0%} threshold)")
            print("→ Saturation generalizes: a property of regulatory-enriched selection.")
        else:
            print(f"Gate 3: FAIL (<{GATE3_THRESHOLD:.0%} threshold)")
            print("→ Saturation is Tewhey-specific. Frame paper accordingly.")
    else:
        print("No blood trait scores available. Check cache and API.")

    # Step 6: Figure
    print("\nGenerating Phase 3 saturation CDF...")
    make_phase3_cdf(blood_scores)
    print("\nPhase 3 complete. Figure in figures/phase3_saturation_cdf.png")


if __name__ == "__main__":
    main()
