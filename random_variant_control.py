"""Random variant negative control for saturation analysis.

Samples 1000 common variants (MAF > 0.01) from five chromosomal regions
with no regulatory pre-selection. Scores them with AlphaGenome's
quantile-normalized expression_subscore and caches results.

Expected result: scores distribute approximately uniformly — the
predictor does not saturate on unselected common variation, only on
regulatory-enriched sets. This is the negative control that demonstrates
saturation is a property of selection, not of the predictor itself.

Run:
    python random_variant_control.py
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
CACHE_DB = DIR / "random_variant_cache.db"
CONFIG_PATH = DIR / "config.yaml"

VARIANT_TIMEOUT_SECS = 60
TARGET_N = 1000
MAF_THRESHOLD = 0.01
BATCH_SIZE = 20       # small batches keep per-request time under 5s
SEED = 123

VALID_BASES = re.compile(r"^[ACGTNacgtn]+$")

# Five regions spread across the genome, away from well-known GWAS loci.
# 100kb windows; each yields ~4000+ biallelic SNVs, ~85% common.
# 50kb windows (Ensembl overlap limit); 10 regions for chromosomal diversity
SAMPLE_REGIONS = [
    ("1",  45_000_000,  45_050_000),
    ("2",  88_000_000,  88_050_000),
    ("5", 100_000_000, 100_050_000),
    ("7",  90_000_000,  90_050_000),
    ("10", 80_000_000,  80_050_000),
    ("12", 70_000_000,  70_050_000),
    ("15", 60_000_000,  60_050_000),
    ("17", 35_000_000,  35_050_000),
    ("19", 20_000_000,  20_050_000),
    ("20", 30_000_000,  30_050_000),
]
VARIANTS_PER_REGION = TARGET_N // len(SAMPLE_REGIONS)  # 100 each

ENSEMBL_OVERLAP = "https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}"
ENSEMBL_POST    = "https://rest.ensembl.org/variation/human?pops=1"

# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

def _init_cache() -> None:
    conn = sqlite3.connect(CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scores (
            rsid                TEXT PRIMARY KEY,
            chrom               TEXT,
            pos                 INTEGER,
            ref                 TEXT,
            alt                 TEXT,
            maf                 REAL,
            expression_subscore REAL,
            error               TEXT,
            scored_at           INTEGER
        )
    """)
    conn.commit()
    conn.close()


def _load_cache() -> dict[str, dict]:
    if not CACHE_DB.exists():
        return {}
    conn = sqlite3.connect(CACHE_DB)
    rows = conn.execute(
        "SELECT rsid, expression_subscore, error FROM scores"
    ).fetchall()
    conn.close()
    return {r[0]: {"expression_subscore": r[1], "error": r[2]} for r in rows}


def _save(rsid: str, v: dict, score: float | None, error: str | None) -> None:
    conn = sqlite3.connect(CACHE_DB)
    conn.execute(
        "INSERT OR REPLACE INTO scores "
        "(rsid,chrom,pos,ref,alt,maf,expression_subscore,error,scored_at) "
        "VALUES (?,?,?,?,?,?,?,?,?)",
        (rsid, v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"),
         v.get("maf"), score, error, int(time.time())),
    )
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Variant fetching + MAF resolution
# ---------------------------------------------------------------------------

def _fetch_region_rsids(chrom: str, start: int, end: int) -> list[str]:
    """Return rsIDs of biallelic SNVs in the region."""
    url = ENSEMBL_OVERLAP.format(chrom=chrom, start=start, end=end)
    resp = requests.get(url, params={"feature": "variation"},
                        headers={"Content-Type": "application/json",
                                 "Accept": "application/json"}, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    rsids = [
        v["id"] for v in data
        if str(v.get("id", "")).startswith("rs")
        and len(v.get("alleles", [])) == 2
        and all(len(a) == 1 for a in v.get("alleles", []))
    ]
    return rsids


def _resolve_batch(rsids: list[str]) -> list[dict]:
    """Resolve a small batch with MAF; returns list of variant dicts."""
    resp = requests.post(
        ENSEMBL_POST,
        json={"ids": rsids},
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    variants = []
    for rsid in rsids:
        info = data.get(rsid)
        if not info:
            continue

        # Get GRCh38 mapping
        mappings = info.get("mappings", [])
        grch38 = [m for m in mappings if m.get("assembly_name") == "GRCh38"]
        mapping = grch38[0] if grch38 else (mappings[0] if mappings else None)
        if not mapping:
            continue

        chrom_raw = str(mapping.get("seq_region_name", ""))
        chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
        pos = mapping.get("start")
        allele_str = mapping.get("allele_string", "")
        if not allele_str or "/" not in allele_str:
            continue
        parts = allele_str.split("/")
        ref = parts[0]
        alts = [a for a in parts[1:] if a and a != ref and len(a) == 1]
        if not alts or len(ref) != 1:
            continue
        if not re.match(r"^[ACGTacgt]$", ref) or not re.match(r"^[ACGTacgt]$", alts[0]):
            continue

        # MAF from populations
        pops = info.get("populations", [])
        maf = max((p.get("frequency", 0.0) for p in pops if p.get("allele") != ref),
                  default=0.0)
        if maf < MAF_THRESHOLD:
            continue

        alt = alts[0]
        variants.append({"rsid": rsid, "chrom": chrom, "pos": int(pos),
                         "ref": ref, "alt": alt, "maf": maf})
    return variants


def _sample_region(chrom: str, start: int, end: int,
                   n_target: int) -> list[dict]:
    """Fetch and resolve up to n_target common SNVs from one region."""
    print(f"  Fetching {chrom}:{start:,}-{end:,}...")
    rsids = _fetch_region_rsids(chrom, start, end)
    rng = np.random.default_rng(SEED)
    rng.shuffle(rsids)  # type: ignore[arg-type]

    variants: list[dict] = []
    for i in range(0, len(rsids), BATCH_SIZE):
        if len(variants) >= n_target:
            break
        batch = rsids[i:i + BATCH_SIZE]
        try:
            time.sleep(1.0 / 15)
            resolved = _resolve_batch(batch)
            variants.extend(resolved)
        except Exception as e:
            print(f"    batch error: {e}")

    result = variants[:n_target]
    print(f"  {chrom}: {len(result)} common SNVs resolved")
    return result


# ---------------------------------------------------------------------------
# Scoring
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
    mb = result.get("modality_breakdown")
    if mb is None or mb.empty:
        return None
    expr = mb[mb["modality"] == "expression"]
    if expr.empty:
        return None
    val = expr.iloc[0].get("signed_max_score")
    return float(val) if val is not None else None


def _score_one(model, vi, profile) -> dict:
    from scoring.composite import score_single_variant
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    future = executor.submit(score_single_variant, model, vi, profile)
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


def _score_variants(variants: list[dict], model, profile) -> list[float]:
    from scoring.composite import VariantInput

    cache = _load_cache()
    n_total = len(variants)
    n_cached = sum(1 for v in variants if v["rsid"] in cache
                   and cache[v["rsid"]]["expression_subscore"] is not None)
    print(f"  Scoring {n_total} variants ({n_cached} cached)...")

    scores = []
    counters = {"ok": 0, "error": 0, "timeout": 0}

    for i, v in enumerate(variants):
        rsid = v["rsid"]
        if rsid in cache:
            s = cache[rsid]["expression_subscore"]
            if s is not None:
                scores.append(s)
            continue

        if not (VALID_BASES.match(v["ref"]) and VALID_BASES.match(v["alt"])):
            _save(rsid, v, None, "invalid_bases")
            continue

        try:
            vi = VariantInput(rsid=rsid, chrom=v["chrom"], pos=v["pos"],
                              ref=v["ref"], alt=v["alt"], maf=v.get("maf"))
        except Exception as e:
            _save(rsid, v, None, f"input_error:{e}")
            continue

        result = _score_one(model, vi, profile)
        error = result.get("error")
        score = result.get("expression_subscore")

        _save(rsid, v, score, error)

        if error == "api_timeout":
            counters["timeout"] += 1
        elif error:
            counters["error"] += 1
        else:
            counters["ok"] += 1
            if score is not None:
                scores.append(score)

        if (i + 1) % 100 == 0 or i + 1 == n_total:
            pct_sat = (np.abs(scores) > 0.9).mean() if scores else 0.0
            print(f"  [{i+1}/{n_total}] ok={counters['ok']} "
                  f"saturation so far: {pct_sat:.1%} |score|>0.9")

    print(f"  Done: ok={counters['ok']} error={counters['error']} "
          f"timeout={counters['timeout']}")
    return scores


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _ecdf(vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.abs(vals))
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def make_updated_saturation_figure(random_scores: np.ndarray) -> None:
    """Regenerate saturation_cdf.png with random variant control added."""
    import pandas as pd

    tewhey_df = pd.read_parquet(DIR / "tewhey_mpra.parquet")
    tewhey_arr = tewhey_df["expression_subscore"].dropna().to_numpy(float)

    import sqlite3 as sq
    conn = sq.connect(DIR / "scored_variants.db")
    rows = conn.execute(
        "SELECT signed_max_score FROM modality_scores "
        "WHERE modality='expression' AND signed_max_score IS NOT NULL"
    ).fetchall()
    conn.close()
    gwas_arr = np.array([r[0] for r in rows], dtype=float)

    t_ref = np.linspace(0, 1, 300)
    C_TEWHEY = "#2166ac"
    C_GWAS   = "#d6604d"
    C_RANDOM = "#4dac26"   # green — negative control
    C_UNIF   = "#999999"
    LW = 1.8
    ALPHA = 0.08

    fig, ax = plt.subplots(figsize=(8, 5.5))
    fig.suptitle(
        "Quantile-normalized AlphaGenome scores saturate on regulatory-enriched variant sets\n"
        "but not on random common variation",
        fontsize=10.5, fontweight="bold",
    )

    ax.plot(t_ref, t_ref, color=C_UNIF, lw=LW, ls="--",
            label="Theoretical uniform (AlphaGenome calibration construction)",
            zorder=1)

    xr, yr = _ecdf(random_scores)
    ax.plot(xr, yr, color=C_RANDOM, lw=LW,
            label=f"Random common variants — negative control\n"
                  f"(n={len(random_scores)}, genome-sampled, no regulatory pre-selection)",
            zorder=2)
    ax.fill_between(xr, 0, yr, color=C_RANDOM, alpha=ALPHA)

    xg, yg = _ecdf(gwas_arr)
    ax.plot(xg, yg, color=C_GWAS, lw=LW,
            label=f"Disease GWAS variants (n={len(gwas_arr)} pairs, AD/T2D/SCZ/PD)",
            zorder=3)
    ax.fill_between(xg, 0, yg, color=C_GWAS, alpha=ALPHA)

    xt, yt = _ecdf(tewhey_arr)
    ax.plot(xt, yt, color=C_TEWHEY, lw=LW,
            label=f"Tewhey 2016 MPRA regulatory loci (n={len(tewhey_arr)})",
            zorder=4)
    ax.fill_between(xt, 0, yt, color=C_TEWHEY, alpha=ALPHA)

    ax.axvline(0.9, color="#888888", lw=0.8, ls=":", zorder=0)
    ax.text(0.905, 0.08, "|score| = 0.9", fontsize=7.5, color="#555555", rotation=90)

    # Saturation annotations
    pct_rand = (np.abs(random_scores) > 0.9).mean()
    ax.text(0.97, 0.55,
            f"Random: {pct_rand:.0%} >0.9\n"
            f"GWAS: {(np.abs(gwas_arr)>0.9).mean():.0%} >0.9\n"
            f"Tewhey: {(np.abs(tewhey_arr)>0.9).mean():.0%} >0.9",
            transform=ax.transAxes, ha="right", va="center",
            fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#aaaaaa", alpha=0.9))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("|Expression subscore| (quantile-normalized)", fontsize=10)
    ax.set_ylabel("Cumulative fraction of variants", fontsize=10)
    ax.grid(axis="both", lw=0.3, alpha=0.5)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

    out = FIG_DIR / "saturation_cdf.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved updated → {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    _init_cache()
    print("=== Random Variant Negative Control ===\n")

    # Step 1: Collect common variants from random regions
    print("Sampling variants from random genomic regions...")
    all_variants: list[dict] = []
    for chrom, start, end in SAMPLE_REGIONS:
        region_variants = _sample_region(chrom, start, end, VARIANTS_PER_REGION)
        all_variants.extend(region_variants)
    print(f"\nTotal collected: {len(all_variants)} common SNVs\n")

    # Step 2: Load model and score
    print("Loading AlphaGenome model...")
    from alphagenome.models import dna_client
    model = dna_client.create(api_key=os.environ["ALPHAGENOME_API_KEY"])
    profile = _load_k562_profile()
    print("Model loaded.\n")

    scores = _score_variants(all_variants, model, profile)
    arr = np.array([s for s in scores if s is not None], dtype=float)

    print(f"\n=== Result ===")
    print(f"n={len(arr)}, |score|>0.9: {(np.abs(arr)>0.9).mean():.1%}, "
          f"median|score|={np.median(np.abs(arr)):.4f}")

    if (np.abs(arr) > 0.9).mean() > 0.50:
        print("WARNING: saturation rate unexpectedly high for random variants.")
        print("Check that regions are not near known regulatory loci.")
    else:
        print("Negative control holds: random variants do not saturate.")

    print("\nUpdating saturation figure with negative control...")
    make_updated_saturation_figure(arr)
    print("\nDone.")


if __name__ == "__main__":
    main()
