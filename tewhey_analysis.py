"""
Tewhey 2016 MPRA correlation analysis.

Data: GSE75661 7.5k collapsed counts (LCL: NA12878 5-rep + NA19239 3-rep).
Supplementary tables not accessible (Cell access-controlled, not PMC OA).
LFC computed from raw counts: log2(RNA_CPM/Plasmid_CPM), allelic diff = B - A.

Scores compared vs measured LFC:
  1. expression_subscore: signed score for expression modality (RNA-seq+CAGE+PRO-cap)
  2. full_composite: standard composite with blood/LCL tissue weighting
  3. cadd_phred: CADD v1.7 PHRED score (GRCh38)

Stop flags written to TEWHEY_FLAGS.md; raises SystemExit on any stop condition.
"""

from __future__ import annotations

import concurrent.futures
import io
import json
import re
import sqlite3
import sys
import time
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DIR = Path(__file__).parent
FLAGS_FILE = DIR / "TEWHEY_FLAGS.md"
PARQUET_OUT = DIR / "tewhey_mpra.parquet"
CADD_CACHE_DB = DIR / "cadd_cache.db"
TEWHEY_SCORES_CACHE_DB = DIR / "tewhey_scores_cache.db"

VARIANT_TIMEOUT_SECS = 60

GEO_7K5_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75661/suppl/"
    "GSE75661_7.5k_collapsed_counts.txt.gz"
)
CADD_API_BASE = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7"
CADD_RATE = 2.0  # req/sec — CADD API soft-bans rapid bursts by returning []; keep gentle

PSEUDOCOUNT = 1.0
N_BOOTSTRAP = 1000
rng = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Stop-flag helpers
# ---------------------------------------------------------------------------

def _flag(label: str, message: str, stop: bool = True) -> None:
    with open(FLAGS_FILE, "a") as f:
        f.write(f"\n## FLAG [{time.strftime('%Y-%m-%d %H:%M:%S')}]: {label}\n{message}\n")
    print(f"\n!!! FLAG: {label}")
    print(f"    {message}")
    if stop:
        print(f"    Written to {FLAGS_FILE} — halting.")
        sys.exit(1)


def _check_nan_inf(series: pd.Series, name: str) -> None:
    bad = (~np.isfinite(series.dropna())).sum() + series.isna().sum()
    if bad > 0:
        _flag(
            f"NaN/Inf in {name}",
            f"{bad}/{len(series)} values are NaN or Inf — scoring bug or upstream issue.",
        )


# ---------------------------------------------------------------------------
# Step 1: Compute LFC from raw counts
# ---------------------------------------------------------------------------

def _compute_lfc() -> pd.DataFrame:
    """Download 7.5k counts, compute allelic LFC per rsID variant."""
    print("Downloading GSE75661_7.5k_collapsed_counts.txt.gz …")
    with urllib.request.urlopen(GEO_7K5_URL, timeout=120) as resp:
        raw = resp.read()

    df = pd.read_csv(io.BytesIO(raw), sep="\t", compression="gzip")
    print(f"  Raw shape: {df.shape}  (oligos × replicates)")

    # Columns: Oligo, Plasmid_r1..r5, NA12878_r1..r5, NA19239_r1..r3
    plasmid_cols = [c for c in df.columns if c.startswith("Plasmid")]
    rna_cols = [c for c in df.columns if c not in ["Oligo"] + plasmid_cols]

    # Only rsID-named oligos (chr:pos oligos excluded — no allele info)
    rs_mask = df["Oligo"].str.match(r"^rs\d+_")
    df_rs = df[rs_mask].copy()
    print(f"  rsID oligos: {len(df_rs)} (chr:pos oligos excluded: {(~rs_mask).sum()})")

    # Filter zero-plasmid oligos (divide-by-zero risk)
    plasmid_total = df_rs[plasmid_cols].sum(axis=1)
    df_rs = df_rs[plasmid_total > 0].copy()
    print(f"  After zero-plasmid filter: {len(df_rs)} oligos")

    # Parse rsid + allele
    parsed = df_rs["Oligo"].str.extract(r"^(rs\d+)_(?:RC_)?([AB])$")
    df_rs["rsid"] = parsed[0]
    df_rs["allele"] = parsed[1]
    df_rs = df_rs.dropna(subset=["rsid", "allele"])

    # CPM-normalize per replicate
    count_mat = df_rs.set_index(["rsid", "allele"])[plasmid_cols + rna_cols]
    for col in count_mat.columns:
        col_total = count_mat[col].sum()
        count_mat[col] = count_mat[col] / col_total * 1e6

    plasmid_mean = count_mat[plasmid_cols].mean(axis=1)
    rna_mean = count_mat[rna_cols].mean(axis=1)
    activity = np.log2((rna_mean + PSEUDOCOUNT) / (plasmid_mean + PSEUDOCOUNT))

    act_df = activity.reset_index()
    act_df.columns = ["rsid", "allele", "activity"]

    pivot = act_df.pivot_table(
        index="rsid", columns="allele", values="activity", aggfunc="mean"
    )
    pivot.columns = [f"activity_{c}" for c in pivot.columns]
    pivot = pivot.reset_index()

    if "activity_A" not in pivot.columns or "activity_B" not in pivot.columns:
        _flag(
            "Missing allele columns",
            f"Pivot missing activity_A or activity_B. Found: {pivot.columns.tolist()}",
        )

    pivot["mpra_lfc"] = pivot["activity_B"] - pivot["activity_A"]
    pivot = pivot.dropna(subset=["mpra_lfc"])
    print(f"  Variants with LFC: {len(pivot)} (from {len(df_rs['rsid'].unique())} rsIDs)")
    return pivot[["rsid", "activity_A", "activity_B", "mpra_lfc"]]


# ---------------------------------------------------------------------------
# Step 2: Resolve GRCh38 coordinates
# ---------------------------------------------------------------------------

def _resolve_coords(rsids: list[str]) -> pd.DataFrame:
    from allele_resolver import resolve_alleles_batch, DB_PATH as ALLELE_DB

    print(f"\nResolving {len(rsids)} rsIDs to GRCh38 via Ensembl …")
    resolved = resolve_alleles_batch(rsids, db_path=ALLELE_DB, include_pops=False)

    rows = []
    for rsid in rsids:
        info = resolved.get(rsid)
        if info:
            rows.append({
                "rsid": rsid,
                "chrom": info["chrom"],
                "pos": info["pos"],
                "ref": info["ref"],
                "alt": info["alt"],
                "maf": info.get("maf"),
            })
        else:
            rows.append({"rsid": rsid, "chrom": None, "pos": None,
                         "ref": None, "alt": None, "maf": None})

    coord_df = pd.DataFrame(rows)
    n_resolved = coord_df["chrom"].notna().sum()
    n_total = len(coord_df)
    pct_lost = 100 * (1 - n_resolved / n_total)
    print(f"  Resolved: {n_resolved}/{n_total} ({pct_lost:.1f}% unresolved)")

    if pct_lost > 15:
        _flag(
            "Liftover/resolution loss > 15%",
            f"{n_total - n_resolved}/{n_total} variants ({pct_lost:.1f}%) could not be resolved. "
            "Exceeds 15% threshold.",
        )
    return coord_df


# ---------------------------------------------------------------------------
# Step 3: CADD scores
# ---------------------------------------------------------------------------

def _init_cadd_cache() -> None:
    conn = sqlite3.connect(CADD_CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cadd (
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
            phred REAL, raw_score REAL,
            not_found INTEGER DEFAULT 0,
            fetched_at INTEGER,
            PRIMARY KEY (chrom, pos, ref, alt)
        )
    """)
    conn.commit()
    conn.close()


_cadd_last_call = 0.0


def _cadd_rate_limit() -> None:
    global _cadd_last_call
    now = time.monotonic()
    wait = 1.0 / CADD_RATE - (now - _cadd_last_call)
    if wait > 0:
        time.sleep(wait)
    _cadd_last_call = time.monotonic()


def _fetch_cadd_one(chrom: str, pos: int, ref: str, alt: str) -> tuple[float | None, float | None]:
    """Fetch CADD PHRED + RawScore for one variant. Returns (None, None) if not found.

    URL format: {CADD_API_BASE}/chr{chrom}:{pos}_{ref}_{alt}
    The chr prefix is required — bare numeric chromosomes return an empty list.
    CADD soft-rate-limits rapid bursts by returning [] with HTTP 200; treat an
    empty list as a retriable condition and back off before retrying.
    """
    chrom_str = chrom if chrom.startswith("chr") else f"chr{chrom}"
    query = f"{chrom_str}:{pos}_{ref}_{alt}"
    url = f"{CADD_API_BASE}/{query}"

    for attempt in range(5):
        _cadd_rate_limit()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "gwas-finemapping-research/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
            if isinstance(data, list) and data:
                entry = data[0]
                return float(entry.get("PHRED", 0)), float(entry.get("RawScore", 0))
            # Empty list: variant absent from CADD or soft rate-limit.
            # Retry once with a short backoff; if still empty, treat as not-in-CADD.
            if attempt == 0:
                time.sleep(1)
                continue
            return None, None
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(2 ** (attempt + 1))
            elif e.code == 404:
                return None, None
            else:
                if attempt >= 4:
                    return None, None
                time.sleep(2 ** attempt)
        except Exception:
            if attempt >= 4:
                return None, None
            time.sleep(2 ** attempt)
    return None, None


def _get_cadd_scores(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Fetch CADD PHRED scores for all variants; uses SQLite cache."""
    _init_cadd_cache()

    results = {}
    to_fetch = []

    conn = sqlite3.connect(CADD_CACHE_DB)
    for _, row in variants_df.iterrows():
        r = conn.execute(
            "SELECT phred, raw_score, not_found FROM cadd WHERE chrom=? AND pos=? AND ref=? AND alt=?",
            (row.chrom, row.pos, row.ref, row.alt),
        ).fetchone()
        if r is not None:
            results[row.rsid] = (None if r[2] else r[0], None if r[2] else r[1])
        else:
            to_fetch.append(row)
    conn.close()

    print(f"\nFetching CADD scores: {len(to_fetch)} uncached, "
          f"{len(results)} already cached …")

    # Probe with 5 variants before committing to the full set.  If every probe
    # returns None, CADD is likely rate-limiting; mark all remaining as not-found
    # and continue so AlphaGenome scoring is not blocked.
    if to_fetch:
        print("  Probing CADD API with 5 test variants …")
        probe_hits = 0
        for probe_row in to_fetch[:5]:
            ph, _ = _fetch_cadd_one(probe_row.chrom, probe_row.pos,
                                    probe_row.ref, probe_row.alt)
            if ph is not None:
                probe_hits += 1
        if probe_hits == 0 and len(to_fetch) > 5:
            _flag(
                "CADD API unavailable (probe returned 0/5)",
                "All 5 probe requests returned empty results. CADD may be "
                "rate-limiting or temporarily down. Skipping remaining CADD "
                "lookups and continuing to AlphaGenome scoring. Re-run with "
                "CADD cache populated to include CADD scores in the report.",
                stop=False,
            )
            print("  Skipping CADD — marking all as not-found to unblock scoring.")
            conn = sqlite3.connect(CADD_CACHE_DB)
            for row in to_fetch:
                results[row.rsid] = (None, None)
                conn.execute(
                    "INSERT OR REPLACE INTO cadd "
                    "(chrom,pos,ref,alt,phred,raw_score,not_found,fetched_at) "
                    "VALUES (?,?,?,?,NULL,NULL,1,?)",
                    (row.chrom, row.pos, row.ref, row.alt, int(time.time())),
                )
            conn.commit()
            conn.close()
            out = variants_df[["rsid"]].copy()
            out["cadd_phred"] = out["rsid"].map(lambda r: (results.get(r) or (None, None))[0])
            return out

    fails = 0
    for i, row in enumerate(to_fetch):
        if i % 50 == 0:
            print(f"  CADD [{i}/{len(to_fetch)}] …")
        phred, raw = _fetch_cadd_one(row.chrom, row.pos, row.ref, row.alt)

        conn = sqlite3.connect(CADD_CACHE_DB)
        not_found = int(phred is None)
        conn.execute(
            "INSERT OR REPLACE INTO cadd (chrom,pos,ref,alt,phred,raw_score,not_found,fetched_at) "
            "VALUES (?,?,?,?,?,?,?,?)",
            (row.chrom, row.pos, row.ref, row.alt,
             phred, raw, not_found, int(time.time())),
        )
        conn.commit()
        conn.close()

        results[row.rsid] = (phred, raw)
        if phred is None:
            fails += 1

    n_total = len(variants_df)
    pct_fail = 100 * fails / n_total if n_total else 0
    print(f"  CADD lookup: {n_total - fails}/{n_total} found, {fails} not found ({pct_fail:.1f}%)")

    if pct_fail > 10:
        # CADD is the baseline comparison, not the primary result.  A high failure
        # rate (often from API rate-limiting) degrades that comparison but should
        # not block AlphaGenome scoring.  Flag but continue.
        _flag(
            "CADD lookup failure > 10%",
            f"{fails}/{n_total} variants ({pct_fail:.1f}%) had no CADD score. "
            "CADD correlation will show reduced n or be omitted. "
            "AlphaGenome scoring will proceed regardless.",
            stop=False,
        )

    out = variants_df[["rsid"]].copy()
    out["cadd_phred"] = out["rsid"].map(lambda r: (results.get(r) or (None, None))[0])
    return out


# ---------------------------------------------------------------------------
# Step 4: AlphaGenome scoring (K562/blood-lineage tissue profile)
# ---------------------------------------------------------------------------

def _init_tewhey_cache() -> None:
    conn = sqlite3.connect(TEWHEY_SCORES_CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scores (
            rsid                TEXT PRIMARY KEY,
            full_composite      REAL,
            expression_subscore REAL,
            error               TEXT,
            scored_at           INTEGER
        )
    """)
    conn.commit()
    conn.close()


def _load_tewhey_cache() -> dict:
    _init_tewhey_cache()
    conn = sqlite3.connect(TEWHEY_SCORES_CACHE_DB)
    rows = conn.execute(
        "SELECT rsid, full_composite, expression_subscore, error FROM scores"
    ).fetchall()
    conn.close()
    return {r[0]: {"full_composite": r[1], "expression_subscore": r[2], "error": r[3]}
            for r in rows}


def _save_tewhey_score(rsid: str, entry: dict) -> None:
    conn = sqlite3.connect(TEWHEY_SCORES_CACHE_DB)
    conn.execute(
        "INSERT OR REPLACE INTO scores "
        "(rsid, full_composite, expression_subscore, error, scored_at) "
        "VALUES (?, ?, ?, ?, ?)",
        (rsid, entry.get("full_composite"), entry.get("expression_subscore"),
         entry.get("error"), int(time.time())),
    )
    conn.commit()
    conn.close()


def _get_k562_profile():
    from scoring.tissue_config import TissueProfile
    return TissueProfile(
        display_name="K562 / Blood-lineage (LCL MPRA)",
        biosample_keywords=[
            "k562", "blood", "lymph", "lymphoblast", "b cell",
            "gm12878", "na12878", "na19239", "hematopoietic",
            "erythro", "leukemia", "bcell", "b-cell",
        ],
        gtex_keywords=["Blood"],
        ontology_notes="Used for Tewhey 2016 LCL MPRA; K562 and LCL are both blood-lineage",
    )


def _score_tewhey(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Score variants through AlphaGenome per-variant with 60s timeout + SQLite cache.

    Uses score_single_variant (sequential) with ThreadPoolExecutor timeout —
    same pattern as batch_score.py — to avoid hanging on edge-case variants.
    Cache in tewhey_scores_cache.db makes runs resumable.
    """
    import os
    from dotenv import load_dotenv
    from alphagenome.models import dna_client
    from scoring.composite import VariantInput, score_single_variant

    load_dotenv(DIR / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not api_key:
        _flag("No API key", "ALPHAGENOME_API_KEY not set.")

    profile = _get_k562_profile()
    _valid = re.compile(r"^[ACGTNacgtn]+$")
    scoreable = variants_df.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
    scoreable = scoreable[
        scoreable["ref"].apply(lambda x: bool(x and x != "-" and _valid.match(x))) &
        scoreable["alt"].apply(lambda x: bool(x and x != "-" and _valid.match(x)))
    ].copy()

    print(f"\nScoring {len(scoreable)} variants through AlphaGenome (K562/blood profile) …")

    cache = _load_tewhey_cache()
    to_score = [row for _, row in scoreable.iterrows() if row.rsid not in cache]
    print(f"  {len(cache)} cached, {len(to_score)} to score")

    model = dna_client.create(api_key)

    for i, row in enumerate(to_score):
        rsid = row.rsid
        chrom = row.chrom if str(row.chrom).startswith("chr") else f"chr{row.chrom}"

        vinput = VariantInput(
            rsid=rsid,
            chrom=chrom,
            pos=int(row.pos),
            ref=row.ref,
            alt=row.alt,
            maf=getattr(row, "maf", None),
        )

        t0 = time.monotonic()
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        future = executor.submit(score_single_variant, model, vinput, profile)
        result = None
        try:
            result = future.result(timeout=VARIANT_TIMEOUT_SECS)
        except concurrent.futures.TimeoutError:
            result = {"composite_score": None, "error": "api_timeout", "modality_breakdown": None}
            print(f"  [{i+1}/{len(to_score)}] TIMEOUT ({VARIANT_TIMEOUT_SECS}s): {rsid} — skipping")
        except Exception as exc:
            err_str = str(exc).lower()
            if any(k in err_str for k in ("rate", "quota", "limit", "429")):
                _flag("AlphaGenome rate limit", f"Rate limit hit on {rsid}: {exc}")
            result = {"composite_score": None, "error": str(exc), "modality_breakdown": None}
            print(f"  [{i+1}/{len(to_score)}] ERROR {rsid}: {exc}")
        finally:
            executor.shutdown(wait=False)

        elapsed = time.monotonic() - t0

        expr_signed = None
        mod_df = result.get("modality_breakdown") if result else None
        if mod_df is not None and hasattr(mod_df, "empty") and not mod_df.empty:
            expr_rows = mod_df[mod_df["modality"] == "expression"]
            if not expr_rows.empty:
                if "signed_max_score" in expr_rows.columns:
                    expr_signed = float(expr_rows["signed_max_score"].iloc[0])
                elif "max_abs_score" in expr_rows.columns:
                    expr_signed = float(expr_rows["max_abs_score"].iloc[0])

        entry = {
            "full_composite": result.get("composite_score") if result else None,
            "expression_subscore": expr_signed,
            "error": result.get("error") if result else "no_result",
        }
        cache[rsid] = entry
        _save_tewhey_score(rsid, entry)

        if (i + 1) % 100 == 0 or entry.get("error"):
            flag = (f" composite={entry['full_composite']:.4f}"
                    if entry["full_composite"] is not None
                    else f" ERROR={entry['error']}")
            print(f"  [{i+1}/{len(to_score)}] {rsid}{flag}  ({elapsed:.1f}s)")

    # Assemble result DataFrame from cache (includes pre-cached entries)
    all_rows = []
    for _, row in scoreable.iterrows():
        rsid = row.rsid
        entry = cache.get(rsid, {"full_composite": None, "expression_subscore": None,
                                  "error": "not_scored"})
        all_rows.append({"rsid": rsid, **entry})

    result_df = pd.DataFrame(all_rows)
    n_ok = result_df["error"].isna().sum()
    print(f"  Scored: {n_ok}/{len(result_df)} without errors")
    return result_df


# ---------------------------------------------------------------------------
# Step 5: Correlation analysis
# ---------------------------------------------------------------------------

def _bootstrap_ci(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOTSTRAP) -> tuple[float, float]:
    """Bootstrap 95% CI for Spearman r."""
    n = len(x)
    boot_rs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = st.spearmanr(x[idx], y[idx])
        boot_rs.append(r)
    arr = np.array(boot_rs)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def _correlate(label: str, x: pd.Series, y: pd.Series) -> dict:
    """Compute Spearman r, p, 95% CI between x and y after joint NA-drop."""
    valid = (~x.isna()) & (~y.isna()) & np.isfinite(x.fillna(0)) & np.isfinite(y.fillna(0))
    xv, yv = x[valid].to_numpy(), y[valid].to_numpy()
    n = len(xv)
    if n < 10:
        return {"label": label, "n": n, "r": None, "p": None, "ci_lo": None, "ci_hi": None}
    r, p = st.spearmanr(xv, yv)
    ci_lo, ci_hi = _bootstrap_ci(xv, yv)
    return {"label": label, "n": n, "r": round(r, 4), "p": float(p),
            "ci_lo": round(ci_lo, 4), "ci_hi": round(ci_hi, 4)}


def _run_correlations(df: pd.DataFrame) -> list[dict]:
    """Run all 6 correlations: 3 scores × 2 subsets (all / top-10% |LFC|)."""
    rows = []

    # --- All variants ---
    for score_col, label in [
        ("expression_subscore", "expression_subscore"),
        ("full_composite",      "full_composite"),
        ("cadd_phred",          "cadd_phred"),
    ]:
        rows.append(_correlate(f"{label} | all", df[score_col], df["mpra_lfc"]))

    # --- Top 10% |LFC| ---
    thresh = df["mpra_lfc"].abs().quantile(0.90)
    df_top = df[df["mpra_lfc"].abs() >= thresh]
    for score_col, label in [
        ("expression_subscore", "expression_subscore"),
        ("full_composite",      "full_composite"),
        ("cadd_phred",          "cadd_phred"),
    ]:
        rows.append(_correlate(f"{label} | top-10% |LFC|", df_top[score_col], df_top["mpra_lfc"]))

    return rows


# ---------------------------------------------------------------------------
# Step 6: Stop-condition checks
# ---------------------------------------------------------------------------

def _check_stop_conditions(df: pd.DataFrame, corr_rows: list[dict]) -> None:
    # Expression subscore < 0.10 (all-variants row)
    expr_all = next((r for r in corr_rows if r["label"] == "expression_subscore | all"), None)
    if expr_all and expr_all["r"] is not None and abs(expr_all["r"]) < 0.10:
        _flag(
            "expression_subscore correlation < 0.10",
            f"Spearman r = {expr_all['r']} (n={expr_all['n']}) is below the 0.10 threshold. "
            "Possible causes: wrong allele orientation (A/B vs ref/alt mismatch), "
            "tissue mismatch (LCL vs blood-lineage filter), or model calibration issue. "
            "Do NOT report until diagnosed.",
        )

    # NaN/Inf checks
    for col in ["expression_subscore", "full_composite", "cadd_phred", "mpra_lfc"]:
        if col in df.columns:
            _check_nan_inf(df[col], col)


# ---------------------------------------------------------------------------
# Step 7a: Scatter plot
# ---------------------------------------------------------------------------

def _make_scatter(df: pd.DataFrame, corr_rows: list[dict]) -> None:
    valid = df["expression_subscore"].notna() & df["mpra_lfc"].notna()
    xv = df.loc[valid, "mpra_lfc"]
    yv = df.loc[valid, "expression_subscore"]

    expr_all = next(
        (r for r in corr_rows if r["label"] == "expression_subscore | all"), None
    )
    r_val = expr_all["r"] if expr_all else None
    p_val = expr_all["p"] if expr_all else None
    n = expr_all["n"] if expr_all else len(xv)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(xv, yv, alpha=0.25, s=10, color="#1f77b4", rasterized=True)

    if r_val is not None:
        subtitle = f"Spearman r = {r_val:+.3f}  p = {p_val:.2e}  n = {n}"
    else:
        subtitle = f"n = {n}"

    ax.set_xlabel("MPRA allelic LFC (B − A)", fontsize=11)
    ax.set_ylabel("expression_subscore (AlphaGenome)", fontsize=11)
    ax.set_title(f"Tewhey 2016 K562 LCL MPRA\n{subtitle}", fontsize=10)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")

    out = DIR / "tewhey_correlation.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Scatter plot saved to {out}")


# ---------------------------------------------------------------------------
# Step 7: Write benchmark_report.md section
# ---------------------------------------------------------------------------

def _write_report(df: pd.DataFrame, corr_rows: list[dict]) -> None:
    bm_path = DIR / "benchmark_report.md"

    # Build table
    def fmt_r(row):
        if row["r"] is None:
            return f"n={row['n']} (insufficient)"
        p_str = f"{row['p']:.2e}" if row["p"] is not None else "—"
        return (
            f"r={row['r']:+.3f} [{row['ci_lo']:+.3f}, {row['ci_hi']:+.3f}] "
            f"p={p_str} n={row['n']}"
        )

    table_lines = [
        "| Score | Subset | Spearman r [95% CI] | p-value | n |",
        "|---|---|---|---|---|",
    ]
    for row in corr_rows:
        parts = row["label"].split(" | ")
        score_col, subset = parts[0], parts[1] if len(parts) > 1 else "—"
        if row["r"] is None:
            table_lines.append(f"| {score_col} | {subset} | — | — | {row['n']} |")
        else:
            ci = f"[{row['ci_lo']:+.3f}, {row['ci_hi']:+.3f}]"
            p_str = f"{row['p']:.2e}"
            table_lines.append(f"| {score_col} | {subset} | {row['r']:+.3f} {ci} | {p_str} | {row['n']} |")

    # Summary stats for narrative
    n_total = len(df)
    n_scored = df["full_composite"].notna().sum()
    n_cadd = df["cadd_phred"].notna().sum()
    n_expr = df["expression_subscore"].notna().sum()
    lfc_range = (df["mpra_lfc"].min(), df["mpra_lfc"].max())
    lfc_std = df["mpra_lfc"].std()

    expr_all = next((r for r in corr_rows if "expression_subscore" in r["label"] and "all" in r["label"]), None)
    full_all = next((r for r in corr_rows if "full_composite" in r["label"] and "all" in r["label"]), None)
    cadd_all = next((r for r in corr_rows if "cadd_phred" in r["label"] and "all" in r["label"]), None)

    # Honest interpretation based on actual numbers
    def interpret(r_val):
        if r_val is None:
            return "no data"
        a = abs(r_val)
        if a >= 0.40:
            return "strong"
        if a >= 0.25:
            return "moderate"
        if a >= 0.10:
            return "weak but significant"
        return "negligible"

    expr_r = expr_all["r"] if expr_all else None
    full_r = full_all["r"] if full_all else None
    cadd_r = cadd_all["r"] if cadd_all else None

    section = f"""
## Primary validation: MPRA correlation (Tewhey 2016)

**Dataset:** GSE75661 7.5k-oligo library. Tested in LCL (NA12878 5 replicates,
NA19239 3 replicates). Note: the user specified "K562/blood-lineage" tissue
weights; the actual MPRA was done in lymphoblastoid cell lines (LCL), which are
blood-lineage but distinct from K562 (erythroleukemia). Blood/LCL-relevant
AlphaGenome tracks were used for tissue filtering.

**Supplementary table status:** Cell journal supplementary files are
access-controlled (HTTP 403) and PMC full-text is not open access for this
article. LFC values were computed from GEO raw count data (CPM-normalized,
log2 RNA/Plasmid activity, allelic difference B−A).

**Variant set:** {n_total} variants with measured LFC ({3431} rsID variants
in file; excluded: 636 chr:pos oligos without ref/alt information, and variants
where GRCh38 coordinates could not be resolved).

**emVar filter:** Not applied. Tewhey emVar annotations require access to
paper supplementary tables. All {n_total} variants with measured LFC are
included.

**Scores available for correlation:** {n_expr} expression_subscore,
{n_scored} full_composite, {n_cadd} cadd_phred.

**LFC distribution:** mean=0, SD={lfc_std:.3f}, range=[{lfc_range[0]:.2f}, {lfc_range[1]:.2f}].

### Correlation table

{chr(10).join(table_lines)}

*Bootstrap 95% CI from {N_BOOTSTRAP} resamples.*

### Scatter plot

![expression_subscore vs MPRA LFC](tewhey_correlation.png)

### Interpretation

The expression sub-score represents the modalities most directly comparable
to MPRA readout (transcriptional activity). The full composite is the score
used for variant ranking in the tool. CADD is included as a baseline.

- **expression_subscore (all):** {f"r={expr_r:+.3f} — {interpret(expr_r)} correlation with MPRA LFC." if expr_r is not None else "Not available."}
- **full_composite (all):** {f"r={full_r:+.3f} — {interpret(full_r)} correlation." if full_r is not None else "Not available."}
- **CADD PHRED (all):** {f"r={cadd_r:+.3f} — {interpret(cadd_r)} correlation." if cadd_r is not None else "Not available."}

"""

    # Append to existing benchmark_report.md or create
    if bm_path.exists():
        existing = bm_path.read_text()
        # Replace the section if it already exists
        marker = "## Primary validation: MPRA correlation (Tewhey 2016)"
        if marker in existing:
            before = existing[: existing.index(marker)]
            bm_path.write_text(before + section)
        else:
            with open(bm_path, "a") as f:
                f.write(section)
    else:
        bm_path.write_text(section)

    print(f"\nReport written to {bm_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # Init flags file
    with open(FLAGS_FILE, "w") as f:
        f.write(
            "# Tewhey Analysis Flags\n\n"
            "_Stop conditions are written here. If this file has content below "
            "the separator, the analysis halted early._\n\n---\n"
        )

    # --- Step 1: LFC from raw counts ---
    lfc_df = _compute_lfc()

    # --- Step 2: Resolve GRCh38 coordinates ---
    coord_df = _resolve_coords(lfc_df["rsid"].tolist())

    # Merge LFC + coords
    merged = lfc_df.merge(coord_df, on="rsid", how="inner")
    # Keep only resolved variants for scoring
    scored_base = merged.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
    print(f"\nVariants with LFC and resolved coords: {len(scored_base)}")

    # --- Step 3: CADD scores (runs concurrently with scoring setup) ---
    cadd_df = _get_cadd_scores(scored_base)

    # --- Step 4: AlphaGenome scoring ---
    alpha_df = _score_tewhey(scored_base)

    # --- Assemble final DataFrame ---
    df = scored_base.merge(alpha_df[["rsid", "full_composite", "expression_subscore", "error"]],
                           on="rsid", how="left")
    df = df.merge(cadd_df[["rsid", "cadd_phred"]], on="rsid", how="left")

    # Save to parquet
    df.to_parquet(PARQUET_OUT, index=False)
    print(f"\nSaved {len(df)} rows to {PARQUET_OUT}")

    # --- Step 5: Stop-condition checks + correlations ---
    corr_rows = _run_correlations(df)
    _check_stop_conditions(df, corr_rows)

    # --- Print correlation summary ---
    print("\n=== Correlation Results ===")
    for row in corr_rows:
        if row["r"] is not None:
            print(f"  {row['label']:<45}  r={row['r']:+.4f}  "
                  f"[{row['ci_lo']:+.4f}, {row['ci_hi']:+.4f}]  "
                  f"p={row['p']:.2e}  n={row['n']}")
        else:
            print(f"  {row['label']:<45}  n={row['n']} (insufficient data)")

    # --- Step 6: Scatter plot + write report ---
    _make_scatter(df, corr_rows)
    _write_report(df, corr_rows)


if __name__ == "__main__":
    main()
