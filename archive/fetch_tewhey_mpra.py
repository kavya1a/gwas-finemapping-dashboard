"""Task 4: Download Tewhey 2016 MPRA data (GSE75661), compute activity scores,
resolve/liftover to GRCh38, and save to tewhey_mpra.parquet.

GEO only provides raw (unnormalized) barcode counts — no processed log2FC table
exists in the supplementary files. This script computes:
  activity_score = log2( mean(RNA_CPM across reps) / mean(plasmid_CPM across reps) )
  mpra_effect    = activity_score_B - activity_score_A  (alt minus ref allele)

Two datasets are processed:
  - 7.5k set: oligos named {rsid}_[RC_]{A|B} — GRCh38 coords via Ensembl
  - 79k set:  oligos named {chr}:{pos}:{type}_{A|B} (hg19) — liftover via pyliftover

Do NOT run AlphaGenome scoring on this set — that is tomorrow's work.
"""

from __future__ import annotations

import io
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BLOCKERS_FILE = Path(__file__).parent / "OVERNIGHT_BLOCKERS.md"
OUT_PARQUET = Path(__file__).parent / "tewhey_mpra.parquet"
CHAIN_FILE = Path(__file__).parent / "hg19ToHg38.over.chain.gz"

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75661/suppl"
FILES = {
    "7k5": f"{GEO_BASE}/GSE75661_7.5k_collapsed_counts.txt.gz",
    "79k": f"{GEO_BASE}/GSE75661_79k_collapsed_counts.txt.gz",
}
UCSC_CHAIN_URL = "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

PSEUDOCOUNT = 1.0  # added before log2 to avoid log(0)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_blocker(msg: str) -> None:
    with open(BLOCKERS_FILE, "a") as f:
        f.write(f"\n## BLOCKER [{time.strftime('%Y-%m-%d %H:%M:%S')}] fetch_tewhey_mpra\n{msg}\n")
    print(f"!! BLOCKER: {msg}")


def _download_gz(url: str) -> pd.DataFrame:
    print(f"  Downloading {url.split('/')[-1]}...")
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    return pd.read_csv(io.BytesIO(data), sep="\t", compression="gzip")


def _ensure_chain_file() -> Path:
    if CHAIN_FILE.exists():
        return CHAIN_FILE
    print(f"  Downloading hg19→hg38 chain file...")
    with urllib.request.urlopen(UCSC_CHAIN_URL, timeout=120) as resp:
        data = resp.read()
    CHAIN_FILE.write_bytes(data)
    print(f"  Chain file saved to {CHAIN_FILE}")
    return CHAIN_FILE


def _compute_activity(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-oligo activity scores from raw count table.

    Input columns: Oligo, Plasmid_r1..rN, <CellLine>_r1..rM
    Returns DataFrame with: oligo, activity_score (log2 RNA/plasmid CPM ratio)
    """
    cols = df.columns.tolist()
    plasmid_cols = [c for c in cols if c.lower().startswith("plasmid")]
    rna_cols = [c for c in cols if not c.lower().startswith("plasmid") and c != "Oligo"]

    if not plasmid_cols or not rna_cols:
        raise ValueError(f"Unexpected columns: {cols[:10]}")

    # Normalize to CPM within each replicate
    count_matrix = df.set_index("Oligo")
    for col in plasmid_cols + rna_cols:
        total = count_matrix[col].sum()
        count_matrix[col] = count_matrix[col] / total * 1e6

    plasmid_mean = count_matrix[plasmid_cols].mean(axis=1)
    rna_mean = count_matrix[rna_cols].mean(axis=1)

    activity = np.log2((rna_mean + PSEUDOCOUNT) / (plasmid_mean + PSEUDOCOUNT))

    result = pd.DataFrame({"oligo": count_matrix.index, "activity_score": activity.values})
    return result


def _parse_oligo_name_7k5(oligo: str) -> tuple[str | None, str | None]:
    """
    Parse 7.5k oligo name → (rsid, allele).
    Formats: rs12345_A, rs12345_B, rs12345_RC_A, rs12345_RC_B
    Returns (rsid, allele) where allele is 'A' or 'B'.
    """
    m = re.match(r"^(rs\d+)(?:_RC)?_([AB])$", oligo)
    if m:
        return m.group(1), m.group(2)
    return None, None


def _parse_oligo_name_79k(oligo: str) -> tuple[str | None, int | None, str | None]:
    """
    Parse 79k oligo name → (chrom_hg19, pos_hg19, allele).
    Format: chr10:109566984:D_A or chr10:110926269:I_B
    Returns (chrom, pos, allele) or (None, None, None).
    """
    m = re.match(r"^(chr[\w]+):(\d+):[A-Za-z]+_([AB])$", oligo)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None, None, None


# ---------------------------------------------------------------------------
# 7.5k dataset (rsID-based → Ensembl for GRCh38)
# ---------------------------------------------------------------------------

def _process_7k5() -> pd.DataFrame:
    print("\n--- 7.5k dataset ---")
    raw = _download_gz(FILES["7k5"])
    print(f"  {len(raw)} oligos, columns: {raw.columns.tolist()[:6]}...")

    activity = _compute_activity(raw)

    # Parse oligo names
    activity[["rsid", "allele"]] = activity["oligo"].apply(
        lambda o: pd.Series(_parse_oligo_name_7k5(o))
    )
    activity = activity.dropna(subset=["rsid", "allele"])

    # Pivot to wide: one row per rsid, columns activity_A and activity_B
    pivot = activity.pivot_table(index="rsid", columns="allele", values="activity_score", aggfunc="mean")
    pivot.columns = [f"activity_{c}" for c in pivot.columns]
    pivot = pivot.reset_index()

    if "activity_A" not in pivot.columns or "activity_B" not in pivot.columns:
        _write_blocker(
            "7.5k pivot missing activity_A or activity_B columns. "
            f"Columns found: {pivot.columns.tolist()}"
        )
        return pd.DataFrame()

    pivot["mpra_effect"] = pivot["activity_B"] - pivot["activity_A"]

    # Resolve GRCh38 coordinates from Ensembl
    from allele_resolver import resolve_alleles_batch, DB_PATH as ALLELE_DB
    rsids = pivot["rsid"].tolist()
    print(f"  Resolving {len(rsids)} rsIDs to GRCh38 via Ensembl...")
    resolved = resolve_alleles_batch(rsids, db_path=ALLELE_DB)

    pivot["chrom_hg38"] = pivot["rsid"].map(lambda r: (resolved.get(r) or {}).get("chrom"))
    pivot["pos_hg38"] = pivot["rsid"].map(lambda r: (resolved.get(r) or {}).get("pos"))
    pivot["ref"] = pivot["rsid"].map(lambda r: (resolved.get(r) or {}).get("ref"))
    pivot["alt"] = pivot["rsid"].map(lambda r: (resolved.get(r) or {}).get("alt"))
    pivot["maf"] = pivot["rsid"].map(lambda r: (resolved.get(r) or {}).get("maf"))
    pivot["multi_allelic"] = pivot["rsid"].map(
        lambda r: (resolved.get(r) or {}).get("multi_allelic", False)
    )

    pivot["dataset"] = "tewhey_7k5"
    n_resolved = pivot["chrom_hg38"].notna().sum()
    print(f"  {len(pivot)} variants; {n_resolved} resolved to GRCh38")
    return pivot


# ---------------------------------------------------------------------------
# 79k dataset (chr:pos hg19 → liftover → hg38)
# ---------------------------------------------------------------------------

def _process_79k() -> pd.DataFrame:
    print("\n--- 79k dataset ---")
    raw = _download_gz(FILES["79k"])
    print(f"  {len(raw)} oligos, columns: {raw.columns.tolist()[:6]}...")

    activity = _compute_activity(raw)

    records = []
    for _, row in activity.iterrows():
        chrom, pos, allele = _parse_oligo_name_79k(row["oligo"])
        if chrom and pos and allele:
            records.append({"chrom_hg19": chrom, "pos_hg19": pos,
                            "allele": allele, "activity_score": row["activity_score"]})
    if not records:
        _write_blocker("79k: could not parse any oligo names — check format")
        return pd.DataFrame()

    act_df = pd.DataFrame(records)

    pivot = act_df.pivot_table(
        index=["chrom_hg19", "pos_hg19"], columns="allele",
        values="activity_score", aggfunc="mean"
    )
    pivot.columns = [f"activity_{c}" for c in pivot.columns]
    pivot = pivot.reset_index()

    if "activity_A" not in pivot.columns or "activity_B" not in pivot.columns:
        _write_blocker(
            f"79k pivot missing activity_A/B. Columns: {pivot.columns.tolist()}"
        )
        return pd.DataFrame()

    pivot["mpra_effect"] = pivot["activity_B"] - pivot["activity_A"]

    # Liftover hg19 → hg38
    chain = _ensure_chain_file()
    try:
        from pyliftover import LiftOver
        lo = LiftOver(str(chain))
    except Exception as exc:
        _write_blocker(f"pyliftover failed to load chain file: {exc}")
        return pd.DataFrame()

    def _lift(row):
        result = lo.convert_coordinate(row["chrom_hg19"], row["pos_hg19"] - 1)  # 0-based
        if result:
            chrom38, pos38_0, strand, _ = result[0]
            return chrom38, pos38_0 + 1  # back to 1-based
        return None, None

    lifted = pivot.apply(_lift, axis=1, result_type="expand")
    pivot["chrom_hg38"] = lifted[0]
    pivot["pos_hg38"] = lifted[1]

    n_lifted = pivot["chrom_hg38"].notna().sum()
    print(f"  {len(pivot)} variants; {n_lifted} lifted to hg38")

    pivot["rsid"] = None  # no rsIDs in 79k oligo names
    pivot["ref"] = None
    pivot["alt"] = None
    pivot["maf"] = None
    pivot["multi_allelic"] = False
    pivot["dataset"] = "tewhey_79k"
    return pivot


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    dfs = []

    try:
        df_7k5 = _process_7k5()
        if not df_7k5.empty:
            dfs.append(df_7k5)
    except Exception as exc:
        _write_blocker(f"7.5k processing failed: {exc}")
        print(f"  7.5k FAILED: {exc}")

    try:
        df_79k = _process_79k()
        if not df_79k.empty:
            dfs.append(df_79k)
    except Exception as exc:
        _write_blocker(f"79k processing failed: {exc}")
        print(f"  79k FAILED: {exc}")

    if not dfs:
        _write_blocker("Both Tewhey datasets failed — tewhey_mpra.parquet not written.")
        return

    combined = pd.concat(dfs, ignore_index=True, sort=False)

    # Canonical column order
    ordered_cols = [
        "dataset", "rsid", "chrom_hg38", "pos_hg38",
        "chrom_hg19", "pos_hg19",
        "ref", "alt", "maf", "multi_allelic",
        "activity_A", "activity_B", "mpra_effect",
    ]
    combined = combined.reindex(columns=[c for c in ordered_cols if c in combined.columns])
    combined.to_parquet(OUT_PARQUET, index=False)
    print(f"\nSaved {len(combined)} rows to {OUT_PARQUET}")
    print(combined.describe())


if __name__ == "__main__":
    main()
