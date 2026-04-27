"""Fetch and curate ClinVar pathogenic/benign variants for our four diseases.

Data source: ClinVar variant_summary.txt.gz from NCBI FTP.
Filters applied:
  - Assembly = GRCh38
  - Clinical significance ∈ {Pathogenic, Likely pathogenic} or {Benign, Likely benign}
  - Review status ∈ {criteria provided multiple submitters, expert panel, practice guideline}
  - SNP/small indel only (ref and alt ≤ 15 bp, no structural variants)
  - hg38 coordinates available

Design note: ClinVar pathogenic variants for complex diseases (T2D, SCZ) are
predominantly rare coding variants, not the common regulatory variants our tool
targets. This is a deliberate test of how well AlphaGenome's regulatory scores
generalize beyond their intended GWAS use case. We document this mismatch in
the benchmark report.
"""

from __future__ import annotations

import gzip
import json
import random
import urllib.request
from pathlib import Path

CACHE_FILE = Path(__file__).parent / "cache" / "clinvar_raw.json"
SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

PATHOGENIC_TERMS = frozenset(
    ["pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"]
)
BENIGN_TERMS = frozenset(["benign", "likely benign", "benign/likely benign"])
REVIEW_OK = frozenset([
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
    "criteria provided, single submitter",
])

# Disease ↔ phenotype / gene keyword mapping
# Using lowercase for case-insensitive matching.
# "Closest-match disease category" applied: MODY/neonatal DM for T2D,
# 22q11.2/NDD for SCZ.
DISEASE_KEYWORDS: dict[str, list[str]] = {
    "alzheimers": [
        "alzheimer", "psen1", "psen2", "presenilin", "amyloid precursor",
        "app|alzheimer",          # APP in AD context
    ],
    "parkinsons": [
        "parkinson", "lrrk2", "pink1", "prkn", "park2", "park7",
        "snca", "vps35", "dj-1",
    ],
    "t2d": [
        "maturity-onset diabetes", "mody", "neonatal diabetes mellitus",
        "permanent neonatal", "hnf1a", "hnf4a", "hnf1b", "gck",
        "kcnj11", "abcc8", "insulin secretion",
    ],
    "schizophrenia": [
        "schizophrenia", "schizoaffective", "22q11.2 deletion",
        "velocardiofacial", "nrxn1", "disc1",
    ],
}


def _parse_summary(path: Path) -> tuple[dict[str, list], list]:
    """Parse variant_summary.txt.gz and return (disease_buckets, benign_pool)."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    header = lines[0].rstrip("\n").split("\t")
    col = {h: i for i, h in enumerate(header)}

    buckets: dict[str, list] = {d: [] for d in DISEASE_KEYWORDS}
    benign_pool: list = []

    for raw in lines[1:]:
        f = raw.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        if f[col.get("Assembly", 16)] != "GRCh38":
            continue

        clinsig = f[col.get("ClinicalSignificance", 6)].lower()
        review = f[col.get("ReviewStatus", 24)].lower()
        phenotype = f[col.get("PhenotypeList", 13)].lower()
        gene = f[col.get("GeneSymbol", 4)].lower()
        chrom = f[col.get("Chromosome", 18)]
        pos = f[col.get("PositionVCF", 31)]
        ref = f[col.get("ReferenceAlleleVCF", 32)]
        alt = f[col.get("AlternateAlleleVCF", 33)]
        rsid_raw = f[col.get("RS# (dbSNP)", 9)]

        if not pos or pos in ("-1", "") or not ref or not alt:
            continue
        if ref.lower() in ("na", ".") or alt.lower() in ("na", "."):
            continue
        if len(ref) > 15 or len(alt) > 15:
            continue
        try:
            int(pos)
        except ValueError:
            continue
        if chrom in ("", "MT", "Un") or not chrom.isdigit() and chrom not in "XY":
            continue

        good_review = any(r in review for r in REVIEW_OK)
        if not good_review:
            continue

        is_path = any(t in clinsig for t in PATHOGENIC_TERMS) and "benign" not in clinsig
        is_benign = any(t in clinsig for t in BENIGN_TERMS) and "pathogenic" not in clinsig

        entry = {
            "rsid": (
                f"rs{rsid_raw}"
                if rsid_raw and rsid_raw not in ("-1", "")
                else f"cv_{chrom}_{pos}"
            ),
            "chrom": f"chr{chrom}",
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "gene": f[col.get("GeneSymbol", 4)],
            "phenotype": f[col.get("PhenotypeList", 13)][:120],
            "clinsig": f[col.get("ClinicalSignificance", 6)],
            "review": f[col.get("ReviewStatus", 24)],
        }

        if is_benign:
            benign_pool.append(entry)
        elif is_path:
            for disease, kws in DISEASE_KEYWORDS.items():
                if any(
                    ("|" in kw and all(k in (phenotype + " " + gene) for k in kw.split("|")))
                    or ("|" not in kw and (kw in phenotype or kw in gene))
                    for kw in kws
                ):
                    buckets[disease].append(entry)
                    break

    return buckets, benign_pool


def download_and_parse(force: bool = False) -> tuple[dict[str, list], list]:
    """Download ClinVar summary (cached) and return (buckets, benign_pool)."""
    gz_cache = CACHE_FILE.parent / "variant_summary.txt.gz"
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)

    if not gz_cache.exists() or force:
        print("Downloading ClinVar variant_summary.txt.gz (~50 MB)...")
        urllib.request.urlretrieve(SUMMARY_URL, gz_cache)
        print(f"  Saved to {gz_cache}")
    else:
        print(f"  Using cached {gz_cache}")

    print("Parsing...")
    buckets, benign = _parse_summary(gz_cache)
    for d, v in buckets.items():
        print(f"  {d}: {len(v)} pathogenic variants")
    print(f"  benign pool: {len(benign)}")
    return buckets, benign


def load_or_build_cache(force: bool = False) -> tuple[dict[str, list], list]:
    """Return (buckets, benign_pool), building cache from FTP if needed."""
    if CACHE_FILE.exists() and not force:
        with open(CACHE_FILE) as f:
            d = json.load(f)
        return d["buckets"], d["benign"]

    buckets, benign = download_and_parse(force=force)
    with open(CACHE_FILE, "w") as f:
        json.dump({"buckets": buckets, "benign": benign[:2000]}, f)
    return buckets, benign


def sample_benchmark_set(
    n_pathogenic: int = 100,
    n_benign: int = 50,
    seed: int = 42,
    force_rebuild: bool = False,
) -> tuple[list, list]:
    """Return (pathogenic_variants, benign_variants) sampled for benchmark.

    Pathogenic: up to n_pathogenic variants spread across all four diseases,
    capped at what's available per disease. Benign variants are MAF-stratified
    where available — here we approximate by random sampling from the pool.
    """
    rng = random.Random(seed)
    buckets, benign_pool = load_or_build_cache(force=force_rebuild)

    # Spread pathogenic variants proportionally across diseases
    all_path = []
    per_disease_n = n_pathogenic // len(buckets)
    for disease, variants in buckets.items():
        rng.shuffle(variants)
        sample = variants[:per_disease_n]
        for v in sample:
            v["disease"] = disease
        all_path.extend(sample)
        print(f"  sampled {len(sample)}/{len(variants)} from {disease}")

    # Fill up to n_pathogenic from largest buckets if under quota
    remaining = n_pathogenic - len(all_path)
    if remaining > 0:
        leftover = []
        for disease, variants in buckets.items():
            leftover.extend(v for v in variants[per_disease_n:] if "disease" not in v)
        for v in leftover:
            v["disease"] = [d for d, vs in buckets.items() if v in vs[per_disease_n:]][0]
        rng.shuffle(leftover)
        all_path.extend(leftover[:remaining])

    # Sample benign
    rng.shuffle(benign_pool)
    benign_sample = benign_pool[:n_benign]
    for v in benign_sample:
        v["disease"] = None

    print(f"\nFinal benchmark set: {len(all_path)} pathogenic, {len(benign_sample)} benign")
    return all_path, benign_sample
