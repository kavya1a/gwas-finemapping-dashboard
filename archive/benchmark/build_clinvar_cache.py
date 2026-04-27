"""Build ClinVar benchmark cache from UCSC genome browser API.

Fetches ClinVar variants for target gene regions using the UCSC REST API,
which returns pre-parsed data with vcfDesc field containing ref/alt alleles.
Much faster than downloading the full ClinVar summary TSV.

Run with:
    /opt/homebrew/bin/python3.11 benchmark/build_clinvar_cache.py
"""

import collections
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

UCSC = (
    "https://api.genome.ucsc.edu/getData/track"
    "?genome=hg38&track=clinvarMain&chrom={chrom}&start={start}&end={end}"
)

PATHOGENIC_TERMS = frozenset(["pathogenic", "likely pathogenic"])
BENIGN_TERMS = frozenset(["benign", "likely benign"])

# Gene regions (GRCh38).  Coordinates from Ensembl GRCh38.
DISEASE_REGIONS: dict[str, list[tuple]] = {
    "alzheimers": [
        ("PSEN1", "chr14", 73136528, 73223691),
        ("PSEN2", "chr1",  226868000, 226903000),
        ("APP",   "chr21",  25880550, 26170969),
    ],
    "parkinsons": [
        ("LRRK2", "chr12", 40224986, 40369285),
        ("PRKN",  "chr6",  161768442, 163148834),
        ("PINK1", "chr1",   20959769,  20978274),
        ("SNCA",  "chr4",   89724098,  89838315),
    ],
    "t2d": [
        ("HNF1A",  "chr12", 120978166, 121002094),
        ("GCK",    "chr7",   44184960,  44236418),
        ("KCNJ11", "chr11",  17407720,  17415057),
        ("ABCC8",  "chr11",  17415059,  17498404),
        ("HNF4A",  "chr20",  43039866,  43107570),
        ("HNF1B",  "chr17",  36046354,  36105069),
    ],
    "schizophrenia": [
        ("NRXN1", "chr2",  49919008,  51259674),
        ("TBX1",  "chr22",  19744225,  19766343),   # 22q11.2 gene
        ("TBR1",  "chr2", 162264006, 162287024),
        ("SHANK3","chr22",  51105436,  51237541),
    ],
}

# Benign controls from well-characterized genes
BENIGN_REGIONS: list[tuple] = [
    ("PSEN1", "chr14", 73136528, 73223691),   # same genes as pathogenic → matched context
    ("LRRK2", "chr12", 40224986, 40369285),
    ("GCK",   "chr7",  44184960,  44236418),
    ("HNF1A", "chr12", 120978166, 121002094),
]

VCF_RE = re.compile(r"^chr\w+:(\d+):([ACGT-]+)>([ACGT-]+)$")


def fetch_region(chrom: str, start: int, end: int) -> list[dict]:
    url = UCSC.format(chrom=chrom, start=start, end=end)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read())
        return data.get("clinvarMain", [])
    except Exception as e:
        print(f"    ERR fetching {chrom}:{start}-{end}: {e}")
        return []


def parse_item(item: dict, gene: str) -> dict | None:
    """Extract a usable variant entry from a UCSC clinvarMain item."""
    vdesc = item.get("vcfDesc", "")
    m = VCF_RE.match(vdesc)
    if not m:
        return None

    pos_raw, ref, alt = m.group(1), m.group(2), m.group(3)
    if len(ref) > 15 or len(alt) > 15 or ref == "-" or alt == "-":
        return None

    chrom = item.get("chrom", "")
    clinsig = item.get("clinSign", "").strip()
    snp_id = item.get("snpId", "").strip()
    star_count = int(item.get("_starCount", 0) or 0)
    pheno = item.get("phenotypeList", "")[:100]

    return {
        "rsid": snp_id if snp_id else f"cv_{chrom}_{pos_raw}",
        "chrom": chrom,
        "pos": int(pos_raw),   # vcfDesc pos is already 1-based
        "ref": ref,
        "alt": alt,
        "gene": gene,
        "phenotype": pheno,
        "clinsig": clinsig,
        "star_count": star_count,
    }


def classify(item: dict) -> str | None:
    """Return 'pathogenic', 'benign', or None."""
    sig = item.get("clinsig", "").lower()
    # Exclude conflicting, VUS, etc.
    if "conflict" in sig or "uncertain" in sig or "not provided" in sig:
        return None
    if any(t in sig for t in PATHOGENIC_TERMS) and "benign" not in sig:
        return "pathogenic"
    if any(t in sig for t in BENIGN_TERMS) and "pathogenic" not in sig:
        return "benign"
    return None


def collect_disease_variants(min_stars: int = 1) -> tuple[dict[str, list], list]:
    buckets: dict[str, list] = collections.defaultdict(list)
    seen: set[str] = set()

    print("Collecting pathogenic variants...")
    for disease, regions in DISEASE_REGIONS.items():
        for gene, chrom, start, end in regions:
            print(f"  {disease}/{gene} {chrom}:{start}-{end}", end="  ")
            items = fetch_region(chrom, start, end)
            added = 0
            for item in items:
                entry = parse_item(item, gene)
                if not entry:
                    continue
                if entry["star_count"] < min_stars:
                    continue
                uid = f"{entry['chrom']}_{entry['pos']}_{entry['ref']}_{entry['alt']}"
                if uid in seen:
                    continue
                label = classify(entry)
                if label == "pathogenic":
                    entry["disease"] = disease
                    buckets[disease].append(entry)
                    seen.add(uid)
                    added += 1
            print(f"→ {added} pathogenic")
            time.sleep(0.2)

    print("\nCollecting benign controls...")
    benign: list = []
    benign_seen: set[str] = set()
    for gene, chrom, start, end in BENIGN_REGIONS:
        print(f"  {gene} {chrom}:{start}-{end}", end="  ")
        items = fetch_region(chrom, start, end)
        added = 0
        for item in items:
            entry = parse_item(item, gene)
            if not entry:
                continue
            if entry["star_count"] < min_stars:
                continue
            uid = f"{entry['chrom']}_{entry['pos']}_{entry['ref']}_{entry['alt']}"
            if uid in benign_seen:
                continue
            label = classify(entry)
            if label == "benign":
                entry["disease"] = None
                benign.append(entry)
                benign_seen.add(uid)
                added += 1
        print(f"→ {added} benign")
        time.sleep(0.2)

    return dict(buckets), benign


def main():
    print("=" * 60)
    print("Building ClinVar benchmark cache via UCSC API")
    print("=" * 60)

    out = Path(__file__).parent / "cache" / "clinvar_raw.json"
    out.parent.mkdir(parents=True, exist_ok=True)

    buckets, benign = collect_disease_variants(min_stars=1)

    print("\n=== Summary ===")
    for d, vs in buckets.items():
        print(f"  {d}: {len(vs)} pathogenic")
        for v in vs[:2]:
            print(f"    {v['rsid']} {v['chrom']}:{v['pos']} "
                  f"{v['ref']}>{v['alt']} gene={v['gene']} "
                  f"★{v['star_count']} {v['clinsig']}")
    print(f"  benign pool: {len(benign)}")

    with open(out, "w") as f:
        json.dump({"buckets": buckets, "benign": benign}, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
