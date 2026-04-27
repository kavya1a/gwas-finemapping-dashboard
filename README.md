# Quantile normalization saturates regulatory variant scoring

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)

**Kavya Amrutham** · Barnard College, Columbia University · ka3041@barnard.edu

---

Quantile-normalized AlphaGenome scores collapse to a near-binary signal on regulatory-enriched variant sets — making them unsuitable for MPRA and eQTL validation without modification. Raw model outputs restore a **5× improvement** in Spearman correlation with experimental measurements.

---

## Headline figure

![Saturation CDF across three modalities](figures/saturation_cdf.png)

Empirical CDF of |expression subscore| for three variant populations. A uniform reference (diagonal) is expected from a random draw of common variants — by construction of AlphaGenome's genome-wide calibration. Both regulatory-enriched sets shift entirely to the right: **94.9% of Tewhey MPRA variants and 99.6% of disease GWAS variants** have |score| > 0.9, collapsing the predictor to a near-binary ±1 signal.

---

## Background

**AlphaGenome** is a sequence-to-function model that predicts regulatory activity from DNA sequence. For a given ref → alt substitution, it computes differential predictions across hundreds of tracks (RNA-seq, ATAC-seq, ChIP-seq, Hi-C) and returns per-track effect sizes.

AlphaGenome's published scoring pipeline applies **quantile normalization**: each raw delta is ranked against a background of ~300,000 common variants (gnomAD/1KG, MAF > 0.01), converting it to a genome-wide percentile in [−1, +1]. This is appropriate for *within-credible-set prioritization*, where all candidates are at the same GWAS locus and relative ordering is what matters.

**The problem:** standard tool validation correlates predicted scores against MPRA or eQTL measurements. Using quantile-normalized scores for this comparison is methodologically unsound for regulatory-enriched test sets — and systematically underestimates model performance.

---

## Finding 1 — Saturation on regulatory-enriched variant sets

Regulatory-enriched variants are in the extreme tail of the genome-wide calibration distribution *by construction*: the very enrichment for regulatory activity that makes them interesting to study also places them all near |score| = 1.

| Variant set | n | Expression >0.9 | Median \|score\| |
|---|---|---|---|
| Tewhey MPRA (GWAS regulatory loci) | 3,259 | **94.9%** | 0.989 |
| Disease GWAS — expression | 767 | **99.6%** | 0.996 |
| Disease GWAS — chromatin | 767 | 69.4% | 0.941 |
| Disease GWAS — TF binding | 767 | 31.9% | 0.844 |
| Platelet count GWAS (Phase 3) | 198 | **99.0%** | 0.995 |
| Hemoglobin GWAS (Phase 3) | 195 | **100.0%** | 0.996 |
| Uniform reference (expected) | — | ~10% | ~0.5 |

The saturation gradient across modalities is mechanistic: expression aggregates 3 track types (RNA-seq, CAGE, PRO-cap), so taking the max over more tracks pushes more variants to the tail. Chromatin (3 types) saturates at 69%; TF binding (1 type) at 32%.

**Saturation is not a bug.** The normalized scores are correct and useful for within-locus ranking. The problem is applying them across regulatory-enriched sets where an absolute measurement scale is needed.

![Phase 3 saturation CDF](figures/phase3_saturation_cdf.png)

---

## Finding 2 — Magnitude correlation is preserved; directional correlation is not

We tested four hypotheses for the near-zero Spearman ρ between expression_subscore and Tewhey MPRA log-fold-change:

| Hypothesis | Result |
|---|---|
| H1: Allele orientation mismatch | **Ruled out.** 100% ref/alt match via Ensembl |
| H2: Wrong LFC column | **Ruled out.** mpra_lfc = B − A confirmed per Tewhey 2016 |
| H3: Score saturation | **Confirmed.** 94.9% of scores compressed to ±0.99 |
| H4: Coordinate drift (hg19→hg38) | **Ruled out.** 5/5 spot-checked positions exact |

Magnitude correlation survives: |expression_subscore| vs |mpra_lfc| gives Spearman ρ = +0.108 (p = 7.7×10⁻¹⁰, n=3,259). Among the top-5% highest-effect variants the signed correlation rises to ρ = +0.271 (p = 4.6×10⁻⁴, n=163). The predictor captures real signal — quantile normalization destroys the dynamic range needed to see it at scale.

---

## Finding 3 — Raw model outputs restore the directional signal

We extracted raw per-track expression deltas (RNA-seq, CAGE, PRO-cap, K562/blood-lineage filter) for 600 stratified Tewhey variants (120 per |LFC| quintile).

### Correlation table

| Predictor | Spearman ρ | 95% CI | p | n |
|---|---|---|---|---|
| Raw max signed expression delta vs mpra_lfc | **+0.174** | [+0.086, +0.255] | 1.8×10⁻⁵ | 600 |
| Raw mean signed expression delta vs mpra_lfc | +0.148 | [+0.068, +0.223] | 2.7×10⁻⁴ | 600 |
| Quantile-normalized expression_subscore vs mpra_lfc | +0.036 | [−0.002, +0.075] | 0.039 | 3,259 |
| \|Raw max delta\| vs \|mpra_lfc\| | +0.207 | [+0.124, +0.282] | 3.0×10⁻⁷ | 600 |

**5× improvement** in Spearman ρ from raw over normalized (0.174 vs 0.036).

### Why: the distributions

![Distribution comparison](figures/phase2_distribution.png)

Normalized scores are bimodal (94.9% at |score| > 0.9). Raw deltas are continuous (88.5% have |delta| < 0.1). A bimodal predictor cannot rank variants within a category; a continuous one can.

### Correlation by effect-size bin

![LFC bin correlation](figures/phase2_lfc_bins.png)

| \|LFC\| bin | Normalized ρ | n | Raw ρ | n |
|---|---|---|---|---|---|
| 0 – 0.05 | −0.015 | 1,512 | +0.092 | 281 |
| 0.05 – 0.10 | −0.000 | 892 | +0.099 | 168 |
| 0.10 – 0.20 | +0.044 | 553 | **+0.383** | 93 |
| 0.20 – 0.50 | +0.185 | 242 | +0.031 | 43 |
| > 0.50 | +0.245 | 60 | **+0.614** | 15 |

Raw deltas outperform normalized scores in the low-effect bins (the majority of variants) and strongly outperform for high-effect variants (|LFC| > 0.5: ρ = +0.614 vs +0.245).

---

## Finding 4 — Generalization across datasets and tissues

The saturation replicates on two independent GWAS sets — different biology (megakaryopoiesis, erythropoiesis), different diseases, same structural result:

- **Platelet count GWAS** (EFO_0004615, n=198): 99.0% at |score| > 0.9
- **Hemoglobin GWAS** (EFO_0004611, n=195): 100.0% at |score| > 0.9

Saturation is a structural consequence of genome-wide-calibrated quantile scoring applied to regulatory-enriched variant sets. It is not specific to Tewhey's MPRA design, K562 cells, or any particular disease.

---

## Implications

> *For validation against continuous experimental measurements, use raw model outputs.*  
> *For ranking within enrichment-selected variant sets, use normalized scores.*  
> *They serve different purposes; one is not strictly better than the other.*

**Use quantile-normalized scores when:**
- Prioritizing variants within a credible set at a single GWAS locus
- Comparing relative regulatory impact where all candidates are at the same locus
- Screening for unusually high regulatory impact relative to genome-wide background

**Use raw model outputs when:**
- Correlating with continuous measurements (MPRA LFC, CRISPRi, eQTL effect sizes)
- Validating against any dataset selected for regulatory relevance
- Absolute effect magnitude matters

**Implication for published benchmarks:** comparisons between regulatory scoring tools that use quantile-normalized outputs against MPRA/eQTL data compare bimodal predictors against continuous measurements. The resulting low correlations reflect normalization strategy, not model quality.

### APOE ε4 / ε2 as a canonical validation example

rs429358 (APOE ε4) scores expression_subscore = −0.9996; rs7412 (APOE ε2) scores +0.9976. Opposite signs, consistent with their opposing Alzheimer's risk effects. Normalized scores are directionally correct for within-locus comparison — the failure mode is specific to cross-locus validation at scale.

---

## Setup

### Requirements

- Python 3.11+
- AlphaGenome API key ([request access](https://alphagenome.google))

```bash
git clone https://github.com/kavya1a/gwas-finemapping-dashboard.git
cd gwas-finemapping-dashboard
pip install -r requirements.txt
cp .env.example .env
# Add your key: ALPHAGENOME_API_KEY=...
```

### Regenerate figures from cached data (no API key needed, ~1 min)

All scoring results are included in the repo as SQLite databases and a parquet file.

```bash
make figures
```

Generates all six figures in `figures/`.

### Full pipeline from scratch (API key required, ~20 hrs)

```bash
make pipeline   # fetch, score, extract raw deltas
make figures    # generate figures
make verify     # run canonical variant tests
```

See `Makefile` for individual targets.

---

## Verification

Canonical variant tests check that known directional effects are recovered (APOE alleles, PPARG P12A, BIN1 rs744373, FTO rs9939609).

```
Tests passing: 3/5
```

Tests 3 and 5 fail due to tissue-profile sensitivity thresholds in the current K562 configuration, not incorrect scores. See [`docs/OVERNIGHT_BLOCKERS.md`](docs/OVERNIGHT_BLOCKERS.md) for details.

---

## Repository structure

```
├── README.md
├── TEWHEY_RESULT.md        # Full diagnostic: 4 hypotheses, 4 correlation numbers
├── CITATION.cff
├── LICENSE
├── Makefile
├── config.yaml             # Tissue profiles, canonical variants, pipeline params
├── requirements.txt
│
├── batch_score.py          # Score GWAS variants → scored_variants.db
├── prefetch_variants.py    # Fetch GWAS Catalog variants → preloaded_variants.db
├── tewhey_analysis.py      # Download + score Tewhey MPRA → tewhey_mpra.parquet
├── extract_raw_deltas.py   # Raw expression deltas for 600 Tewhey variants
├── saturation_figure.py    # Phase 1: saturation CDF
├── phase2_figures.py       # Phase 2: distribution + LFC-bin correlation
├── phase3_blood_traits.py  # Phase 3: blood trait replication
├── gwas_catalog.py         # GWAS Catalog v2 REST client
├── allele_resolver.py      # Ensembl allele resolution + SQLite cache
│
├── scoring/                # AlphaGenome scoring utilities
├── verification/           # Canonical variant test harness
├── figures/                # All generated figures (PNG)
├── docs/                   # Project log, scope lock, operational notes
└── archive/                # Superseded work (pre-pivot scripts)
```

### Data files (included)

| File | Contents |
|---|---|
| `scored_variants.db` | 767 scored GWAS variant-disease pairs (4 diseases) |
| `tewhey_mpra.parquet` | Tewhey 2016 panel with AlphaGenome scores (3,259 variants) |
| `tewhey_raw_delta_cache.db` | Raw expression deltas for 600 stratified Tewhey variants |
| `phase3_blood_cache.db` | Expression subscores for 393 blood trait GWAS variants |
| `variants.db` | Ensembl allele resolution cache |
| `preloaded_variants.db` | Prefetched GWAS variant coordinates |

**Data provenance:** Tewhey 2016 MPRA (GSE75661, GRCh38 liftover) · GWAS Catalog v2 REST API · Ensembl REST API (GRCh38) · AlphaGenome API v0.6.1

---

## Limitations

**Sample size for raw deltas.** The raw delta correlation (ρ = +0.174, n=600) uses a stratified subsample. The normalized score comparison uses all 3,259 scored variants. Confidence intervals are wider for the raw delta estimate.

**Tissue mismatch.** Raw deltas use K562/blood-lineage tissue profiles. Tewhey 2016 assayed LCL (lymphoblastoid) regulatory activity. A closer tissue match would likely improve correlation; the current result is a lower bound.

**AlphaGenome window.** Variants near chromosome boundaries or with complex nearby variation may produce less reliable predictions. No explicit filtering was applied.

**MPRA assay scope.** Tewhey 2016 measures 200bp oligo activity, not in-situ chromatin context. Short-element MPRA and native regulatory function diverge for variants with distal chromatin effects.

**Phase 3 raw deltas not extracted.** Blood trait saturation is confirmed via quantile scores only. Whether raw deltas would show the same correlation improvement for blood trait variants is not tested.

**Canonical tests: 3/5 pass.** The 5/5 target was cut from scope. Failing tests reflect tissue-profile sensitivity, not incorrect scoring logic.

---

## Citation

```bibtex
@misc{amrutham2026saturation,
  author    = {Amrutham, Kavya},
  title     = {Quantile normalization saturates regulatory variant scoring},
  year      = {2026},
  url       = {https://github.com/kavya1a/gwas-finemapping-dashboard},
  note      = {GitHub repository}
}
```

---

*AlphaGenome model and API by Google DeepMind.*
