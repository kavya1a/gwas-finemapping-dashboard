# Quantile normalization saturates regulatory variant scoring

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)

**Kavya Amrutham** · Barnard College, Columbia University · ka3041@barnard.edu

*Independent undergraduate research project. Not peer reviewed. Code, data, and methodology open for inspection.*

---

This repository documents a methodological investigation into how quantile normalization affects validation of regulatory variant prediction tools, using AlphaGenome as the test case.

Tools like AlphaGenome rank variants by their predicted regulatory impact relative to millions of common variants in the genome. That ranking is useful — but it has a blind spot. When you apply it to a set of variants already selected for regulatory relevance (GWAS hits, MPRA panels, eQTL credible sets), nearly all of them land at the extreme end of the scale. The predictor goes effectively binary. Correlation with experimental measurements drops to near zero — not because the model is wrong, but because the normalization step erases the differences between variants that the experiment is designed to detect.

Switching to raw model outputs (un-normalized deltas) recovers a directional signal that quantile normalization discards: Spearman correlation with MPRA measurements rises from ρ = 0.036 to ρ = 0.174, with no changes to the underlying model.

---

![Saturation CDF](figures/saturation_cdf.png)

*Empirical CDF of |expression subscore| for two regulatory-enriched variant sets, plotted against the theoretical uniform distribution (grey dashed) that AlphaGenome's quantile calibration is designed to produce for random common variation. Both empirical curves pile up near 1 while the diagonal stays linear. Tewhey MPRA loci (blue, 94.9% above 0.9) and disease GWAS variants (red, 99.6% above 0.9) are entirely compressed into the same extreme bin, leaving no room to rank them against each other.*

---

## Why this happens

AlphaGenome's published scoring pipeline converts raw per-track effect sizes into genome-wide percentiles, calibrated against ~300,000 common variants from gnomAD and 1000 Genomes. The result is a score in [−1, +1] where 0.9 means "this variant has a larger predicted regulatory effect than 90% of common genetic variation." For prioritizing candidates at a GWAS locus — asking which variant in a credible set has the largest regulatory footprint — this is exactly the right tool.

The problem arises at validation time, through two compounding mechanisms. First, MPRA panels and eQTL datasets are not random: they are assembled from variants that already showed up in GWAS, passed regulatory annotation filters, or were selected specifically because they might be functional. By the same logic that makes them interesting to study, they are enriched for regulatory activity relative to background. Second, the expression_subscore is computed as the maximum signed percentile across the selected tissue tracks (RNA-seq, CAGE, PRO-cap in K562). Taking the max over many tracks inflates saturation even for non-regulatory variants — with 20–30 K562 tracks, the probability that the max percentile exceeds 0.9 approaches 90–95% for any variant. The two effects compound: selection bias pushes most variants toward extreme scores, and the max aggregation fills in the rest. The result: 94.9% of Tewhey MPRA variants score above 0.9 in absolute value. 99.6% of GWAS disease variants do. When everything scores near ±1, the predictor cannot distinguish a variant that modestly nudges transcription from one that dramatically shuts it down.

This is not a bug in AlphaGenome, and it probably affects any tool that uses genome-wide quantile normalization — Enformer, Sei, and others that score in similar frameworks would face the same issue applied to the same kinds of test sets.

---

## What we found

We scored 3,259 variants from the Tewhey 2016 MPRA panel (GSE75661), a set of regulatory loci drawn from GWAS studies across 7.5k individuals. The genome-wide quantile score gives Spearman ρ = +0.036 against MPRA log-fold-change — statistically marginal. Before concluding the model fails, we worked through four alternative explanations:

- Allele orientation mismatch (ref/alt flipped vs. MPRA A/B convention) — ruled out, 100% match via Ensembl
- Wrong LFC column — ruled out, mpra_lfc = B − A matches Tewhey 2016 directly
- Score saturation — confirmed, 94.9% of scores are above |0.9|
- Coordinate drift from hg19 to hg38 — ruled out, 5/5 spot-checked positions are exact

The diagnosis is saturation. Magnitude correlation actually survives: |score| vs |LFC| gives ρ = +0.108 across the full panel, and among the top-5% highest-effect variants the signed correlation rises to ρ = +0.271. The model is tracking the right biology. The normalization step just compresses everyone into the same bin.

Switching to raw per-track expression deltas — extracted before quantile normalization, filtered to K562/blood-lineage tracks (RNA-seq, CAGE, PRO-cap) — recovers the continuous signal:

| Predictor | Spearman ρ | 95% CI | p | n |
|---|---|---|---|---|
| Raw max signed expression delta | **+0.174** | [+0.086, +0.255] | 1.8×10⁻⁵ | 600 |
| Raw mean signed expression delta | +0.148 | [+0.068, +0.223] | 2.7×10⁻⁴ | 600 |
| Quantile-normalized expression_subscore | +0.036 | [−0.002, +0.075] | 0.039 | 3,259 |
| \|Raw max delta\| vs \|mpra_lfc\| | +0.207 | [+0.124, +0.282] | 3.0×10⁻⁷ | 600 |

![Distribution comparison](figures/phase2_distribution.png)

The distribution comparison makes the issue concrete. Normalized scores pile up near ±1 — bimodal, no room to rank. Raw deltas form a smooth continuous distribution centered near zero. Most variants have small effects; a few have large ones. That's what you need to correlate against an MPRA.

![Correlation by LFC bin](figures/phase2_lfc_bins.png)

Breaking down by effect size shows where each predictor succeeds. Among variants with large measured effects (|LFC| > 0.5), raw deltas reach ρ = +0.614; normalized scores get to +0.245. In the low-effect bins where most variants sit, normalized scores are near zero or negative, while raw deltas stay modestly positive. The raw signal is most valuable precisely where it's hardest to separate signal from noise.

| \|LFC\| bin | Normalized ρ | n | Raw ρ | n |
|---|---|---|---|---|---|
| 0 – 0.05 | −0.015 | 1,512 | +0.092 | 281 |
| 0.05 – 0.10 | −0.000 | 892 | +0.099 | 168 |
| 0.10 – 0.20 | +0.044 | 553 | +0.383 | 93 |
| 0.20 – 0.50 | +0.185 | 242 | +0.031 | 43 |
| > 0.50 | +0.245 | 60 | **+0.614** | 15 |

We also confirmed that saturation isn't specific to Tewhey's MPRA design or to neurological and metabolic diseases. Platelet count GWAS variants (n=198) saturate at 99.0%; hemoglobin GWAS variants (n=195) at 100%. Same pattern, different biology, different experimental paradigm.

| Variant set | n | Expression >0.9 |
|---|---|---|
| Tewhey MPRA | 3,259 | 94.9% |
| Disease GWAS — expression | 767 | 99.6% |
| Disease GWAS — chromatin | 767 | 69.4% |
| Disease GWAS — TF binding | 767 | 31.9% |
| Platelet count GWAS | 198 | 99.0% |
| Hemoglobin GWAS | 195 | 100.0% |

The saturation gradient across modalities is mechanistic. Expression scores aggregate across three track types (RNA-seq, CAGE, PRO-cap), so the max over more tracks pushes more variants to the tail. TF binding, using only one track type, saturates the least at 32%.

![Blood trait replication](figures/phase3_saturation_cdf.png)

---

## What this means in practice

The normalized and raw scores answer different questions. Using either one for the wrong purpose gives a misleading answer.

For variant prioritization within a credible set — asking which of five fine-mapped candidates has the largest regulatory footprint — normalized scores are exactly right. APOE illustrates this well: rs429358 (ε4) scores −0.9996 and rs7412 (ε2) scores +0.9976, opposite signs that match their opposing Alzheimer's risk effects. That relative ordering is meaningful, and it's what the normalization is built to produce. The problem is when you take those same scores and try to rank thousands of variants across different loci against each other. Once everything compresses to ±1, you've lost the information about which variants are actually different from each other.

Raw deltas are what you need when the goal is correlation with an experimental measurement. An MPRA, a CRISPRi screen, a set of eQTL effect sizes — all of these are continuous, and you need a continuous predictor to compare against them. The raw delta is AlphaGenome's actual predicted effect size before it gets converted to a percentile. Most variants score near zero; a few score large. That spread is what makes correlation possible.

The practical consequences extend in a few directions. There's an active debate in the field about how well deep learning regulatory models actually generalize — papers benchmarking AlphaGenome, Enformer, Sei, and others against MPRA or eQTL data consistently report correlations in the 0.1–0.3 range, and the field treats these as a ceiling on what the models can do. But if the benchmark is saturating the normalization, the ceiling is artificial. The models may be capturing more real biology than the benchmarks suggest, and raw outputs would give a fairer comparison.

For rare disease interpretation, the distinction matters more directly. A de novo regulatory variant in a patient isn't interesting because it ranks in the 99th percentile of common variation — it's interesting if it's predicted to substantially disrupt expression of a dosage-sensitive gene. A raw delta near zero means the model thinks the substitution changes expression little. A large value is grounds for prioritizing the variant for functional follow-up. The percentile score would call both of them extreme and give you no way to tell them apart.

There's also a design implication for anyone planning an MPRA. If you select candidates from GWAS credible sets and use normalized scores to decide which variants to test, you'll find that most of them score similarly — the normalization doesn't discriminate between them. Raw deltas will separate them, which means you can focus experimental capacity on variants with predicted effects large enough to detect against the noise floor of the assay. This pipeline scores 363 variant-disease pairs across Alzheimer's, type 2 diabetes, schizophrenia, and Parkinson's disease; the per-locus rankings work as a first-pass filter for exactly this kind of prioritization.

---

## Setup

You need Python 3.11+ and an AlphaGenome API key ([request access here](https://alphagenome.google)). The scoring runs are included in the repo as cached databases, so you can reproduce all figures without re-running the API.

```bash
git clone https://github.com/kavya1a/regulatory-score-saturation.git
cd regulatory-score-saturation
pip install -r requirements.txt
cp .env.example .env
# add ALPHAGENOME_API_KEY to .env
```

To regenerate all figures from the cached data (takes about a minute, no API key needed):

```bash
make figures
```

To run the full pipeline from scratch — fetching variants, scoring them, and extracting raw deltas — expect about 20 hours of API calls:

```bash
make pipeline
make figures
make verify
```

`make verify` runs canonical variant tests. Three of five pass. The two failures are documented in [`docs/OVERNIGHT_BLOCKERS.md`](docs/OVERNIGHT_BLOCKERS.md).

---

## Repo structure

```
├── README.md
├── TEWHEY_RESULT.md        # full diagnostic writeup: four hypotheses, four correlation numbers
├── config.yaml             # tissue profiles, canonical variants, pipeline parameters
│
├── batch_score.py          # score GWAS variants → scored_variants.db
├── prefetch_variants.py    # fetch GWAS Catalog variants → preloaded_variants.db
├── tewhey_analysis.py      # download and score Tewhey MPRA panel
├── extract_raw_deltas.py   # extract raw expression deltas for 600 Tewhey variants
├── saturation_figure.py    # saturation CDF figures
├── phase2_figures.py       # distribution and LFC-bin correlation figures
├── phase3_blood_traits.py  # blood trait replication
├── gwas_catalog.py         # GWAS Catalog v2 REST client
├── allele_resolver.py      # Ensembl allele resolution with SQLite caching
│
├── scoring/                # AlphaGenome scoring utilities
├── verification/           # canonical variant test harness
├── figures/                # all generated figures
├── docs/                   # project log and operational notes
└── archive/                # pre-pivot scripts, kept for reference
```

All scoring results are included so you can inspect the data or regenerate figures without API access:

| File | Contents |
|---|---|
| `scored_variants.db` | 767 scored GWAS variant-disease pairs across 4 diseases |
| `tewhey_mpra.parquet` | Tewhey 2016 panel with AlphaGenome scores (3,259 variants) |
| `tewhey_raw_delta_cache.db` | raw expression deltas for 600 stratified Tewhey variants |
| `phase3_blood_cache.db` | expression subscores for 393 blood trait GWAS variants |
| `variants.db` | Ensembl allele resolution cache |

Data sources: Tewhey 2016 MPRA (GSE75661, GRCh38 liftover) · GWAS Catalog v2 · Ensembl REST API · AlphaGenome v0.6.1

---

## Limitations

**Sample size.** The raw delta correlation (ρ = +0.174, n=600) uses a stratified subsample; the normalized comparison uses all 3,259 scored variants. Confidence intervals are wider for the raw estimate. Stratification by |LFC| quintile ensures representation across effect sizes, but the estimate would tighten with a larger sample.

**Tissue mismatch.** Raw deltas use K562/blood-lineage tissue profiles. Tewhey 2016 measured regulatory activity in LCLs (lymphoblastoid cell lines), which have different chromatin accessibility and transcription factor occupancy. This mismatch suppresses correlation — ρ = +0.174 is a lower bound on what a correctly matched tissue profile would give.

**Single model, single normalization scheme.** This analysis tests one tool (AlphaGenome v0.6.1) with one normalization approach. The saturation pattern likely applies to other quantile-normalized regulatory predictors, but we don't test that claim here.

**Raw output access required.** Recovering the directional signal requires access to per-track model outputs before normalization. Some tool deployments only expose the normalized score. This limits applicability in contexts where API access is restricted to the published scoring pipeline.

**Quantile normalization is a design choice, not a flaw.** The normalization makes AlphaGenome scores interpretable as genome-wide percentiles and is appropriate for within-locus variant prioritization. The issue documented here is about using those scores for a purpose they weren't designed for — cross-locus validation against continuous measurements.

**Blood trait replication is partial.** Saturation on platelet count and hemoglobin GWAS variants is confirmed using quantile scores. Whether raw deltas on those variants also improve correlation is not tested.

**Max-over-tracks aggregation inflates saturation.** The expression_subscore is the maximum signed percentile across the selected K562 expression tracks (RNA-seq, CAGE, PRO-cap). For a variant set of n tracks, the max of n independent uniform random variables has P(max > 0.9) = 1 − 0.9^n. With 20–30 K562 tracks, this approaches 90–95% for any variant, regardless of regulatory function. We attempted a random-variant negative control using 1,000 common variants from gene-poor regions; they also saturated at ~97%, confirming that the aggregation step contributes to saturation alongside the regulatory enrichment of the primary datasets. The core finding — that raw deltas better predict MPRA correlation — is unaffected by this confound, since raw deltas avoid the quantile aggregation entirely.

**Canonical tests: 3/5 pass.** The two failures are documented in [`docs/OVERNIGHT_BLOCKERS.md`](docs/OVERNIGHT_BLOCKERS.md).

---

## Citation

```bibtex
@misc{amrutham2026saturation,
  author = {Amrutham, Kavya},
  title  = {Quantile normalization saturates regulatory variant scoring},
  year   = {2026},
  url    = {https://github.com/kavya1a/regulatory-score-saturation}
}
```

---

*AlphaGenome model and API by Google DeepMind.*
