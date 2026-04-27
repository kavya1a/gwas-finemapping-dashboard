# Quantile normalization saturates regulatory variant scoring

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)

**Kavya Amrutham** · Barnard College, Columbia University · ka3041@barnard.edu

---

Tools like AlphaGenome rank variants by their predicted regulatory impact relative to millions of common variants in the genome. That ranking is useful — but it has a blind spot. When you apply it to a set of variants already selected for regulatory relevance (GWAS hits, MPRA panels, eQTL credible sets), nearly all of them land at the extreme end of the scale. The predictor goes effectively binary. Correlation with experimental measurements drops to near zero — not because the model is wrong, but because the normalization step erases the differences between variants that the experiment is designed to detect.

Switching to raw model outputs (un-normalized deltas) brings Spearman correlation with MPRA measurements from ρ = 0.036 to ρ = 0.174 — a fivefold improvement, with no changes to the underlying model.

---

![Saturation CDF across three modalities](figures/saturation_cdf.png)

*Empirical CDF of |expression subscore| for three variant sets. A random draw of common variants would follow the diagonal by construction. 94.9% of Tewhey MPRA variants and 99.6% of disease GWAS variants score above 0.9 — the predictor has nowhere left to go.*

---

## Why this happens

AlphaGenome's published scoring pipeline converts raw per-track effect sizes into genome-wide percentiles, calibrated against ~300,000 common variants from gnomAD and 1000 Genomes. The result is a score in [−1, +1] where 0.9 means "this variant has a larger predicted regulatory effect than 90% of common genetic variation." For prioritizing candidates at a GWAS locus — asking which variant in a credible set has the largest regulatory footprint — this is exactly the right tool.

The problem arises at validation time. MPRA panels and eQTL datasets are not random. They are assembled from variants that already showed up in GWAS, passed regulatory annotation filters, or were selected specifically because they might be functional. By the same logic that makes them interesting to study, they are enriched for regulatory activity relative to background. The calibration treats them all as extreme, because they are. 94.9% of Tewhey MPRA variants score above 0.9 in absolute value. 99.6% of GWAS disease variants do. When everything scores near ±1, the predictor cannot distinguish a variant that modestly nudges transcription from one that dramatically shuts it down.

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

The normalized and raw scores answer different questions, and using either for the wrong one gives misleading results.

Normalized scores are the right choice when you're working within a credible set at a single GWAS locus — asking which of five fine-mapped variants has the largest regulatory footprint relative to the rest of the genome. APOE is a good example: rs429358 (ε4) scores −0.9996 and rs7412 (ε2) scores +0.9976, opposite signs consistent with their opposing Alzheimer's risk. That directional information is meaningful within a locus. The scores fail when you try to use them to rank thousands of variants across different loci against each other, because compression to ±1 removes the differences between them.

Raw model outputs are the right choice when you're comparing against experimental measurements. If you're running an MPRA, doing a CRISPRi screen, or benchmarking against eQTL effect sizes, the raw delta gives you the absolute predicted effect of the substitution — how much AlphaGenome thinks the variant changes expression in relevant tissues. That continuous scale is what enables meaningful correlation.

A few concrete applications:

**Benchmarking regulatory tools.** Most published comparisons of tools like AlphaGenome, Enformer, and Sei against MPRA or eQTL data report correlations in the 0.1–0.3 range and interpret them as reflecting model quality. Some of that signal is real, but part of it is a normalization artifact — the tools are being evaluated on test sets that saturate their calibration. Using raw outputs would give a more accurate picture of which models are actually learning regulatory logic.

**Rare variant interpretation.** For a de novo variant in a patient, you often want to know whether the variant is likely to have a large regulatory effect in absolute terms — not whether it ranks in the 99th percentile of common variation. Raw deltas give you that: a value near zero means the substitution is predicted to change expression little; a large value means it probably matters.

**MPRA panel design.** Before running an MPRA, you might want to prioritize which candidate variants to include based on predicted effect size. Normalized scores will rank most GWAS-derived candidates similarly (they all score near 1). Raw deltas will actually separate them, letting you focus experimental effort on variants with predicted effects large enough to detect.

**Fine-mapping support.** Within a credible set, normalized scores work well and are what the pipeline is built for. The pipeline here scores 363 GWAS variant-disease pairs across Alzheimer's, type 2 diabetes, schizophrenia, and Parkinson's disease, and reports per-locus rankings that can feed directly into fine-mapping workflows alongside PIPs from SuSiE or FINEMAP.

---

## Setup

You need Python 3.11+ and an AlphaGenome API key ([request access here](https://alphagenome.google)). The scoring runs are included in the repo as cached databases, so you can reproduce all figures without re-running the API.

```bash
git clone https://github.com/kavya1a/gwas-finemapping-dashboard.git
cd gwas-finemapping-dashboard
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

`make verify` runs canonical variant tests. Three of five pass; the two failures are documented in [`docs/OVERNIGHT_BLOCKERS.md`](docs/OVERNIGHT_BLOCKERS.md) and reflect tissue-profile thresholds, not incorrect scores.

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

The raw delta correlation (ρ = +0.174, n=600) uses a stratified subsample while the normalized comparison uses all 3,259 scored variants, so the confidence intervals are wider for the raw estimate. The stratification is by |LFC| quintile to ensure representation across effect sizes, but a larger sample would tighten the estimate.

Raw deltas here use K562/blood-lineage tissue profiles. Tewhey 2016 assayed LCL regulatory activity, which is related but not identical. A better-matched tissue profile would likely improve the correlation further — ρ = +0.174 is a lower bound.

Saturation for blood trait GWAS variants (Phase 3) was confirmed using quantile scores only. We didn't extract raw deltas for those variants, so we can't confirm that the correlation improvement generalizes beyond the Tewhey dataset. The saturation finding generalizes; the raw delta fix is demonstrated only on Tewhey.

The canonical variant test suite passes 3/5. The two failures are sensitivity issues with the tissue profile configuration, not wrong scores — but they should be resolved before using this pipeline to make clinical claims.

---

## Citation

```bibtex
@misc{amrutham2026saturation,
  author = {Amrutham, Kavya},
  title  = {Quantile normalization saturates regulatory variant scoring},
  year   = {2026},
  url    = {https://github.com/kavya1a/gwas-finemapping-dashboard}
}
```

---

*AlphaGenome model and API by Google DeepMind.*
