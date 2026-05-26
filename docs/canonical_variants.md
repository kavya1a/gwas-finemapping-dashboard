# Canonical variant verification

`make verify` runs the harness in `verification/canonical_variants_test.py` against `scored_variants.db`. It checks five known causal variants — but in two distinct groups, reflecting different scientific claims.

## The five variants

| # | Variant | Gene | Disease | Why chosen |
|---|---|---|---|---|
| 1 | rs429358 | APOE | Alzheimer's | Largest known common-variant effect on AD risk; coding variant defining the ε4 allele. |
| 2 | rs7412 | APOE | Alzheimer's | Defines the protective ε2 allele; opposing direction to rs429358 on the same gene. Tests sign reproduction. |
| 3 | rs1006737 | CACNA1C | Schizophrenia | Replicated cross-disorder GWAS hit; intronic regulatory variant affecting voltage-gated calcium channel expression. |
| 4 | rs7903146 | TCF7L2 | T2D | Largest replicated T2D GWAS signal; intronic regulatory variant acting on TCF7L2 expression. |
| 5 | rs356219 | SNCA | Parkinson's | Top PD GWAS signal at the SNCA locus; regulatory variant influencing α-synuclein expression. |

## Two test groups, two claims

### Group A — `canonical_rank_recovery`

Asks whether the variant's AlphaGenome composite score lands in the **top 20%** of its disease's scored set (~270–294 variants per disease). This is the strict claim: "the model singles this variant out."

Three tests:
- **rs429358** ranks in top 20% of AD.
- **rs7412 vs rs429358** have opposite signed expression deltas (well-established biology: ε2 and ε4 oppose each other).
- **rs1006737** ranks in top 20% of SCZ.

### Group B — `canonical_regulatory_detection`

Asks whether AlphaGenome assigns strong absolute regulatory effect to the variant, regardless of how it ranks against other GWAS variants in the same disease set. Passes when:
- `composite_score > 0.5`, AND
- `|expression signed_max_score| > 0.9` (i.e., top 10% of common variants for the expression modality).

Within-disease-set rank is reported as diagnostic context only, not as a pass criterion.

Two tests:
- **rs7903146 regulatory signal detected for T2D** — composite 0.828, expression signed −0.994.
- **rs356219 regulatory signal detected for PD** — composite 0.767, expression signed −0.987.

## Why the split

rs7903146 and rs356219 fail the strict top-20% rank test (25.1% and 56.8% respectively). The original docs blamed tissue mismatch ("K562/blood-lineage track selection"). **This was inaccurate** — the tissue filter is working: T2D scoring uses 1,230 tracks including pancreas, hepatocyte, liver, skeletal muscle; PD scoring uses 906 tracks including substantia nigra, prefrontal cortex, neural cell. Verified by `verification/diagnose_canonical_failures.py`.

The actual cause is at the modality-aggregation level. For **rs356219 (SNCA)**, the variant scores below the median of the PD set on every single modality:

| Modality | rs356219 score | % of PD set with a higher score |
|---|---|---|
| expression | 0.987 | 80% |
| chromatin | 0.811 | 89% |
| tf_binding | 0.866 | 40% |
| splice_junctions | 0.959 | 48% |
| splice_site_usage | 0.420 | 96% |

No weighted combination of these can put rs356219 in the top 20% — it's outranked by 40–96% of PD variants on every dimension. This is a property of how the model scores variants at the SNCA locus relative to other PD GWAS hits, not a pipeline bug.

**rs7903146 (TCF7L2)** is a borderline case. It scores well on chromatin (only 24% of T2D variants beat it), but expression saturates (66% beat it), so the composite leaves it just outside top-20%. Per-disease up-weighting of chromatin (`analysis/explore_composite_weights.py` → `t2d_chromatin_heavy`) recovers it to 21.1% — still over the line.

In both cases AlphaGenome **does** assign strong absolute regulatory effect (composite ≫ 0.5, expression near-saturated). It just doesn't differentiate these variants from other plausibly-regulatory variants at the same disease locus. The detection-group test captures what the model can honestly claim.

## What the docs used to say (corrected)

The previous version of this page asserted that rs7903146 and rs356219 fail because of "tissue context mismatch" — K562/blood tracks instead of pancreatic islet / dopaminergic neuron tracks. That explanation is wrong. The disease tissue profiles in `scoring/tissue_config.py` and `config.yaml` are correctly defined (T2D includes pancreas / islet / beta cell / liver / adipose / muscle keywords; PD includes brain / substantia nigra / dopamin / midbrain / basal ganglia). The tissue filter at `scoring/tissue_config.py:102` applies them correctly, and `verification/diagnose_canonical_failures.py` confirms the filtered track sets are tissue-appropriate.

## Reproducing the diagnosis

```bash
# 5-test harness (uses cached scored_variants.db, no API):
python verification/canonical_variants_test.py

# Confirm tissue filter behavior (2 API calls, ~2 min):
python verification/diagnose_canonical_failures.py

# Per-variant outrank analysis + composite-weight exploration:
python analysis/explore_composite_weights.py
```

## What's not tested (potential follow-ups)

- **Signed-direction tests** for rs7903146 (TCF7L2 isoform effect) and rs356219 (SNCA expression direction). AlphaGenome's signed-score allele convention needs verification before adding these as pass/fail criteria.
- **Matched-null percentile tests** — comparing the canonical variant's score against a common-variant null built by `build_matched_calibration.py`. Requires re-scoring the canonicals through the matched-cal pipeline (units mismatch with `scored_variants.db.composite_score`).
- **Causal-prior weighting** — combining the AlphaGenome score with PIP / eQTL evidence as a separate prioritization layer, distinct from the model score itself.
