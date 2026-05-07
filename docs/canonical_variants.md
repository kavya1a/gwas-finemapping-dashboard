# Canonical variant verification

`make verify` runs five tests against `scored_variants.db` to check that the scoring pipeline recovers known causal variants at biologically expected ranks. Three pass, two fail. This page documents what each test asks, why those specific variants were chosen, and what the failures mean.

## The five variants

Each variant is a textbook causal hit for its disease — the kind of result a working regulatory predictor should rank highly. They span four diseases (Alzheimer's, T2D, schizophrenia, Parkinson's) and three causal mechanisms (coding, regulatory cis-eQTL, intronic enhancer).

| # | Variant | Gene | Disease | Why chosen |
|---|---|---|---|---|
| 1 | rs429358 | APOE | Alzheimer's | Largest known common-variant effect on AD risk; coding variant defining the ε4 allele. |
| 2 | rs7412 | APOE | Alzheimer's | Defines the protective ε2 allele; opposing direction to rs429358 on the same gene. Tests sign reproduction. |
| 3 | rs7903146 | TCF7L2 | T2D | Largest replicated T2D GWAS signal; intronic regulatory variant acting on TCF7L2 expression. |
| 4 | rs1006737 | CACNA1C | Schizophrenia | Replicated cross-disorder GWAS hit; intronic regulatory variant affecting voltage-gated calcium channel expression. |
| 5 | rs356219 | SNCA | Parkinson's | Top PD GWAS signal at the SNCA locus; regulatory variant influencing α-synuclein expression. |

Each rank-based test passes if the variant lands in the top 20% of its disease's scored set (300 variants per disease). Test 2 also requires that rs7412 and rs429358 have opposite signs on their expression subscore.

## Results

| # | Test | Result | Score / detail |
|---|---|---|---|
| 1 | rs429358 in AD top 20% | **PASS** | rank 40 / 294 (13.6%); composite score 0.843 |
| 2 | APOE ε4 / ε2 expression direction inversion | **PASS** | rs429358 expression signed = −0.9996; rs7412 = +0.9976 |
| 3 | rs7903146 in T2D top 20% | **FAIL** | rank 69 / 275 (25.1%); composite score 0.828 |
| 4 | rs1006737 in SCZ top 20% | **PASS** | rank 49 / 285 (17.2%); composite score 0.844 |
| 5 | rs356219 in PD top 20% | **FAIL** | rank 154 / 271 (56.8%); composite score 0.767 |

The two passing rank tests (rs429358, rs1006737) are both variants whose regulatory or coding effect is detectable in the K562/blood-lineage tracks the pipeline uses. The APOE sign test passes because the ε4 / ε2 substitutions are direct coding changes — well captured regardless of tissue context.

## Why the two rank tests fail

Both failures share a single mechanism: tissue context mismatch.

**Test 3 — rs7903146 (TCF7L2, rank 69/275 = 25.1%).** The variant is a near-miss: composite score 0.828 is high in absolute terms, but other K562-detected variants saturate ahead of it (consistent with the saturation pattern documented in the README). TCF7L2's causal regulatory effect for T2D is in pancreatic β-islet cells, where the variant disrupts an enhancer driving islet-specific expression. K562 tracks (chronic myelogenous leukemia, erythroid lineage) provide no signal for islet enhancers, so the model can only score the variant via generic regulatory features rather than its disease-relevant context.

**Test 5 — rs356219 (SNCA, rank 154/271 = 56.8%).** A more substantial miss: composite score 0.767 places this variant near the median of the PD set rather than the top. SNCA's regulatory effect operates in dopaminergic neurons of the substantia nigra, a cell type entirely absent from the K562/blood-lineage track selection. The pipeline has no way to detect a brain-specific cis-regulatory effect through hematopoietic tracks. The miss is the predicted consequence of the chosen tissue profile, not a model failure.

Both failures are consistent with the broader tissue-mismatch limitation noted in the README: scoring with K562 tracks works for variants whose mechanism manifests in hematopoietic cells (APOE coding effects, CACNA1C regulatory activity in K562) and fails for variants whose mechanism is restricted to tissues outside that track set (pancreatic islets, dopaminergic neurons). Re-running these two tests with matched-tissue track sets is the obvious next step and would either recover the expected ranks (confirming the diagnosis) or persist (refuting it).
