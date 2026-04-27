# Draft findings — AlphaGenome GWAS scoring pipeline

_Status: draft for review. Sections marked [VERIFY] require additional data or
re-scoring runs before inclusion in the preprint. All statistics cited below are
drawn directly from `scored_variants.db` and `preloaded_variants.db`; no values
are imputed or estimated._

---

## Finding 1: Cross-disease recovery of canonical causal variants

**Status: partially verified. See [VERIFY] notes before including verbatim.**

We scored the top-ranked GWAS variants for four diseases — Alzheimer's disease
(AD, 95 variants), type 2 diabetes (T2D, 76 variants), schizophrenia (SCZ, 99
variants), and Parkinson's disease (PD, 93 variants) — totalling 363 unique
variant-disease pairs. Variants were selected from the GWAS Catalog using
MONDO disease ontology identifiers and ranked by association p-value; allele
resolution used the Ensembl REST API with an SQLite cache.

The most consistently high-scoring variant across diseases was **rs429358
(APOE ε4, chr19:44908684 T→C; CADD PHRED = 16.6)**, which appeared in the AD
fetch (composite score 0.846, rank 12/95), the PD fetch (composite 0.824, rank
18/93), and the T2D fetch (composite 0.825, rank 21/76). Across all three
disease contexts this variant scored in the top decile of expression
(max\_abs\_score 0.9996–0.9997) and chromatin remodelling (0.9914–0.9962)
modalities, consistent with the established APOE promoter and enhancer activity
at this locus. The near-identical scores across different tissue-weight profiles
(AD: neural/microglia; PD: substantia nigra; T2D: pancreatic β-cell) reflect the
ubiquitous regulatory impact of this variant rather than tissue-specific scoring
artefacts.

> **[VERIFY — ranks]** The GWAS Catalog returns variants ranked by p-value; the
> specific ranks above (AD #12, PD #18, T2D #21) are computed from composite
> AlphaGenome scores within each disease cohort and may differ from GWAS
> significance ranking. The manuscript should clarify which ranking scheme is
> cited.

> **[VERIFY — TCF7L2, CACNA1C, SNCA]** rs7903146 (TCF7L2/T2D), rs1006737
> (CACNA1C/SCZ), and rs356219 (SNCA/PD) were **not recovered** in the GWAS
> Catalog top-100 returns for their respective diseases in the current pipeline
> run. These three variants require either: (a) fetching more than 100 variants
> per disease (increase `VARIANTS_PER_DISEASE` in `prefetch_variants.py`) or
> (b) confirming that our MONDO ID mapping returns them in a larger pull. Until
> these variants are scored, the "canonical causal variant recovery" claim
> cannot be made for T2D, SCZ, and PD with the current data. **Do not assert
> these rankings without running the extended fetch.**

---

## Finding 2: Direction-of-effect preservation at allelic series

**Status: methodologically justified; numeric verification table pending.**

The AlphaGenome scoring pipeline produces **signed** per-modality effect
scores (`signed_max_score` in `scoring/null_distribution.py` line 155) that
encode the direction of regulatory effect for the scored alternate allele. A
positive value indicates the alternate allele increases the relevant signal
(e.g., RNA-seq, chromatin accessibility); a negative value indicates a
decrease. This directional encoding is scientifically critical for
interpreting GWAS hits: two variants at the same locus can have opposite
allelic effects, and a tool that only reports unsigned magnitude cannot
distinguish between them.

The APOE locus provides the natural test case. The ε4 allele is defined by
the T→C substitution at rs429358 (chr19:44908684), which increases AD risk;
the ε2 allele is partly defined by rs7412 (chr19:44908822), which reduces
AD risk. If our signed scores correctly capture allele-specific directionality,
the expression and chromatin `signed_max_score` values at these two positions
should have opposite sign, reflecting the opposing regulatory effects of
risk-increasing vs. protective substitutions.

> **[VERIFY — comparison table]** `signed_max_score` values are computed in
> memory during scoring (via `scoring/null_distribution.py::summarize_modality_scores`)
> but are **not persisted** to `scored_variants.db` — only `max_abs_score` is
> stored. rs7412 was also not in the current GWAS fetch. To generate the
> direction-of-effect comparison table for this section:
>
> 1. Add a `signed_max_score` column to the `modality_scores` table in
>    `batch_score.py::_save_single_result` (extract from `modality_breakdown`
>    before discarding it).
> 2. Score rs429358 and rs7412 directly via `batch_score.py alzheimers` after
>    adding them to `preloaded_variants.db` if not already present.
> 3. Expected result: rs429358 T→C should show a positive expression
>    `signed_max_score` (ε4 increases APOE expression → increased AD risk via
>    amyloid pathway); rs7412 T→C at the protective ε2 haplotype should show
>    a negative or attenuated signed score. Any inversion between the two
>    supports direction-of-effect preservation.
>
> This verification requires one additional scoring run and a one-line DB
> schema addition. It is strongly recommended before making this claim in print.

---

## Pipeline statistics (for Methods section)

| Metric | Value |
|---|---|
| Diseases scored | 4 (AD, T2D, SCZ, PD) |
| Variants fetched (total) | 366 |
| Filtered at prefetch (indel) | 2 |
| API errors during scoring | 1 (rs9940149, telomere-proximal) |
| Variants scored successfully | 363 (99.2% yield) |
| Scoring framework | AlphaGenome (DeepMind, GRCh38), per-variant sequential with 60 s timeout |
| Allele resolution | Ensembl REST API POST, GRCh38, SQLite-cached |
| Disease ontology | MONDO (alzheimers: MONDO_0004975, t2d: MONDO_0005148, scz: MONDO_0005090, pd: MONDO_0005180) |
| CADD baseline | v1.7 GRCh38, local tabix query (whole\_genome\_SNVs.tsv.gz) |
