# Data dictionary

Every parquet / SQLite file committed to the repo, column by column. All artifacts cached so figures and analyses regenerate without API access.

---

## `matched_calibration_cache.db` → table `scores`

Per-variant K562 expression summary statistics for the 5,993 random common autosomal SNVs scored for the matched-calibration null. 5,933 rows have a clean (error-free) score; the other 60 are API timeouts (dropped from the null).

| Column | Type | Description |
|---|---|---|
| `rsid` | TEXT PK | dbSNP identifier; primary key |
| `chrom` | TEXT | GRCh38 chromosome (`chr1`–`chr22`) |
| `pos` | INTEGER | 1-based GRCh38 position |
| `ref` / `alt` | TEXT | Reference / alternate allele (single base) |
| `maf` | REAL | gnomAD v3 allele frequency (filter ≥ 0.01) |
| `raw_max_signed_delta` | REAL | `sign(argmax|raw|) · max|raw|` across all K562 expression tracks for this variant; the "peak track" signal AlphaGenome's pipeline would aggregate |
| `single_track_quantile` | REAL | `quantile_score` of the *same peak-track row* — i.e. AlphaGenome's published quantile for the track that produced the max raw delta |
| `raw_mean_signed_delta` | REAL | Arithmetic mean of signed `raw_score` across the K562 expression tracks (added 2026-05-21) |
| `raw_median_signed_delta` | REAL | Median of signed `raw_score` across the K562 expression tracks (added 2026-05-21) |
| `n_expr_tracks` | INTEGER | Number of K562-filtered expression tracks contributing to the aggregation (typically 300+) |
| `error` | TEXT | Failure mode if scoring did not complete (`api_timeout`, etc.); NULL on success |
| `scored_at` | INTEGER | Unix timestamp of the most recent scoring attempt |

`matched_calibration_null.parquet` is the clean subset of this table (filtered to `error IS NULL AND raw_max_signed_delta IS NOT NULL`), sorted ascending by `raw_max_signed_delta`. Used directly as the empirical CDF for re-quantiling Tewhey variants.

---

## `tewhey_raw_delta_cache.db` → table `raw_deltas`

Per-variant K562 expression summary statistics for the 3,301-variant Tewhey 2016 panel. 3,275 rows have a clean score.

| Column | Type | Description |
|---|---|---|
| `rsid` | TEXT PK | dbSNP identifier from the Tewhey 2016 panel |
| `max_signed_raw` | REAL | Max-over-tracks signed raw delta (same definition as in the matched cache) |
| `mean_signed_raw` | REAL | Mean-over-tracks signed raw delta |
| `n_expr_tracks` | INTEGER | Tracks contributing (median 339; min 16; max 3,178) |
| `error` | TEXT | NULL on success |
| `scored_at` | INTEGER | Unix timestamp |

Tewhey median is **not** cached here — only max and mean. To extract median over Tewhey, re-score via `extract_raw_deltas.py` with a small modification (~20 hours of API).

---

## `tewhey_mpra.parquet`

Tewhey 2016 MPRA panel joined with AlphaGenome composite scores. 3,301 rows.

| Column | Type | Description |
|---|---|---|
| `rsid` | str | dbSNP id |
| `activity_A`, `activity_B` | float | Raw MPRA activity for the two alleles (B = alternate) |
| `mpra_lfc` | float | `log2(activity_B / activity_A)` — the headline measurement we correlate against |
| `chrom`, `pos`, `ref`, `alt` | str / int | GRCh38 coordinates |
| `maf` | float | gnomAD MAF |
| `full_composite` | float | AlphaGenome multi-modality composite score (weights from `config.yaml`) |
| `expression_subscore` | float | **The published quantile-calibrated expression score with max-over-tracks aggregation applied.** This is the saturated predictor (94.9% above \|0.9\|). |
| `error` | str / NaN | Scoring failure flag |
| `cadd_phred` | float | CADD v1.7 phred score (independent reference predictor) |

---

## `matched_calibration_null.parquet`

Clean subset of `matched_calibration_cache.db.scores` ready for re-quantiling. Same columns; 5,933 rows; sorted by `raw_max_signed_delta` ascending so binary searches over the empirical CDF are O(log n).

---

## `scored_variants.db`

GWAS variant scoring across four diseases (Alzheimer's, T2D, schizophrenia, Parkinson's). 1,125 variant–disease pairs.

| Table | Key columns | Notes |
|---|---|---|
| `scores` | `(rsid, disease)` | composite_score, pip_weighted_score, error, rank |
| `modality_scores` | `(rsid, disease, modality)` | per-modality max/mean/signed scores; 767 pairs have `signed_max_score` cached |
| `run_metadata` | `(run_id, disease)` | provenance for each scoring run |

Used by the saturation-by-modality figures and by `make verify` for the canonical-variant tests.

---

## `phase3_blood_cache.db` → table `scores`

Blood-trait GWAS variant scoring (platelet count + hemoglobin panels). 393 variants.

| Column | Type | Description |
|---|---|---|
| `rsid` | TEXT | dbSNP id |
| `trait` | TEXT | Trait label |
| `chrom`, `pos`, `ref`, `alt` | TEXT / INT | GRCh38 coords |
| `expression_subscore` | REAL | Published quantile-calibrated expression score (max-aggregated) — same definition as Tewhey |
| `error` | TEXT | NULL on success |
| `scored_at` | INTEGER | Unix timestamp |

Used in the blood-trait replication of the saturation pattern (99.0–100% above \|0.9\|). Matched-calibration recipe has **not** been applied to this panel yet — known limitation.

---

## `variants.db` → table `allele_cache`

Ensembl REST API cache for allele resolution. Used by all scoring scripts to canonicalize ref/alt for an rsid before submission.

| Column | Type | Description |
|---|---|---|
| `rsid` | TEXT PK | dbSNP id |
| `chrom`, `pos`, `ref`, `alt`, `maf` | TEXT / INT / REAL | Resolved canonical record from Ensembl |
| `multi_allelic` | INTEGER (0/1) | Flag for variants with >2 alleles in dbSNP |
| `not_found` | INTEGER (0/1) | Flag for rsids that Ensembl didn't return |
| `fetched_at` | INTEGER | Unix timestamp of fetch |

---

## CSV outputs (not under version control discipline; regenerable)

| File | Produced by | Contents |
|---|---|---|
| `matched_calibration_comparison.csv` | `analyze_matched_calibration.py` | Four-predictor Spearman + CI on Tewhey (n = 3,246) |
| `matched_recipes_comparison.csv` | `analyze_matched_calibration_recipes.py` | Three-recipe (max / mean / median) saturation + Spearman comparison |
| `mean_aggregation_comparison.csv` | `analyze_mean_aggregation.py` | Tewhey raw-side max vs mean |

---

## Reproducibility notes

- Window-sampling RNG seed is **`2026`** (see `config.yaml` → `matched_calibration.seed`). The 66 windows, the 6,000-variant pool cap, and the gnomAD region fetches are all deterministic from this seed.
- AlphaGenome SDK version pinned to **v0.6.1**. Newer SDK versions may expose different track sets. The `phred_empirical` scale used here is derived locally from the matched null (see `docs/matched_calibration.md`); a grep across SDK source, all ten tagged releases, the issue tracker, the docs, and the Nature supplement found no phred-scaled output in v0.6.1.
- The K562 tissue profile (`config.yaml` → `tissue_profiles.tewhey_k562`) is a keyword filter; per-variant track counts vary 16 → 3,178 depending on which biosample tags AlphaGenome returns.
