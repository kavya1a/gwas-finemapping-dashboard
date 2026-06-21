# Matched-statistic calibration

This document specifies the procedure used to re-derive AlphaGenome's quantile
calibration so that the calibration's summary statistic matches the one used
at test time.

## Background

AlphaGenome's published `quantile_score` is computed by ranking each variant's
**single-track** raw score against an empirical background of ~300,000 common
variants (gnomAD / 1000 Genomes, MAF > 0.01). The calibration is correct for
its design purpose: per-track, within-locus ranking.

In a tissue-level pipeline a different summary statistic is typically applied
on top — e.g., the maximum signed value across all K562/blood-lineage
expression tracks. Composing that order statistic with single-track-calibrated
quantiles is what produces saturation: a large fraction of *random* common
variants gets pushed against the ±1 ceiling because the max-of-many-quantiles
sits in the tail of the single-track null even when no underlying variant is
unusual.

The fix is methodological: re-derive the quantile calibration with the **same**
summary statistic that will be applied at test time.

## Procedure

The procedure is parameterized in `config.yaml` under `matched_calibration:` and
implemented in `build_matched_calibration.py`. It is reproducible from
`seed = 2026` alone.

### 1. Sample windows

* 22 autosomes (chr1-22; sex chromosomes excluded).
* `n_windows = 66` total, allocated proportional to chromosome length.
* Each window is `window_bp = 50,000` bp wide.
* Window starts are drawn uniformly at random within
  `[min_telomere_bp, length - min_telomere_bp - window_bp]` per chromosome
  (i.e., at least 1 Mb from each chromosome end), seed = 2026.

This yields ~3.3 Mb of genome sampled, allocated:
chr1 → 7 windows, chr2 → 6, chr3..chr7 → 4-5 each, chr8..chr13 → 3 each,
chr14..chr18 → 2 each, chr19..chr22 → 1 each.

### 2. Variant fetch and filter

* Within each window, all variants from the region were filtered to
  MAF > 0.01 via gnomAD v3 GraphQL (the same allele-frequency source
  AlphaGenome calibrates against), keeping biallelic SNVs only.
* Drop any variant whose rsID appears in `tewhey_mpra.parquet`. Three
  variants were dropped in the present run.

A deterministic seed-based subset capped at **6,000** candidates
(`config.yaml` → `matched_calibration.pool_cap`) is sent to the scoring
stage so that overall runtime is bounded; the subset is reproducible from
`seed` and `pool_cap`.

### 3. Score with the matched summary statistic

Each variant is scored through AlphaGenome v0.6.1 with the K562/blood-lineage
expression pipeline used elsewhere in this repo (see
`extract_raw_deltas.py` and `scoring/composite.py`). For each variant the
score retained is

```
raw_max_signed_delta = sign(argmax(|raw_score|)) * max(|raw_score|)
```

over all RNA_SEQ + CAGE + PROCAP tracks that pass the K562/blood tissue
filter from `config.yaml`. The peak track's published `quantile_score` is
also retained as `single_track_quantile` for the side-by-side saturation
comparison below.

The scoring is parallelized 4-way over variants. Each individual call wraps
the existing 60-second per-variant timeout from `batch_score.py`. Suspected
rate-limit responses (`"rate"`, `"quota"`, `"limit"`, `"429"`) abort the run.
Writes to `matched_calibration_cache.db` go through a `threading.Lock` and
are committed per row, so the run is fully resumable on interruption.

### 4. Empirical null

After scoring, all rows with non-null `raw_max_signed_delta` are written to
`matched_calibration_null.parquet` and form the empirical null distribution
used at test time.

In the committed run (`matched_calibration_cache.db`): **5,993 variants
scored, 60 API timeouts dropped, 0 non-timeout errors → 5,933 clean rows in
`matched_calibration_null.parquet`.**

### 5. Apply at test time

For each test variant (e.g., a Tewhey panel variant) with summary statistic
`x = raw_max_signed_delta` against the null array `N`:

```
n_lt = #{n in N : n < x}
n_eq = #{n in N : n == x}
q_unsigned = (n_lt + 0.5 * n_eq) / |N|     in (0, 1)
matched_quantile = 2 * q_unsigned - 1       in (-1, +1)
```

This mirrors AlphaGenome's signed-quantile convention. A monotone-transform
phred form is also computed for downstream use:

```
phred_empirical = -10 * log10(1 - matched_quantile + 1e-6)
```

`phred_empirical` is a strict monotone function of `matched_quantile` and
therefore preserves rank correlation against any third variable. It is
reported only as a diagnostic and a possible interface for downstream
software that prefers PHRED-style scaling. Note that AlphaGenome SDK v0.6.1
does not expose a phred-scaled output (verified across the source tree, all
ten tagged releases, the issue tracker, and the Nature paper supplement).
This `phred_empirical` is derived locally from the matched null and is not
identical to any internal phred convention DeepMind may use.

## The smoking-gun comparison

On the **same** 5,933 random common autosomal variants, two scores are
computed:

* `raw_max_signed_delta` — the matched (max-over-tracks) raw delta. Its
  empirical distribution is the new calibration null.
* `single_track_quantile` — the published AlphaGenome quantile of the
  *same* peak track (single-track calibrated, max-over-tracks applied).

| | Matched null (raw max-over-tracks Δ) | Published single-track quantile (peak track) |
|---|---|---|
| Saturation `|score| > 0.9` | **0.42 %** | **41.4 %** |
| Exact ceiling `≈ ±1` | 0.000 % | 1.1 % |
| Shape | Sharp spike at 0, heavy tails (kurtosis 81.8) | U-shape pushed against ±1 |
| Mean / median / std | +0.004 / +0.009 / 0.131 | (saturated; descriptive stats not informative) |
| Range | [-2.135, +2.637] | [-1, +1] (closed by construction) |

The histogram is in `figures/matched_calibration_histogram.png`.

The ~100× gap on identical variants is direct evidence that saturation in
the published-quantile pipeline is produced by the calibration-statistic
mismatch, not by any property of the variants being scored.

## Effect on Tewhey

Applied to the full Tewhey 2016 panel (n = 3,246 with all four predictors
valid):

| Predictor | Spearman ρ vs MPRA LFC | 95% CI | p | n |
|---|---|---|---|---|
| Original quantile (single-track calibration, max applied) | +0.0367 | [-0.0025, +0.0722] | 3.6e-2 | 3,246 |
| Matched-calibration quantile | +0.1225 | [+0.0876, +0.1573] | 2.5e-12 | 3,246 |
| Phred empirical | +0.1225 | [+0.0876, +0.1573] | 2.5e-12 | 3,246 |
| Raw max signed delta (no normalization) | +0.1225 | [+0.0876, +0.1573] | 2.6e-12 | 3,246 |

Saturation snapshot on the Tewhey panel:

| | `|·| > 0.9` |
|---|---|
| Original quantile | 94.9 % |
| Matched quantile | 12.9 % |
| Raw max signed Δ | 1.1 % |

`matched_quantile`, `phred_empirical`, and `raw max signed delta` produce
identical Spearman to numerical precision (|Δρ| = 0 to 8 decimals) because
they are monotone transforms of each other on the relevant Tewhey range.
The differences between them are scale and distribution shape, not rank.
The published quantile is a different rank ordering — that's where
correlation is lost.

## Files

| File | Source / contents |
|---|---|
| `build_matched_calibration.py` | Variant sampling, scoring, null parquet, histogram. |
| `analyze_matched_calibration.py` | Apply the null to Tewhey, four-row Spearman + figure. |
| `matched_calibration_cache.db` | Per-variant raw delta + single-track quantile (resumable). |
| `matched_calibration_null.parquet` | Clean null distribution used for re-quantiling. |
| `matched_calibration_comparison.csv` | Four-row Spearman table from `analyze_matched_calibration.py`. |
| `figures/matched_calibration_histogram.png` | Side-by-side null comparison (smoking gun). |
| `figures/three_way_comparison.png` | Tewhey distribution under each predictor + phred diagnostic. |
| `config.yaml` (`matched_calibration:` block) | All sampling and runtime parameters. |
