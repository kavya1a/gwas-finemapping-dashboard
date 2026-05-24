.PHONY: install figures phase1 phase2 phase3 raw-stratified raw-full matched_calibration verify pipeline help

PYTHON := python3.11

help:
	@echo "Targets:"
	@echo "  install          Install dependencies"
	@echo "  figures          Regenerate all figures from cached data (fast)"
	@echo "  phase1           Saturation CDF figure"
	@echo "  phase2           Distribution + LFC-bin correlation figures"
	@echo "  phase3           Blood trait replication figure"
	@echo "  raw-stratified   Re-run stratified 600-variant raw delta extraction"
	@echo "  raw-full         Score full 3,259-variant Tewhey panel (requires API, ~8 hrs)"
	@echo "  matched_calibration  Build matched null + run 4-way Tewhey comparison (requires API, ~70 min)"
	@echo "  verify           Run canonical variant tests"
	@echo "  pipeline         Full scoring pipeline (requires API key, ~20 hrs)"

install:
	pip install -r requirements.txt

# ── Figure targets (use cached DB/parquet; no API calls) ──────────────────────

phase1: figures/saturation_cdf.png

figures/saturation_cdf.png: tewhey_mpra.parquet scored_variants.db tewhey_raw_delta_cache.db
	$(PYTHON) saturation_figure.py

phase2: figures/phase2_distribution.png figures/phase2_lfc_bins.png

figures/phase2_distribution.png figures/phase2_lfc_bins.png: tewhey_mpra.parquet tewhey_raw_delta_cache.db
	$(PYTHON) phase2_figures.py

phase3: figures/phase3_saturation_cdf.png

figures/phase3_saturation_cdf.png: phase3_blood_cache.db scored_variants.db tewhey_mpra.parquet
	$(PYTHON) phase3_blood_traits.py --from-cache

figures: phase1 phase2 phase3
	@echo "All figures generated."

# ── Raw delta runs (AlphaGenome API required) ────────────────────────────────

raw-stratified: tewhey_raw_delta_results.png TEWHEY_RESULT.md

tewhey_raw_delta_results.png TEWHEY_RESULT.md: tewhey_mpra.parquet
	$(PYTHON) extract_raw_deltas.py

raw-full: figures/tewhey_raw_delta_full_results.png

figures/tewhey_raw_delta_full_results.png: tewhey_mpra.parquet
	$(PYTHON) extract_raw_deltas.py --full

# ── Matched-statistic calibration (Components 2 + 3) ─────────────────────────
# Component 2: build_matched_calibration.py samples 5,000+ random common SNVs
#   and scores them through the same K562 max-over-tracks pipeline as Tewhey.
#   It is resume-able via matched_calibration_cache.db, so re-running is safe.
# Component 3: analyze_matched_calibration.py applies the matched null to the
#   Tewhey panel and produces the four-row Spearman comparison + figure.

matched_calibration:
	$(PYTHON) build_matched_calibration.py
	$(PYTHON) analyze_matched_calibration.py

# ── Verification ──────────────────────────────────────────────────────────────

verify:
	$(PYTHON) verification/canonical_variants_test.py

# ── Full pipeline (AlphaGenome API required) ──────────────────────────────────

pipeline:
	@echo "Step 1/4: Prefetch GWAS variants"
	$(PYTHON) prefetch_variants.py
	@echo "Step 2/4: Score disease GWAS variants (~10 hrs)"
	$(PYTHON) batch_score.py
	@echo "Step 3/4: Score Tewhey MPRA panel (~8 hrs)"
	$(PYTHON) tewhey_analysis.py
	@echo "Step 4/4: Extract raw expression deltas"
	$(PYTHON) extract_raw_deltas.py
