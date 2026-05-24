.PHONY: install figures phase1 phase2 phase3 raw-full matched_calibration verify help

PYTHON := python3.11

help:
	@echo "Targets:"
	@echo "  install              Install dependencies"
	@echo "  figures              Regenerate all figures from cached data (fast)"
	@echo "  phase1               Saturation CDF figure"
	@echo "  phase2               Distribution + LFC-bin correlation figures"
	@echo "  phase3               Blood trait replication figure"
	@echo "  raw-full             Score full 3,259-variant Tewhey panel (requires API, ~8 hrs)"
	@echo "  matched_calibration  Build matched null + run 4-way Tewhey comparison (requires API, ~70 min)"
	@echo "  verify               Run canonical variant tests"

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

# ── Raw delta extraction (AlphaGenome API required) ──────────────────────────

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
