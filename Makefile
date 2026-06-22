.PHONY: install reproduce figures saturation distribution blood hero \
        matched-figures analysis test verify \
        raw-full matched_calibration pipeline \
        dnase-atac phred-check help

PYTHON := python3.11

# Determinism: window sampling uses seed = 2026 (config.yaml -> matched_calibration.seed);
# every bootstrap CI uses SEED = 42 inside the analysis scripts. No wall-clock or
# network input enters the reproduce path, so `make reproduce` is byte-stable for the
# CSV numbers (figure PNGs may differ only in matplotlib rendering metadata).

help:
	@echo "Targets (cache-only, no API key needed):"
	@echo "  reproduce            Regenerate every figure + every README number from cache, then check them"
	@echo "  figures              Regenerate all figures from cached data"
	@echo "  test                 Run the smoke suite (asserts headline numbers + shapes)"
	@echo "  verify               Run canonical variant tests (5/5)"
	@echo ""
	@echo "Targets (AlphaGenome API key required; hours of runtime):"
	@echo "  raw-full             Score the full Tewhey panel from scratch (~8 hrs)"
	@echo "  matched_calibration  Build the matched null from scratch + Tewhey re-analysis (~70 min)"
	@echo "  pipeline             Full from-scratch GWAS + Tewhey pipeline (~20 hrs)"
	@echo ""
	@echo "Author-suggested next experiments (scaffolded, NOT yet run; fail loud without inputs):"
	@echo "  dnase-atac           Matched re-calibration under DNase/ATAC tracks"
	@echo "  phred-check          Compare SDK phred values vs single-track / matched transforms"

install:
	pip install -r requirements.txt

# ════════════════════════════════════════════════════════════════════════════
#  REPRODUCE — the one entry point. No API key. Deterministic.
#  Regenerates every figure and every CSV the README cites from committed cache,
#  then runs the smoke suite to confirm the regenerated numbers match the text.
# ════════════════════════════════════════════════════════════════════════════

reproduce: figures test
	@echo ""
	@echo "✓ reproduce complete — all figures + numbers regenerated from committed cache."

# ── Figures + CSV tables (all cache-only) ────────────────────────────────────

figures: saturation distribution blood hero matched-figures
	@echo "All figures regenerated from cache."

saturation:
	$(PYTHON) saturation_figure.py

distribution:
	$(PYTHON) distribution_figures.py

blood:
	$(PYTHON) blood_trait_replication.py --from-cache

hero:
	$(PYTHON) make_hero_figure.py

# Matched-calibration figures + the four-row / three-recipe / mean-vs-max CSV tables.
# build --post rebuilds the null parquet + histogram from matched_calibration_cache.db
# (no API); the analyze_* scripts then re-quantile Tewhey against that null.
matched-figures:
	$(PYTHON) build_matched_calibration.py --post
	$(PYTHON) analyze_matched_calibration.py
	$(PYTHON) analyze_matched_calibration_recipes.py
	$(PYTHON) analyze_mean_aggregation.py

analysis: matched-figures

# ── Smoke suite + canonical tests ────────────────────────────────────────────

test:
	$(PYTHON) -m pytest tests/ -q

verify:
	$(PYTHON) verification/canonical_variants_test.py

# ════════════════════════════════════════════════════════════════════════════
#  FROM-SCRATCH PATHS — require ALPHAGENOME_API_KEY and hours of API time.
# ════════════════════════════════════════════════════════════════════════════

# Full Tewhey panel raw-delta extraction (writes tewhey_raw_delta_cache.db + figure).
raw-full:
	$(PYTHON) extract_raw_deltas.py --full

# Build the matched null on common variants, then run the 4-way Tewhey comparison.
# Resume-able via matched_calibration_cache.db.
matched_calibration:
	$(PYTHON) build_matched_calibration.py
	$(PYTHON) analyze_matched_calibration.py

# Full from-scratch: prefetch GWAS -> score GWAS -> Tewhey raw deltas (~20 hrs API).
#   prefetch_variants.py    GWAS Catalog -> preloaded_variants.db
#   batch_score.py          AlphaGenome  -> scored_variants.db
#   extract_raw_deltas.py   Tewhey raw   -> tewhey_raw_delta_cache.db
pipeline: scored_variants.db figures/tewhey_raw_delta_full_results.png

preloaded_variants.db:
	$(PYTHON) prefetch_variants.py

scored_variants.db: preloaded_variants.db
	$(PYTHON) batch_score.py

figures/tewhey_raw_delta_full_results.png: tewhey_mpra.parquet
	$(PYTHON) extract_raw_deltas.py --full

# ════════════════════════════════════════════════════════════════════════════
#  AUTHOR-SUGGESTED NEXT EXPERIMENTS — scaffolded, not yet run.
#  Both fail loudly with an explanatory error if their required inputs are
#  absent (API key / a newer SDK's phred values). They never emit simulated
#  numbers. See the README "Next experiments" subsection.
# ════════════════════════════════════════════════════════════════════════════

dnase-atac:
	$(PYTHON) dnase_atac_recalibration.py

phred-check:
	$(PYTHON) phred_scale_check.py
