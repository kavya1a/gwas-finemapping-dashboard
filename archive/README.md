# Archive

Files moved here when the project pivoted to the calibration-statistic / matched-recalibration finding.

**Nothing in this directory is part of the current result or the `make reproduce` path.** It is pre-pivot work, kept for provenance only — earlier scripts, superseded write-ups (e.g. `TEWHEY_RESULT.md`, the raw-delta-primary framing; `OVERNIGHT_BLOCKERS.md`, a historical run-log), and the old random-variant control (superseded by the matched-statistic null in the main repo). Numbers in these files reflect earlier runs and are not maintained against the committed caches.

| File / Dir | What it was | Why cut |
|---|---|---|
| `app.py` | Streamlit dashboard UI | UI explicitly out of scope (Path 2 deliverable is README + repo, no UI) |
| `alphagenome_scorer.py` | Early scoring wrapper, predates batch_score.py | Superseded by batch_score.py; not referenced by any active script |
| `ranking.py` | GWAS PIP-weighted ranking prototype | GWAS PIP correlation cut from scope |
| `run_ad_test.py` | Ad-hoc AD scoring test script | One-off diagnostic, not part of reproducible pipeline |
| `run_generalization_test.py` | Early generalization test scaffold | Superseded by Phase 3 plan in PROJECT_LOCK.md |
| `fetch_tewhey_mpra.py` | Standalone Tewhey downloader | Logic absorbed into tewhey_analysis.py |
| `run_overnight.sh` | Overnight batch runner shell script | Pipeline now run directly via Python scripts |
| `findings.md` | Early draft findings (pre-pivot) | Replaced by TEWHEY_RESULT.md and eventual README; contained unverified claims |
| `benchmark/` | Benchmark harness directory | Benchmark work cut from Path 2 scope |
| `benchmark_report.md` | Benchmark summary | Cut with benchmark work |
| `benchmark_results.json` | Benchmark data | Cut with benchmark work |
