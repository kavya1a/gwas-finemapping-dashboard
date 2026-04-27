#!/usr/bin/env bash
# Overnight orchestration: Tasks 1–4.
# Stops at the first task that writes to OVERNIGHT_BLOCKERS.md.

set -euo pipefail

PYTHON=/opt/homebrew/bin/python3.11
DIR="$(cd "$(dirname "$0")" && pwd)"
BLOCKERS="$DIR/OVERNIGHT_BLOCKERS.md"
LOG="$DIR/overnight.log"

cd "$DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

check_blockers() {
    # Only halt if the blockers file has content beyond the header sentinel
    if grep -q "^## BLOCKER" "$BLOCKERS" 2>/dev/null; then
        log "BLOCKER detected — halting overnight run. See $BLOCKERS"
        exit 1
    fi
}

# Clear previous blockers file so we start clean
cat > "$BLOCKERS" << 'HEADER'
# Overnight Blockers

_This file is cleared at the start of each run_overnight.sh execution.
Blockers are appended here by individual scripts when they hit unexpected
conditions requiring human judgment._

_If this file has content below the separator, the overnight run stopped early._

---
HEADER

log "=== Overnight run started ==="

# ── Task 2: Prefetch (Task 1 allele resolution runs inline) ──────────────────
log "--- Task 2: Prefetch top 100 GWAS variants per disease ---"
$PYTHON prefetch_variants.py 2>&1 | tee -a "$LOG"
check_blockers

# ── Task 3: Batch-score all diseases (indels already filtered in DB) ─────────
log "--- Task 3: Batch-score all diseases through AlphaGenome ---"
$PYTHON batch_score.py 2>&1 | tee -a "$LOG"
check_blockers

# ── Pipeline yield report ─────────────────────────────────────────────────────
log "--- Pipeline yield report ---"
$PYTHON pipeline_yield_report.py 2>&1 | tee -a "$LOG"

# ── Task 4: Tewhey MPRA (only if Tasks 1–3 complete cleanly) ─────────────────
log "--- Task 4: Download and process Tewhey MPRA data ---"
$PYTHON fetch_tewhey_mpra.py 2>&1 | tee -a "$LOG"
check_blockers

log "=== Overnight run complete — no blockers ==="
