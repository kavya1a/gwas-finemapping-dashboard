# Tewhey Analysis Flags

_Stop conditions are written here. If this file has content below the separator, the analysis halted early._

---

## FLAG [2026-04-21 16:50:17]: CADD API unavailable (probe returned 0/5)
All 5 probe requests returned empty results. CADD may be rate-limiting or temporarily down. Skipping remaining CADD lookups and continuing to AlphaGenome scoring. Re-run with CADD cache populated to include CADD scores in the report.

## FLAG [2026-04-22 00:31:15]: expression_subscore correlation < 0.10
Spearman r = 0.0361 (n=3259) is below the 0.10 threshold. Possible causes: wrong allele orientation (A/B vs ref/alt mismatch), tissue mismatch (LCL vs blood-lineage filter), or model calibration issue. Do NOT report until diagnosed.
