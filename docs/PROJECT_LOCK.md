PROJECT_LOCK.md — Path 2: Saturation finding pivot
Goal: Demonstrate that quantile-normalized AlphaGenome-based variant scoring saturates on regulatory-enriched variant sets, with implications for how such tools should be validated. Frame project around this novel methodological finding.
Timeline: 5-7 days from today. Hard stop at day 7 regardless of state.
Deliverable: GitHub repo with rigorous README as primary writeup. No formal preprint. No UI.
Headline claim (working): "Quantile-normalized regulatory variant scoring saturates on enrichment-selected variant sets, with magnitude correlation preserved but signed correlation collapsed. We demonstrate this on AlphaGenome-derived scores against Tewhey 2016 MPRA, and show that raw (un-normalized) deltas restore the directional signal."

What's locked (don't change)

Pipeline architecture, scoring code, config
363 scored GWAS variants across 4 diseases
Verification harness (3/5 pass, accept as-is)
APOE ε4/ε2 sign inversion result
Tewhey diagnostic findings (the four hypotheses + saturation conclusion)
Decision to use Option C (raw deltas) for primary correlation

What's in scope for the pivot (5-7 days of work)
Phase 1: Lock the saturation evidence (2 days)
The core finding needs to be airtight. Three concrete tasks:

Raw vs normalized comparison on Tewhey (in progress per Option C prompt). Lock these four numbers.
Quantify the saturation directly. Plot the empirical CDF of expression_subscore for: (a) all Tewhey variants, (b) random genome-wide common variants, (c) your 363 GWAS variants. Show that Tewhey is shifted entirely to the tail. This is the figure of the paper.
Show this isn't unique to AlphaGenome's expression modality. Run the same saturation analysis on chromatin and TF binding subscores. If they all saturate, the finding generalizes. If only expression saturates, that's also a finding (and an interesting one).

Phase 2: Demonstrate the fix (1-2 days)
Show that the raw-delta approach actually produces a continuous predictor, not just better correlation:

Distribution of raw deltas vs distribution of normalized scores on Tewhey — visual proof that raw deltas are continuous while normalized scores are bimodal.
Correlation by |LFC| bin for both: show that raw deltas correlate even in the low-|LFC| range where normalized scores fail.
The methodological recommendation: "For validation against continuous experimental measurements, use raw model outputs. For ranking within enrichment-selected variant sets, use normalized scores. They serve different purposes and one is not strictly better than the other."

Phase 3: Generalize the claim (1-2 days, optional but lifts ceiling)
This is what could push the project from "interesting finding" to "actually novel contribution":

Demonstrate the saturation occurs in other regulatory-enriched datasets, not just Tewhey. Suggested: a second public MPRA dataset (Ulirsch 2016 or Inoue 2019), or eQTL credible sets from GTEx, or fine-mapped disease variants from any source.
If it does, the claim becomes: "saturation is a property of how researchers select test sets for regulatory predictors, not a property of any specific tool."
If saturation doesn't replicate on a second dataset, that's also informative — the finding becomes specific to Tewhey-style screens, and you frame the paper accordingly.

This phase is the difference between "we found something on our data" and "we found something with general implications." Time-box at 2 days.
Phase 4: Writeup (1-2 days)
README structure, locked:

Title + 1-sentence summary
Headline figure: the empirical CDF showing saturation
Background: what AlphaGenome is, what regulatory variant prioritization tools do, how they're typically validated
Finding 1: Saturation effect on regulatory-enriched datasets (Tewhey)
Finding 2: Magnitude correlation is preserved, signed correlation is not
Finding 3: Raw deltas restore the directional signal
Finding 4 (if Phase 3 succeeds): Generalization to other datasets
Implications: Methodological recommendations for tool validation
Reproducibility: clear instructions to regenerate every figure and number
Limitations: honest catalog
Code & data: what's in the repo

The pipeline + canonical variant validation stays in the repo as supporting work. The README leads with the saturation finding because that's the novel claim.

What's cut (do not work on)

UI / dashboard
GWAS PIP correlation
CADD baseline beyond what local lookup gives in 1hr
SNCA tissue retuning
TCF7L2 investigation
5/5 canonical variant target
Indel support
Allele-frequency-matched null
Anything not directly serving the saturation framing


Decision gates
Gate 1 (end of Phase 1): Does the saturation evidence hold across all three modalities, and does the raw-delta correlation come back at ρ ≥ 0.15? If yes → proceed. If no → fall back to Path 1 framing, ship in 2 more days.
Gate 2 (end of Phase 2): Does the raw-delta vs normalized comparison produce visually compelling figures? Show them to one outside person — does it land? If yes → proceed to Phase 3. If no → simplify framing or fall back to Path 1.
Gate 3 (end of Phase 3 or day 5, whichever first): Does the finding generalize to a second dataset? If yes → strong paper. If no → frame paper as Tewhey-specific. Either way, proceed to writeup.
Gate 4 (day 6): Is the README 70% drafted? If not, cut Phase 3 entirely and write what you have.
Day 7: Ship regardless of state.
