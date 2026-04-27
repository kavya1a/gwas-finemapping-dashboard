"""Per-modality variant scorer configurations.

Maps scorer names to SDK objects with the correct aggregation type and window
size, as recommended by the AlphaGenome paper and SDK docs.
"""

from alphagenome.models import variant_scorers as vs
from alphagenome.models.dna_model import Organism

# --- Differential scorers (ALT vs REF) ---
# Each entry is the SDK scorer object with window and aggregation pre-set.
DIFFERENTIAL_SCORERS: dict[str, vs.VariantScorerTypes] = {
    "rna_seq": vs.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
    # CAGE: log2-ratio of summed signals in 501-bp window (captures TSS activity)
    "cage": vs.RECOMMENDED_VARIANT_SCORERS["CAGE"],
    # PRO-cap: same math as CAGE; human-only
    "procap": vs.RECOMMENDED_VARIANT_SCORERS["PROCAP"],
    # ATAC: log2-ratio of summed signals in 501-bp window
    "atac": vs.RECOMMENDED_VARIANT_SCORERS["ATAC"],
    # DNase: same formula as ATAC
    "dnase": vs.RECOMMENDED_VARIANT_SCORERS["DNASE"],
    # ChIP-TF: 501-bp window
    "chip_tf": vs.RECOMMENDED_VARIANT_SCORERS["CHIP_TF"],
    # ChIP-Histone: 2001-bp window (marks spread further than TF binding)
    "chip_histone": vs.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE"],
    # Splicing — three separate scorers kept distinct per spec
    "splice_sites": vs.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
    "splice_site_usage": vs.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
    "splice_junctions": vs.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    # 3D genome contacts: mean |ALT - REF| contact freq in 1 Mb
    "contact_maps": vs.RECOMMENDED_VARIANT_SCORERS["CONTACT_MAPS"],
    # Polyadenylation: max |log-FC(isoform ratio)| in 400-bp window; human-only
    "polyadenylation": vs.RECOMMENDED_VARIANT_SCORERS["POLYADENYLATION"],
}

# --- Active scorers: max(REF, ALT) — not differential ---
# Used as activity priors: a large differential score in a silent region
# is suspect. Active score gates trustworthiness of differential scores.
ACTIVE_SCORERS: dict[str, vs.VariantScorerTypes] = {
    "rna_seq_active": vs.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ_ACTIVE"],
    "atac_active": vs.RECOMMENDED_VARIANT_SCORERS["ATAC_ACTIVE"],
    "dnase_active": vs.RECOMMENDED_VARIANT_SCORERS["DNASE_ACTIVE"],
    "chip_tf_active": vs.RECOMMENDED_VARIANT_SCORERS["CHIP_TF_ACTIVE"],
    "chip_histone_active": vs.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE_ACTIVE"],
    "cage_active": vs.RECOMMENDED_VARIANT_SCORERS["CAGE_ACTIVE"],
    "procap_active": vs.RECOMMENDED_VARIANT_SCORERS["PROCAP_ACTIVE"],
}

# --- Modality groups for composite scoring ---
# Splicing sub-scores are kept separate internally and combined with equal
# weight into a single "splicing_impact" group in composite.py.
MODALITY_GROUPS: dict[str, list[str]] = {
    "expression": ["rna_seq", "cage", "procap"],
    "chromatin": ["atac", "dnase", "chip_histone"],
    "tf_binding": ["chip_tf"],
    "splicing": ["splice_sites", "splice_site_usage", "splice_junctions"],
    "contact": ["contact_maps"],
    "polyadenylation": ["polyadenylation"],
}

# Human-only scorers (PROCAP and POLYADENYLATION are not available for mouse)
HUMAN_ONLY_SCORERS = {"procap", "polyadenylation"}


def get_scorers_for_organism(organism: Organism = Organism.HOMO_SAPIENS) -> list:
    """Returns differential scorer list filtered to the given organism."""
    if organism == Organism.HOMO_SAPIENS:
        return list(DIFFERENTIAL_SCORERS.values())
    return [
        scorer
        for name, scorer in DIFFERENTIAL_SCORERS.items()
        if name not in HUMAN_ONLY_SCORERS
    ]


def get_active_scorers_for_organism(organism: Organism = Organism.HOMO_SAPIENS) -> list:
    if organism == Organism.HOMO_SAPIENS:
        return list(ACTIVE_SCORERS.values())
    return [
        scorer
        for name, scorer in ACTIVE_SCORERS.items()
        if name.replace("_active", "") not in HUMAN_ONLY_SCORERS
    ]
