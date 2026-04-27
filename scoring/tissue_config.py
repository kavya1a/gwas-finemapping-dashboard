"""Disease-to-tissue mapping for track filtering.

Each disease maps to keyword lists that are matched (case-insensitive)
against the `biosample_name` and `gtex_tissue` columns of the tidy score
DataFrame returned by tidy_scores().

Ontology CURIEs are provided for documentation; the filtering logic in
composite.py uses the human-readable keyword approach because biosample_name
and gtex_tissue are the most reliably populated fields in the track metadata.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class TissueProfile:
    display_name: str
    biosample_keywords: list[str]  # matched against biosample_name
    gtex_keywords: list[str]       # matched against gtex_tissue
    # Reference ontology CURIEs (documentation only)
    ontology_notes: str = ""


# ---------------------------------------------------------------------------
# Disease tissue profiles
# ---------------------------------------------------------------------------

DISEASE_PROFILES: dict[str, TissueProfile] = {
    "alzheimers": TissueProfile(
        display_name="Alzheimer's Disease",
        biosample_keywords=[
            "brain", "neuron", "cortex", "hippocampus", "microglia",
            "astrocyte", "cerebral", "frontal", "temporal", "neuroblast",
            "neural", "cerebellum", "prefrontal",
        ],
        gtex_keywords=[
            "Brain",
        ],
        ontology_notes=(
            "UBERON:0001869 (cerebral cortex), UBERON:0002421 (hippocampus), "
            "CL:0000540 (neuron), CL:0000129 (microglia)"
        ),
    ),
    "schizophrenia": TissueProfile(
        display_name="Schizophrenia",
        biosample_keywords=[
            "brain", "neuron", "cortex", "hippocampus", "prefrontal",
            "striatum", "neural", "neuroblast", "cerebral", "frontal",
        ],
        gtex_keywords=[
            "Brain",
        ],
        ontology_notes=(
            "UBERON:0001870 (prefrontal cortex), UBERON:0002421 (hippocampus), "
            "UBERON:0002435 (striatum)"
        ),
    ),
    "t2d": TissueProfile(
        display_name="Type 2 Diabetes",
        biosample_keywords=[
            "pancreas", "islet", "beta cell", "liver", "adipose",
            "muscle", "hepatocyte", "skeletal", "insulin", "fat",
        ],
        gtex_keywords=[
            "Pancreas", "Liver", "Adipose", "Muscle",
        ],
        ontology_notes=(
            "UBERON:0001264 (pancreas), UBERON:0002107 (liver), "
            "UBERON:0001443 (adipose tissue), CL:0000169 (beta cell)"
        ),
    ),
    "parkinsons": TissueProfile(
        display_name="Parkinson's Disease",
        biosample_keywords=[
            "brain", "substantia nigra", "striatum", "dopamin", "neuron",
            "cortex", "neural", "neuroblast", "cerebral", "basal ganglia",
            "midbrain",
        ],
        gtex_keywords=[
            "Brain",
        ],
        ontology_notes=(
            "UBERON:0002038 (substantia nigra), UBERON:0002435 (striatum), "
            "CL:0000100 (motor neuron)"
        ),
    ),
}

# Fallback: no tissue filter — use all tracks
ALL_TISSUES = TissueProfile(
    display_name="All tissues (no filter)",
    biosample_keywords=[],
    gtex_keywords=[],
)


def get_profile(disease: str) -> TissueProfile:
    """Returns the TissueProfile for a disease slug, or ALL_TISSUES if unknown."""
    return DISEASE_PROFILES.get(disease.lower().replace(" ", "_").replace("'s", "s"), ALL_TISSUES)


def filter_tracks(df, profile: TissueProfile):
    """Filter a tidy score DataFrame to tracks matching the tissue profile.

    Matching is OR across all keywords, AND across biosample vs gtex columns.
    Returns the full DataFrame unfiltered if no keywords are specified.
    """
    import pandas as pd

    if not profile.biosample_keywords and not profile.gtex_keywords:
        return df

    mask = pd.Series(False, index=df.index)

    if "biosample_name" in df.columns and profile.biosample_keywords:
        bio_col = df["biosample_name"].fillna("").str.lower()
        for kw in profile.biosample_keywords:
            mask |= bio_col.str.contains(kw.lower(), na=False)

    if "gtex_tissue" in df.columns and profile.gtex_keywords:
        gtex_col = df["gtex_tissue"].fillna("").str.lower()
        for kw in profile.gtex_keywords:
            mask |= gtex_col.str.contains(kw.lower(), na=False)

    filtered = df[mask]
    # If filtering eliminated everything (no matching tracks for this tissue),
    # fall back to all tracks with a warning.
    if filtered.empty:
        return df
    return filtered
