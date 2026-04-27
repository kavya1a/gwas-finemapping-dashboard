import pandas as pd


_WEIGHTS = {"expression_score": 0.5, "splicing_score": 0.3, "chromatin_score": 0.2}


def rank_variants(scored: list[dict]) -> pd.DataFrame:
    """Compute composite causal impact score and return sorted DataFrame."""
    df = pd.DataFrame(scored)

    for col in _WEIGHTS:
        if col not in df.columns:
            df[col] = None
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).abs()

    df["impact_score"] = sum(df[col] * w for col, w in _WEIGHTS.items())
    df = df.sort_values("impact_score", ascending=False).reset_index(drop=True)
    return df
