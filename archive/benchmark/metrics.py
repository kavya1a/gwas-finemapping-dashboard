"""Benchmark metrics: recall@K, auROC, AUPRC, negative control specificity.

Recall@K (credible-set framing)
--------------------------------
For a synthetic credible set of size N containing 1 pathogenic variant and
N-1 background (benign) variants, Recall@K is the fraction of pathogenic
variants that rank in the top K when sorted by our score.

Equivalently (without constructing every credible set explicitly):
  Recall@K in credible set of N = fraction of pathogenic variants with
  score exceeding the (1 - K/N)-quantile of the benign distribution.

We report for N=10 credible sets:
  Recall@1  = fraction above 90th percentile of benign
  Recall@5  = fraction above 50th percentile of benign
  Recall@10 = fraction above 10th percentile of benign (≈ above min benign)

auROC and AUPRC
---------------
Standard binary classification: pathogenic=1, benign=0.
Score = composite_score (or CADD PHRED, or -log10(GWAS p-value)).

Negative control specificity
-----------------------------
Fraction of benign variants scoring above the 75th percentile of the
pathogenic score distribution. Should be low (near 0.25 by chance) if
the tool is specific; higher = more false positives.
"""

from __future__ import annotations

import math
import numpy as np
import pandas as pd


def _safe_quantile(arr: list[float], q: float) -> float:
    if not arr:
        return 0.0
    return float(np.quantile(arr, q))


def recall_at_k(
    path_scores: list[float],
    benign_scores: list[float],
    credible_set_size: int = 10,
    k_values: tuple[int, ...] = (1, 5, 10),
) -> dict[str, float]:
    """Compute Recall@K for each K in k_values using the quantile shortcut."""
    results = {}
    for k in k_values:
        if k > credible_set_size:
            continue
        # threshold = (1 - k/N) quantile of benign distribution
        q = 1.0 - k / credible_set_size
        threshold = _safe_quantile(benign_scores, q)
        recall = (
            sum(1 for s in path_scores if s > threshold) / len(path_scores)
            if path_scores else 0.0
        )
        results[f"recall_at_{k}"] = round(recall, 4)
    return results


def auroc(
    path_scores: list[float],
    benign_scores: list[float],
) -> float:
    """Compute area under ROC curve (Wilcoxon-Mann-Whitney estimator)."""
    if not path_scores or not benign_scores:
        return 0.5
    # P(pathogenic score > benign score)
    wins = sum(
        1 for s_p in path_scores for s_b in benign_scores if s_p > s_b
    )
    ties = sum(
        0.5 for s_p in path_scores for s_b in benign_scores if s_p == s_b
    )
    return (wins + ties) / (len(path_scores) * len(benign_scores))


def auprc(
    path_scores: list[float],
    benign_scores: list[float],
) -> float:
    """Compute area under precision-recall curve via trapezoidal integration."""
    labels = [1] * len(path_scores) + [0] * len(benign_scores)
    scores = path_scores + benign_scores
    if not scores:
        return 0.0

    # Sort by score descending
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    n_pos = len(path_scores)
    n_total = len(labels)

    precisions, recalls = [], []
    tp = fp = 0
    for _, label in pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
        precisions.append(tp / (tp + fp))
        recalls.append(tp / n_pos if n_pos else 0)

    # Trapezoid integration
    area = 0.0
    for i in range(1, len(recalls)):
        area += (recalls[i] - recalls[i - 1]) * (precisions[i] + precisions[i - 1]) / 2
    return abs(area)


def negative_control_specificity(
    path_scores: list[float],
    benign_scores: list[float],
    percentile: float = 75.0,
) -> float:
    """Fraction of benign variants above the Nth percentile of pathogenic scores."""
    if not path_scores or not benign_scores:
        return float("nan")
    threshold = _safe_quantile(path_scores, percentile / 100.0)
    return sum(1 for s in benign_scores if s > threshold) / len(benign_scores)


def compute_all_metrics(
    path_scores: list[float],
    benign_scores: list[float],
    label: str = "",
) -> dict:
    """Run the full metric suite and return a single dict."""
    recall = recall_at_k(path_scores, benign_scores)
    auc = auroc(path_scores, benign_scores)
    auprc_val = auprc(path_scores, benign_scores)
    specificity = negative_control_specificity(path_scores, benign_scores)
    median_path = float(np.median(path_scores)) if path_scores else float("nan")
    median_benign = float(np.median(benign_scores)) if benign_scores else float("nan")

    return {
        "label": label,
        "n_pathogenic": len(path_scores),
        "n_benign": len(benign_scores),
        "median_pathogenic_score": round(median_path, 4),
        "median_benign_score": round(median_benign, 4),
        **recall,
        "auROC": round(auc, 4),
        "AUPRC": round(auprc_val, 4),
        "neg_ctrl_frac_above_p75": round(specificity, 4),
    }


def scores_from_df(
    df: pd.DataFrame,
    score_col: str,
    is_pathogenic_col: str = "is_pathogenic",
) -> tuple[list[float], list[float]]:
    """Split a DataFrame's scores into (path_scores, benign_scores)."""
    path = df[df[is_pathogenic_col] == 1][score_col].dropna().tolist()
    benign = df[df[is_pathogenic_col] == 0][score_col].dropna().tolist()
    return path, benign
