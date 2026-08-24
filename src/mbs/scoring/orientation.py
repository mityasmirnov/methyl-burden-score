"""ADR 0008 score-orientation anchor (predictive MBS, not constraint/LOEUF)."""

from __future__ import annotations

from typing import Any

import numpy as np

ANCHOR_RECIPE = "pearson_vs_signed_gene_mean_mvalue"


def polarity_from_correlation(mbs: np.ndarray, signed_m: np.ndarray) -> str:
    """Return ``hyper_aligned`` or ``flipped`` from Pearson r of gene-mean scores."""
    a = np.asarray(mbs, dtype=np.float64).reshape(-1)
    b = np.asarray(signed_m, dtype=np.float64).reshape(-1)
    if a.shape != b.shape or a.size < 2:
        raise ValueError("mbs and signed_m must be 1-D and same length >= 2")
    if float(np.std(a)) == 0 or float(np.std(b)) == 0:
        return "hyper_aligned"
    r = float(np.corrcoef(a, b)[0, 1])
    return "flipped" if r < 0 else "hyper_aligned"


def apply_orientation(
    mbs: np.ndarray,
    *,
    signed_m: np.ndarray,
    head_weights: np.ndarray | None = None,
) -> dict[str, Any]:
    """Flip MBS (and optional head weights) so higher MBS tracks hypermethylation."""
    polarity = polarity_from_correlation(mbs.mean(axis=0) if mbs.ndim == 2 else mbs, signed_m)
    scores = np.asarray(mbs, dtype=np.float64)
    weights = None if head_weights is None else np.asarray(head_weights, dtype=np.float64).copy()
    if polarity == "flipped":
        scores = 1.0 - scores
        if weights is not None:
            weights = -weights
    return {
        "mbs": scores.astype(np.float32, copy=False),
        "head_weights": weights,
        "score_polarity": polarity,
        "anchor_recipe": ANCHOR_RECIPE,
        "score_family": "predictive_mbs",
    }


def score_manifest(
    *,
    score_polarity: str,
    fold_id: str | None = None,
    restart_id: str | None = None,
) -> dict[str, Any]:
    return {
        "artifact_version": "1",
        "score_family": "predictive_mbs",
        "score_polarity": score_polarity,
        "anchor_recipe": ANCHOR_RECIPE,
        "fold_id": fold_id,
        "restart_id": restart_id,
        "notes": "Predictive sample x gene MBS; not a methylation constraint or LOEUF analogue.",
    }
