"""Late-fuse independently trained branch scores into linear phenotype heads."""

from __future__ import annotations

from typing import Any

import numpy as np

from mbs.training.transparent_baselines import (
    evaluate_multitask_predictions,
    fit_linear_multitask,
    predict_linear_multitask,
)


def concatenate_score_blocks(
    blocks: list[np.ndarray],
    *,
    sample_ids: list[str] | None = None,
) -> np.ndarray:
    """Column-stack score matrices ``[n_samples, n_scores_i]`` with matching rows."""
    if not blocks:
        raise ValueError("need at least one score block")
    mats = [np.asarray(b, dtype=np.float32) for b in blocks]
    n = mats[0].shape[0]
    for m in mats[1:]:
        if m.shape[0] != n:
            raise ValueError(f"score block row mismatch: {m.shape[0]} vs {n}")
    if sample_ids is not None and len(sample_ids) != n:
        raise ValueError("sample_ids length must match score rows")
    return np.concatenate(mats, axis=1)


def fit_late_fusion_heads(
    scores_train: np.ndarray,
    *,
    age: np.ndarray | None,
    age_mask: np.ndarray | None,
    tissue: np.ndarray | None,
    tissue_mask: np.ndarray | None,
    sex: np.ndarray | None,
    sex_mask: np.ndarray | None,
) -> dict[str, Any]:
    """Train-fold linear heads on concatenated branch scores."""
    return fit_linear_multitask(
        scores_train,
        age=age,
        age_mask=age_mask,
        tissue=tissue,
        tissue_mask=tissue_mask,
        sex=sex,
        sex_mask=sex_mask,
    )


def evaluate_late_fusion(
    *,
    scores_train: np.ndarray,
    scores_test: np.ndarray,
    age_train: np.ndarray | None,
    age_mask_train: np.ndarray | None,
    tissue_train: np.ndarray | None,
    tissue_mask_train: np.ndarray | None,
    sex_train: np.ndarray | None,
    sex_mask_train: np.ndarray | None,
    age_test: np.ndarray | None,
    age_mask_test: np.ndarray | None,
    tissue_test: np.ndarray | None,
    tissue_mask_test: np.ndarray | None,
    sex_test: np.ndarray | None,
    sex_mask_test: np.ndarray | None,
    study_ids_test: np.ndarray | None = None,
    platforms_test: np.ndarray | None = None,
    tissue_class_names: list[str] | None = None,
) -> dict[str, Any]:
    """Fit on train scores, evaluate on held-out scores."""
    models = fit_late_fusion_heads(
        scores_train,
        age=age_train,
        age_mask=age_mask_train,
        tissue=tissue_train,
        tissue_mask=tissue_mask_train,
        sex=sex_train,
        sex_mask=sex_mask_train,
    )
    preds = predict_linear_multitask(models, scores_test)
    metrics = evaluate_multitask_predictions(
        preds=preds,
        age=age_test,
        age_mask=age_mask_test,
        tissue=tissue_test,
        tissue_mask=tissue_mask_test,
        sex=sex_test,
        sex_mask=sex_mask_test,
        study_ids=study_ids_test,
        platforms=platforms_test,
        tissue_class_names=tissue_class_names,
    )
    return {
        "metrics": metrics,
        "n_score_features": int(scores_train.shape[1]),
        "n_train": int(scores_train.shape[0]),
        "n_test": int(scores_test.shape[0]),
    }
