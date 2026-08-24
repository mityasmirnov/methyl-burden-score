"""Sparse direct CpG branch: elastic-net over centered M-values."""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.linear_model import ElasticNet


def locus_cross_study_mask(
    observed: np.ndarray,
    study_ids: np.ndarray,
    *,
    min_studies: int = 2,
) -> np.ndarray:
    """``observed`` is bool [n_samples, n_loci]; keep loci seen in >= min_studies."""
    obs = np.asarray(observed, dtype=bool)
    studies = np.asarray(study_ids)
    keep = np.zeros(obs.shape[1], dtype=bool)
    for j in range(obs.shape[1]):
        keep[j] = len(set(studies[obs[:, j]].tolist())) >= int(min_studies)
    return keep


def fit_direct_elasticnet(
    z: np.ndarray,
    observed: np.ndarray,
    y: np.ndarray,
    study_ids: np.ndarray,
    *,
    min_studies: int = 2,
    alpha: float = 0.1,
    l1_ratio: float = 0.5,
) -> dict[str, Any]:
    """Fit D(s)=sum_c w_c z_s,c on observed centered M; missing z treated as 0."""
    z_arr = np.asarray(z, dtype=np.float64)
    obs = np.asarray(observed, dtype=bool)
    if z_arr.shape != obs.shape:
        raise ValueError("z and observed shape mismatch")
    keep = locus_cross_study_mask(obs, study_ids, min_studies=min_studies)
    if not keep.any():
        raise ValueError("no loci pass min_studies coverage")
    x = np.where(obs, z_arr, 0.0)[:, keep]
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=2000)
    model.fit(x, y_arr)
    pred = model.predict(x)
    return {
        "weights": model.coef_.astype(np.float32),
        "intercept": float(model.intercept_),
        "keep_mask": keep,
        "pred": pred.astype(np.float32),
        "n_loci": int(keep.sum()),
    }
