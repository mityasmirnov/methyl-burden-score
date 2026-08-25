"""Sparse direct CpG branch: elastic-net over Level-1 z or centered M."""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.linear_model import ElasticNet

from mbs.training.level1_norm import Level1NormParams, apply_level1_robust_z, fit_level1_robust_z


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


def centered_m_features(
    m: np.ndarray,
    observed: np.ndarray,
) -> np.ndarray:
    """Column-center finite M on observed cells; missing → 0 (7C placeholder)."""
    m_arr = np.asarray(m, dtype=np.float64)
    obs = np.asarray(observed, dtype=bool)
    if m_arr.shape != obs.shape:
        raise ValueError("m and observed shape mismatch")
    out = np.zeros_like(m_arr, dtype=np.float64)
    for j in range(m_arr.shape[1]):
        vals = m_arr[obs[:, j], j]
        if vals.size == 0:
            continue
        mu = float(np.mean(vals))
        out[obs[:, j], j] = m_arr[obs[:, j], j] - mu
    return out


def level1_z_features(
    m: np.ndarray,
    observed: np.ndarray,
    *,
    sigma_min: float = 1e-6,
    epsilon: float = 0.001,
    params: Level1NormParams | None = None,
) -> tuple[np.ndarray, Level1NormParams]:
    """Fit (or apply) Level-1 robust-z; novel/unestimated → 0 with no presence."""
    m_arr = np.asarray(m, dtype=np.float64)
    obs = np.asarray(observed, dtype=bool)
    if m_arr.shape != obs.shape:
        raise ValueError("m and observed shape mismatch")
    if params is None:
        params = fit_level1_robust_z(
            m_arr,
            observed=obs,
            sigma_min=sigma_min,
            epsilon=epsilon,
        )
    z_rows: list[np.ndarray] = []
    for i in range(m_arr.shape[0]):
        z_i, _present = apply_level1_robust_z(m_arr[i], params)
        z_rows.append(np.asarray(z_i, dtype=np.float64))
    return np.stack(z_rows, axis=0), params


def direct_cpg_design_matrix(
    m: np.ndarray,
    observed: np.ndarray,
    *,
    use_level1: bool,
    sigma_min: float = 1e-6,
    epsilon: float = 0.001,
    level1_params: Level1NormParams | None = None,
) -> tuple[np.ndarray, Level1NormParams | None]:
    """Build the z matrix for direct elastic-net (Level-1 when channel B)."""
    if use_level1:
        z, params = level1_z_features(
            m,
            observed,
            sigma_min=sigma_min,
            epsilon=epsilon,
            params=level1_params,
        )
        return z, params
    return centered_m_features(m, observed), None


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
    """Fit D(s)=sum_c w_c z_s,c; missing z treated as 0.

    ``z`` is Level-1 robust-z when channel B is on, else centered M (7C).
    """
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
