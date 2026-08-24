"""Negative-control transforms (static-only, coverage-only, label permutation)."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from random import Random

import numpy as np
from numpy.typing import NDArray
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import OneHotEncoder

from mbs.evaluation.metrics import multiclass_metrics, regression_metrics
from mbs.training.phenotypes import SamplePhenotype


def apply_feature_control(
    features: NDArray[np.float32],
    *,
    mode: str,
    include_m_value: bool = True,
) -> NDArray[np.float32]:
    """Zero methylation and/or static channels. Layout: beta, [M], static..., static_present."""
    out = np.asarray(features, dtype=np.float32).copy()
    if out.ndim != 2 or out.shape[0] == 0:
        return out
    n_methyl = 2 if include_m_value else 1
    if mode in {"none", "off", ""}:
        return out
    if mode == "static_only":
        out[:, :n_methyl] = 0.0
        return out
    if mode == "coverage_only":
        out[:, :-1] = 0.0
        return out
    raise ValueError(f"unknown feature control: {mode}")


def permute_labels_within_study(
    phenotypes: Sequence[SamplePhenotype],
    *,
    seed: int,
) -> list[SamplePhenotype]:
    """Shuffle class_index / age / sex within study strata."""
    rng = Random(seed)  # noqa: S311
    by_study: dict[str, list[int]] = defaultdict(list)
    for i, ph in enumerate(phenotypes):
        by_study[str(ph.study_id or ph.sample_id)].append(i)
    ages = [p.age for p in phenotypes]
    classes = [p.class_index for p in phenotypes]
    sexes = [p.sex_class_index for p in phenotypes]
    for idxs in by_study.values():
        order = list(range(len(idxs)))
        rng.shuffle(order)
        src_ages = [ages[i] for i in idxs]
        src_cls = [classes[i] for i in idxs]
        src_sex = [sexes[i] for i in idxs]
        for k, i in enumerate(idxs):
            ages[i] = src_ages[order[k]]
            classes[i] = src_cls[order[k]]
            sexes[i] = src_sex[order[k]]
    return [
        SamplePhenotype(
            sample_id=ph.sample_id,
            cell_type=ph.cell_type,
            donor_id=ph.donor_id,
            title=ph.title,
            class_index=classes[i],
            study_id=ph.study_id,
            age=ages[i],
            platform=ph.platform,
            age_mask=ph.age_mask,
            tissue_mask=ph.tissue_mask,
            sex_mask=ph.sex_mask,
            sex_class_index=sexes[i],
        )
        for i, ph in enumerate(phenotypes)
    ]


def fit_metadata_only(
    *,
    study_ids: Sequence[str],
    platforms: Sequence[str | None],
    tissues: Sequence[str | None],
    y: np.ndarray,
    task: str,
) -> dict[str, float]:
    """Linear/logistic ceiling from study/platform/tissue one-hots (not the encoder)."""
    cats = np.column_stack(
        [
            np.asarray(study_ids, dtype=object),
            np.asarray([p or "unknown" for p in platforms], dtype=object),
            np.asarray([t or "unknown" for t in tissues], dtype=object),
        ]
    )
    x = OneHotEncoder(handle_unknown="ignore").fit_transform(cats)
    y_arr = np.asarray(y)
    if task == "regression":
        pred = Ridge(alpha=1.0).fit(x, y_arr).predict(x)
        return regression_metrics(y_arr, pred)
    pred = LogisticRegression(max_iter=200).fit(x, y_arr).predict(x)
    return {
        k: v
        for k, v in multiclass_metrics(y_arr, pred).items()
        if k in {"macro_f1", "balanced_accuracy"}
    }
