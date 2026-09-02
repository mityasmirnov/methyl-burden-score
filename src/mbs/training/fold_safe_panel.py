"""Fold-safe probe panel selection for 7G′ Stage B (C-mvalue-enetS)."""

from __future__ import annotations

import hashlib
from typing import Any, Literal

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import SGDClassifier, SGDRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mbs.training.cascade_assign import CascadeAssignment

TaskKind = Literal["age", "sex", "tissue"]
ENET_ALPHA_GRID = (1e-5, 1e-4, 1e-3)
ENET_L1_GRID = (0.25, 0.5, 0.75)


def _enet_classifier_pipeline(*, alpha: float, l1_ratio: float) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler(with_mean=True)),
            (
                "sgd",
                SGDClassifier(
                    loss="log_loss",
                    penalty="elasticnet",
                    alpha=alpha,
                    l1_ratio=l1_ratio,
                    max_iter=60,
                    tol=1e-3,
                    random_state=42,
                    n_jobs=-1,
                ),
            ),
        ]
    )


def _enet_regressor_pipeline(*, alpha: float, l1_ratio: float) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler(with_mean=True)),
            (
                "sgd",
                SGDRegressor(
                    loss="squared_error",
                    penalty="elasticnet",
                    alpha=alpha,
                    l1_ratio=l1_ratio,
                    max_iter=60,
                    tol=1e-3,
                    random_state=42,
                ),
            ),
        ]
    )


def _study_inner_folds(
    study_ids: np.ndarray,
    *,
    n_inner_folds: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Study-grouped inner splits: no study appears in both inner-train and inner-val."""
    studies = np.unique(np.asarray(study_ids, dtype=object))
    if studies.size < 2:
        idx = np.arange(len(study_ids), dtype=np.int64)
        return [(idx, idx)]
    rng = np.random.default_rng(seed)
    order = studies.copy()
    rng.shuffle(order)
    fold_of_study = {str(s): i % n_inner_folds for i, s in enumerate(order.tolist())}
    groups = np.asarray([fold_of_study[str(s)] for s in study_ids.tolist()], dtype=np.int64)
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for holdout in range(n_inner_folds):
        val_idx = np.flatnonzero(groups == holdout).astype(np.int64)
        train_idx = np.flatnonzero(groups != holdout).astype(np.int64)
        if train_idx.size == 0 or val_idx.size == 0:
            continue
        splits.append((train_idx, val_idx))
    return splits or [(np.arange(len(study_ids), dtype=np.int64), np.arange(len(study_ids), dtype=np.int64))]


def _coef_abs(model: Pipeline, n_cols: int) -> np.ndarray:
    raw_coef = model.named_steps["sgd"].coef_
    if raw_coef.ndim == 2:
        return np.max(np.abs(raw_coef), axis=0)
    return np.abs(raw_coef.ravel())


def stability_select_columns(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    study_ids: np.ndarray | None = None,
    task: TaskKind = "tissue",
    max_seeds: int = 10_000,
    n_inner_folds: int = 3,
    n_repeats: int = 5,
    min_frequency: float = 0.34,
    seed: int = 42,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Repeated study-grouped inner-CV elastic-net; rank columns by selection frequency."""
    if x_train.shape[0] < 4:
        picked = np.arange(min(max_seeds, x_train.shape[1]), dtype=np.int64)
        return picked, {"n_runs": 0, "frequency": {}, "mean_abs_coef": {}}
    n_cols = x_train.shape[1]
    counts = np.zeros(n_cols, dtype=np.int64)
    coef_sum = np.zeros(n_cols, dtype=np.float64)
    n_runs = 0
    if study_ids is None:
        study_ids = np.array(["NA"] * x_train.shape[0], dtype=object)
    for repeat in range(n_repeats):
        inner_splits = _study_inner_folds(
            study_ids,
            n_inner_folds=n_inner_folds,
            seed=seed + repeat,
        )
        for train_idx, _val_idx in inner_splits:
            x_tr = x_train[train_idx]
            y_tr = y_train[train_idx]
            if x_tr.shape[0] < 2:
                continue
            if task != "age" and len(np.unique(y_tr)) < 2:
                continue
            for alpha in ENET_ALPHA_GRID:
                for l1_ratio in ENET_L1_GRID:
                    if task == "age":
                        model = _enet_regressor_pipeline(alpha=alpha, l1_ratio=l1_ratio)
                    else:
                        model = _enet_classifier_pipeline(alpha=alpha, l1_ratio=l1_ratio)
                    model.fit(x_tr, y_tr)
                    coef = _coef_abs(model, n_cols)
                    if coef.size != n_cols:
                        raise ValueError(f"coef size {coef.size} != n_cols {n_cols}")
                    selected = coef > 1e-8
                    counts[selected] += 1
                    coef_sum += coef
                    n_runs += 1
    if n_runs == 0:
        picked = np.arange(min(max_seeds, n_cols), dtype=np.int64)
        return picked, {"n_runs": 0, "frequency": {}, "mean_abs_coef": {}}
    freq = counts.astype(np.float64) / float(n_runs)
    order = np.argsort(-freq, kind="stable")
    picked = order[freq[order] >= min_frequency]
    if picked.size == 0:
        picked = order[: min(max_seeds, n_cols)]
    else:
        picked = picked[: min(max_seeds, picked.size)]
    meta = {
        "n_runs": n_runs,
        "frequency": {int(i): float(freq[i]) for i in picked.tolist()},
        "mean_abs_coef": {int(i): float(coef_sum[i] / n_runs) for i in picked.tolist()},
    }
    return picked.astype(np.int64), meta


def expand_panel_columns(
    seed_cols: np.ndarray,
    assignment: CascadeAssignment,
) -> np.ndarray:
    """Expand seeds: gene siblings + same-region for non-gene seeds."""
    seeds = set(int(c) for c in np.asarray(seed_cols, dtype=np.int64).tolist())
    expanded = set(seeds)
    if assignment.edge_col_index.size == 0:
        return np.asarray(sorted(expanded), dtype=np.int64)
    col_to_regions: dict[int, set[int]] = {}
    for col, reg in zip(
        assignment.edge_col_index.tolist(),
        assignment.edge_region_index.tolist(),
        strict=True,
    ):
        col_to_regions.setdefault(int(col), set()).add(int(reg))
    gene_to_cols: dict[int, set[int]] = {}
    for col, reg in zip(
        assignment.edge_col_index.tolist(),
        assignment.edge_region_index.tolist(),
        strict=True,
    ):
        g = int(assignment.region_to_gene[int(reg)])
        if g >= 0:
            gene_to_cols.setdefault(g, set()).add(int(col))
    for col in list(seeds):
        for reg in col_to_regions.get(col, ()):
            g = int(assignment.region_to_gene[reg])
            if g >= 0:
                expanded.update(gene_to_cols.get(g, ()))
            else:
                rid = assignment.region_ids[reg]
                for c2, r2 in zip(
                    assignment.edge_col_index.tolist(),
                    assignment.edge_region_index.tolist(),
                    strict=True,
                ):
                    if assignment.region_ids[int(r2)] == rid:
                        expanded.add(int(c2))
    return np.asarray(sorted(expanded), dtype=np.int64)


def select_multitask_fold_panel(
    *,
    x_train: np.ndarray,
    age: np.ndarray,
    age_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
    study_ids: np.ndarray,
    assignment: CascadeAssignment,
    max_seeds: int = 10_000,
    per_task_quota: int | None = None,
    matrix_id: str | None = None,
    graph_id: str | None = None,
    graph_content_hash: str | None = None,
) -> dict[str, Any]:
    """Outer-train multitask stability selection + graph expansion (canonical artifact)."""
    quota = per_task_quota or max(1, max_seeds // 3)
    seed_union: set[int] = set()
    seed_by_task: dict[str, list[int]] = {}
    meta_by_task: dict[str, Any] = {}

    age_m = np.asarray(age_mask, dtype=bool)
    if age_m.any():
        cols, meta = stability_select_columns(
            x_train[age_m],
            age[age_m],
            study_ids=study_ids[age_m],
            task="age",
            max_seeds=quota,
        )
        seed_by_task["age"] = cols.tolist()
        meta_by_task["age"] = meta
        seed_union.update(int(c) for c in cols.tolist())

    sex_m = np.asarray(sex_mask, dtype=bool)
    if sex_m.any():
        cols, meta = stability_select_columns(
            x_train[sex_m],
            sex[sex_m],
            study_ids=study_ids[sex_m],
            task="sex",
            max_seeds=quota,
        )
        seed_by_task["sex"] = cols.tolist()
        meta_by_task["sex"] = meta
        seed_union.update(int(c) for c in cols.tolist())

    tissue_m = np.asarray(tissue_mask, dtype=bool)
    if tissue_m.any():
        cols, meta = stability_select_columns(
            x_train[tissue_m],
            tissue[tissue_m],
            study_ids=study_ids[tissue_m],
            task="tissue",
            max_seeds=quota,
        )
        seed_by_task["tissue"] = cols.tolist()
        meta_by_task["tissue"] = meta
        seed_union.update(int(c) for c in cols.tolist())

    if not seed_union:
        raise ValueError("multitask panel selection found no labeled training samples")

    seeds = np.asarray(sorted(seed_union), dtype=np.int64)
    if seeds.size > max_seeds:
        seeds = seeds[:max_seeds]
    expanded = expand_panel_columns(seeds, assignment)
    train_studies = sorted({str(s) for s in study_ids.tolist()})
    imputer_hash = hashlib.sha256(b"SimpleImputer:median").hexdigest()[:16]
    scaler_hash = hashlib.sha256(b"StandardScaler:with_mean").hexdigest()[:16]
    return {
        "seed_cols_by_task": seed_by_task,
        "selection_meta_by_task": meta_by_task,
        "seed_cols": seeds.tolist(),
        "panel_cols": expanded.tolist(),
        "n_seed": int(seeds.size),
        "n_panel": int(expanded.size),
        "matrix_id": matrix_id,
        "graph_id": graph_id,
        "graph_content_hash": graph_content_hash,
        "train_study_ids": train_studies,
        "n_train_samples": int(x_train.shape[0]),
        "normalizer_hashes": {
            "imputer": imputer_hash,
            "scaler": scaler_hash,
        },
        "selector": "study_grouped_multitask_enet_stability",
    }


def select_fold_panel(
    *,
    x_train: np.ndarray,
    y_train: np.ndarray,
    assignment: CascadeAssignment,
    max_seeds: int = 10_000,
    study_ids: np.ndarray | None = None,
) -> dict[str, Any]:
    """Backward-compatible tissue-only wrapper."""
    cols, meta = stability_select_columns(
        x_train,
        y_train,
        study_ids=study_ids,
        task="tissue",
        max_seeds=max_seeds,
    )
    expanded = expand_panel_columns(cols, assignment)
    return {
        "seed_cols": cols.tolist(),
        "panel_cols": expanded.tolist(),
        "n_seed": int(cols.size),
        "n_panel": int(expanded.size),
        "selection_meta": meta,
    }
