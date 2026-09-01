"""Fold-safe probe panel selection for 7G′ Stage B (C-mvalue-enetS)."""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mbs.training.cascade_assign import CascadeAssignment


def _enet_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler(with_mean=True)),
            (
                "sgd",
                SGDClassifier(
                    loss="log_loss",
                    penalty="elasticnet",
                    alpha=1e-4,
                    l1_ratio=0.5,
                    max_iter=40,
                    tol=1e-3,
                    random_state=42,
                    n_jobs=-1,
                ),
            ),
        ]
    )


def stability_select_columns(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    max_seeds: int = 10_000,
    n_inner_folds: int = 3,
    min_frequency: float = 0.34,
) -> np.ndarray:
    """Repeated inner-CV elastic-net; rank columns by selection frequency."""
    if x_train.shape[0] < n_inner_folds * 2:
        n_inner_folds = max(2, min(n_inner_folds, x_train.shape[0] // 2))
    n_cols = x_train.shape[1]
    counts = np.zeros(n_cols, dtype=np.int64)
    n_runs = 0
    skf = StratifiedKFold(n_splits=n_inner_folds, shuffle=True, random_state=42)
    imp = SimpleImputer(strategy="median")
    x_imp = imp.fit_transform(x_train)
    for train_idx, _ in skf.split(x_imp, y_train):
        model = _enet_pipeline()
        model.fit(x_imp[train_idx], y_train[train_idx])
        raw_coef = model.named_steps["sgd"].coef_
        if raw_coef.ndim == 2:
            coef = np.max(np.abs(raw_coef), axis=0)
        else:
            coef = np.abs(raw_coef.ravel())
        if coef.size != n_cols:
            raise ValueError(f"coef size {coef.size} != n_cols {n_cols}")
        selected = coef > 1e-8
        counts[selected] += 1
        n_runs += 1
    if n_runs == 0:
        return np.arange(min(max_seeds, n_cols), dtype=np.int64)
    freq = counts.astype(np.float64) / float(n_runs)
    order = np.argsort(-freq, kind="stable")
    picked = order[freq[order] >= min_frequency]
    if picked.size == 0:
        picked = order[: min(max_seeds, n_cols)]
    else:
        picked = picked[: min(max_seeds, picked.size)]
    return picked.astype(np.int64)


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


def select_fold_panel(
    *,
    x_train: np.ndarray,
    y_train: np.ndarray,
    assignment: CascadeAssignment,
    max_seeds: int = 10_000,
) -> dict[str, Any]:
    """Outer-train stability selection + graph expansion."""
    seeds = stability_select_columns(x_train, y_train, max_seeds=max_seeds)
    expanded = expand_panel_columns(seeds, assignment)
    return {
        "seed_cols": seeds.tolist(),
        "panel_cols": expanded.tolist(),
        "n_seed": int(seeds.size),
        "n_panel": int(expanded.size),
    }
