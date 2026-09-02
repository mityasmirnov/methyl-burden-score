"""Unit checks for fold-safe panel selection."""

from __future__ import annotations

import numpy as np

from mbs.training.cascade_loop import make_synthetic_cascade_tables
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.fold_safe_panel import _study_inner_folds, select_multitask_fold_panel


def test_study_inner_folds_no_study_leakage() -> None:
    study_ids = np.array(["A", "A", "B", "B", "C", "C"], dtype=object)
    splits = _study_inner_folds(study_ids, n_inner_folds=3, seed=0)
    assert splits
    for train_idx, val_idx in splits:
        train_studies = set(study_ids[train_idx].tolist())
        val_studies = set(study_ids[val_idx].tolist())
        assert train_studies.isdisjoint(val_studies)


def test_select_multitask_fold_panel_smoke() -> None:
    tables = make_synthetic_cascade_tables(seed=1)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = tables["betas"].shape[0]
    train_idx = np.arange(0, n - 2, dtype=np.int64)
    x = tables["betas"][train_idx].astype(np.float32)
    studies = tables["study_ids"][train_idx]
    panel = select_multitask_fold_panel(
        x_train=x,
        age=tables["ages"][train_idx],
        age_mask=np.ones(train_idx.size, dtype=bool),
        sex=tables["sex"][train_idx],
        sex_mask=np.ones(train_idx.size, dtype=bool),
        tissue=tables["tissue"][train_idx],
        tissue_mask=np.ones(train_idx.size, dtype=bool),
        study_ids=studies,
        assignment=assignment,
        max_seeds=20,
    )
    assert panel["n_panel"] >= panel["n_seed"] > 0
    assert "seed_cols_by_task" in panel
    assert panel["selector"] == "study_grouped_multitask_enet_stability"
