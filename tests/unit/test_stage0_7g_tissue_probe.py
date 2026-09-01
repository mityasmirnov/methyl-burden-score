"""Milestone 7G cascade tissue probe unit checks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mbs.training.cascade_loop as cascade_loop
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import (
    make_synthetic_cascade_tables,
    train_cascade_on_arrays,
)
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.transparent_baselines import fit_linear_multitask, run_mean_baseline


def test_balanced_logistic_fusion_runs() -> None:
    rng = np.random.default_rng(1)
    n_tr, d = 60, 12
    x_tr = rng.normal(size=(n_tr, d))
    tissue_tr = rng.integers(0, 4, size=n_tr)
    models = fit_linear_multitask(
        x_tr,
        age=None,
        age_mask=None,
        tissue=tissue_tr,
        tissue_mask=np.ones(n_tr, dtype=bool),
        sex=None,
        sex_mask=None,
        tissue_solver="balanced_logistic",
        fusion_pca_components=4,
    )
    assert "tissue" in models
    assert "_pca" in models


def test_evaluate_late_fusion_with_fusion_kwargs() -> None:
    rng = np.random.default_rng(2)
    n_tr, n_te, d = 50, 25, 8
    x_tr = rng.normal(size=(n_tr, d))
    x_te = rng.normal(size=(n_te, d))
    tissue_tr = (x_tr[:, 0] > 0).astype(np.int64)
    tissue_te = (x_te[:, 0] > 0).astype(np.int64)
    out = evaluate_late_fusion(
        scores_train=x_tr,
        scores_test=x_te,
        age_train=None,
        age_mask_train=None,
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_tr, dtype=bool),
        sex_train=None,
        sex_mask_train=None,
        age_test=None,
        age_mask_test=None,
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_te, dtype=bool),
        sex_test=None,
        sex_mask_test=None,
        fusion={"tissue_solver": "balanced_logistic", "pca_components": 4},
    )
    assert "metrics" in out
    assert out["metrics"]["tissue"]["macro_f1"] is not None


def test_cascade_tissue_loss_weight_smoke(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=3)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(3, (n * 2) // 3), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()
    out = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=tmp_path / "p2fold",
        max_epochs=2,
        seed=0,
        device_str="cpu",
        tissue_loss_weight=3.0,
        age_loss_weight=0.3,
        sex_loss_weight=1.0,
    )
    assert "metrics" in out
    assert (tmp_path / "p2fold" / "scores" / "score_manifest.json").is_file()


def test_cascade_mean_pooling_smoke(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=5)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(3, (n * 2) // 3), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()
    out = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=tmp_path / "p4fold",
        max_epochs=2,
        seed=0,
        device_str="cpu",
        cpg_pool="mean",
        region_pool="mean",
        tissue_loss_weight=3.0,
        age_loss_weight=0.3,
    )
    assert out.get("pooling") == {
        "cpg_to_region": "mean",
        "region_to_gene": "mean",
    }


def test_cascade_early_stop_smoke(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=6)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, n // 2), dtype=np.int64)
    val_idx = np.arange(train_idx[-1] + 1, min(train_idx[-1] + 3, n), dtype=np.int64)
    test_idx = np.arange(val_idx[-1] + 1, n, dtype=np.int64) if val_idx[-1] + 1 < n else val_idx
    out = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=tmp_path / "p5fold",
        max_epochs=20,
        seed=0,
        device_str="cpu",
        val_idx=val_idx,
        early_stopping_patience=2,
        tissue_loss_weight=3.0,
    )
    sel = out.get("checkpoint_selection") or {}
    assert sel.get("has_validation") is True
    best = sel.get("best_epoch")
    assert best is not None and int(best) <= 20


def test_region_mean_baseline_smoke() -> None:
    rng = np.random.default_rng(4)
    n_tr, n_te, g = 40, 20, 6
    x_tr = rng.random((n_tr, g)).astype(np.float32)
    x_te = rng.random((n_te, g)).astype(np.float32)
    tissue_tr = rng.integers(0, 3, size=n_tr)
    tissue_te = rng.integers(0, 3, size=n_te)
    out = run_mean_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=None,
        age_mask_train=None,
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_tr, dtype=bool),
        sex_train=None,
        sex_mask_train=None,
        age_test=None,
        age_mask_test=None,
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_te, dtype=bool),
        sex_test=None,
        sex_mask_test=None,
        kind="region",
    )
    assert out["kind"] == "region"
    assert "tissue" in out["metrics"]


def test_mean_pooling_and_tissue_early_stop_are_wired(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tables = make_synthetic_cascade_tables(seed=8)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, 6, dtype=np.int64)
    val_idx = np.arange(6, 9, dtype=np.int64)
    test_idx = np.arange(9, n, dtype=np.int64)

    def _constant_validation(*args: object, **kwargs: object) -> dict[str, float]:
        del args, kwargs
        return {"tissue_macro_f1": 0.25, "age_mae": 10.0}

    monkeypatch.setattr(
        cascade_loop,
        "_evaluate_cascade_validation",
        _constant_validation,
    )
    out = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        val_idx=val_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=tmp_path / "phase2",
        max_epochs=8,
        seed=0,
        device_str="cpu",
        cpg_pool="mean",
        region_pool="mean",
        early_stopping_patience=2,
        early_stopping_min_delta=0.0,
    )
    selection = out["checkpoint_selection"]
    assert selection["epochs_completed"] == 3
    assert selection["early_stopping"]["stopped_early"] is True
    assert selection["early_stopping"]["stop_epoch"] == 3
    assert selection["best_epoch"] == 1
    assert out["pooling"] == {
        "cpg_to_region": "mean",
        "region_to_gene": "mean",
    }
