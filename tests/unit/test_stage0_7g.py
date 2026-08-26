"""Milestone 7G unit checks: masks, ROC, skip-if-done, ranking hygiene."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import (
    make_synthetic_cascade_tables,
    train_cascade_on_arrays,
)
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.transparent_baselines import run_mean_baseline

# Keep in sync with scripts/write_7g_methylation_eval_report.py RANKING_ARMS.
RANKING_ARMS = (
    "N-cascade-l1",
    "T-mean-gene",
    "T-mean-region",
    "T-enet",
    "C-mvalue-ridge",
    "C-mvalue-enet",
    "C-mvalue-hgb",
    "C-mvalue-sva",
)


def test_masked_loss_skips_unlabeled_sex(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=7)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, (n * 2) // 3), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()
    sex_mask = np.zeros(n, dtype=bool)
    age_mask = np.ones(n, dtype=bool)
    tissue_mask = np.ones(n, dtype=bool)
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
        out_dir=tmp_path / "fold0",
        max_epochs=2,
        seed=0,
        device_str="cpu",
        age_mask=age_mask,
        tissue_mask=tissue_mask,
        sex_mask=sex_mask,
        cpg_hidden_dim=8,
        region_hidden_dim=4,
        dropout=0.0,
    )
    assert "metrics" in out
    assert (tmp_path / "fold0" / "scores" / "score_manifest.json").is_file()
    # Unlabeled sex → no sex metrics block.
    assert "sex" not in out["metrics"]


def test_fusion_emits_sex_auroc_and_tissue_ovr() -> None:
    rng = np.random.default_rng(0)
    n_tr, n_te, d = 40, 20, 6
    x_tr = rng.normal(size=(n_tr, d)).astype(np.float64)
    x_te = rng.normal(size=(n_te, d)).astype(np.float64)
    sex_tr = (x_tr[:, 0] > 0).astype(np.int64)
    sex_te = (x_te[:, 0] > 0).astype(np.int64)
    tissue_tr = (x_tr[:, 1] > 0).astype(np.int64)
    tissue_te = (x_te[:, 1] > 0).astype(np.int64)
    age_tr = x_tr[:, 2] * 10 + 40
    age_te = x_te[:, 2] * 10 + 40
    fused = evaluate_late_fusion(
        scores_train=x_tr,
        scores_test=x_te,
        age_train=age_tr,
        age_mask_train=np.ones(n_tr, dtype=bool),
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_tr, dtype=bool),
        sex_train=sex_tr,
        sex_mask_train=np.ones(n_tr, dtype=bool),
        age_test=age_te,
        age_mask_test=np.ones(n_te, dtype=bool),
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_te, dtype=bool),
        sex_test=sex_te,
        sex_mask_test=np.ones(n_te, dtype=bool),
        tissue_class_names=["A", "B"],
    )
    sex = fused["metrics"]["sex"]
    assert "auroc" in sex
    assert sex["auroc"] > 0.7
    assert "fpr" in sex and "tpr" in sex
    assert "tissue_roc" in fused["metrics"]


def test_skip_if_done_reuses_metrics(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=3)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, n // 2), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()
    kwargs = dict(
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
        out_dir=tmp_path / "fold",
        max_epochs=1,
        seed=1,
        device_str="cpu",
        cpg_hidden_dim=8,
        region_hidden_dim=4,
        dropout=0.0,
    )
    first = train_cascade_on_arrays(**kwargs, skip_if_done=False)
    assert first.get("skipped") is False
    second = train_cascade_on_arrays(**kwargs, skip_if_done=True)
    assert second.get("skipped") is True
    assert second["score_dir"] == first["score_dir"]


def test_t_mean_region_is_distinct_arm_name() -> None:
    rng = np.random.default_rng(1)
    x_tr = rng.normal(size=(20, 4)).astype(np.float32)
    x_te = rng.normal(size=(10, 4)).astype(np.float32)
    out = run_mean_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=np.linspace(20, 60, 20),
        age_mask_train=np.ones(20, dtype=bool),
        tissue_train=np.array([0, 1] * 10),
        tissue_mask_train=np.ones(20, dtype=bool),
        sex_train=np.array([0, 1] * 10),
        sex_mask_train=np.ones(20, dtype=bool),
        age_test=np.linspace(20, 60, 10),
        age_mask_test=np.ones(10, dtype=bool),
        tissue_test=np.array([0, 1] * 5),
        tissue_mask_test=np.ones(10, dtype=bool),
        sex_test=np.array([0, 1] * 5),
        sex_mask_test=np.ones(10, dtype=bool),
        kind="region",
    )
    assert out["kind"] == "region"


def test_ranking_arms_exclude_metadata_only() -> None:
    assert "C-metadata" not in RANKING_ARMS
    assert "T-mean-region" in RANKING_ARMS
    assert "N-cascade-l1" in RANKING_ARMS
