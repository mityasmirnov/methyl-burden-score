"""Milestone 7E development CV unit tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mbs.evaluation.splits import assert_no_study_leakage, build_outer_study_grouped_folds
from mbs.models import FlatDeepSet, HierarchicalDeepSet, SharedMLP
from mbs.training.dev_cv import arm_run_id, freeze_outer_folds, inject_fold_into_config
from mbs.training.encoder_config import resolve_encoder
from mbs.training.late_fusion import concatenate_score_blocks, evaluate_late_fusion
from mbs.training.transparent_baselines import presence_aware_means, run_mean_baseline


def test_outer_folds_study_disjoint(tmp_path: Path) -> None:
    samples = [
        {"sample_id": f"s{i}", "study_id": f"GSE{i % 9}", "platform": "HM450"} for i in range(90)
    ]
    pack = freeze_outer_folds(samples, out_dir=tmp_path, n_folds=3, seed=7)
    assert pack["n_folds"] == 3
    assert (tmp_path / "folds.json").is_file()
    assert (tmp_path / "folds.sha256.json").is_file()
    test_studies: set[str] = set()
    for fold in pack["folds"]:
        assert_no_study_leakage(fold)
        assert set(fold["external_test_studies"]).isdisjoint(fold["train_studies"])
        assert set(fold["external_test_studies"]).isdisjoint(fold["validation_studies"])
        for sid in fold["external_test_studies"]:
            assert sid not in test_studies
            test_studies.add(sid)
    assert len(test_studies) == 9


def test_identical_folds_stable() -> None:
    samples = [
        {"sample_id": f"s{i}", "study_id": f"GSE{i % 6}", "platform": "HM450"} for i in range(60)
    ]
    a = build_outer_study_grouped_folds(samples, n_folds=3, seed=42)
    b = build_outer_study_grouped_folds(samples, n_folds=3, seed=42)
    assert a["folds"][0]["train_studies"] == b["folds"][0]["train_studies"]
    assert a["folds"][1]["external_test_studies"] == b["folds"][1]["external_test_studies"]


def test_matched_encoder_flat_hier_rho_activation() -> None:
    enc = resolve_encoder(
        {
            "encoder": {
                "activation": "gelu",
                "dropout": 0.1,
                "layer_norm": True,
                "cpg_hidden_dim": 64,
            }
        }
    )
    flat = FlatDeepSet(
        4,
        phi_hidden_dim=enc["cpg_hidden_dim"],
        dropout=enc["dropout"],
        activation=enc["activation"],
        layer_norm=enc["layer_norm"],
    )
    hier = HierarchicalDeepSet(
        4,
        5,
        cpg_hidden_dim=enc["cpg_hidden_dim"],
        dropout=enc["dropout"],
        activation=enc["activation"],
        layer_norm=enc["layer_norm"],
    )
    assert type(flat.phi) is SharedMLP
    assert type(hier.rho) is SharedMLP
    assert type(hier.residual_rho) is SharedMLP
    # Activation modules: SharedMLP stores activation name indirectly via layers;
    # gelu path uses nn.GELU — ensure rho is not stuck on leaky_relu default.
    assert any(type(m).__name__ == "GELU" for m in hier.rho.modules())
    assert any(type(m).__name__ == "GELU" for m in hier.cpg_encoder.modules())


def test_transparent_mean_baseline_smoke() -> None:
    rng = np.random.default_rng(0)
    n_train, n_test, n_loci, n_genes = 20, 8, 12, 4
    gene_index = np.array([i % n_genes for i in range(n_loci)], dtype=np.int64)
    betas_tr = rng.uniform(0.1, 0.9, size=(n_train, n_loci))
    obs_tr = np.ones((n_train, n_loci), dtype=bool)
    betas_te = rng.uniform(0.1, 0.9, size=(n_test, n_loci))
    obs_te = np.ones((n_test, n_loci), dtype=bool)
    x_tr = presence_aware_means(betas_tr, obs_tr, gene_index, n_groups=n_genes)
    x_te = presence_aware_means(betas_te, obs_te, gene_index, n_groups=n_genes)
    age_tr = 30.0 + x_tr.mean(axis=1)
    age_te = 30.0 + x_te.mean(axis=1)
    tissue_tr = (x_tr[:, 0] > 0.5).astype(np.int64)
    tissue_te = (x_te[:, 0] > 0.5).astype(np.int64)
    out = run_mean_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=age_tr,
        age_mask_train=np.ones(n_train, dtype=bool),
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_train, dtype=bool),
        sex_train=None,
        sex_mask_train=None,
        age_test=age_te,
        age_mask_test=np.ones(n_test, dtype=bool),
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_test, dtype=bool),
        sex_test=None,
        sex_mask_test=None,
        kind="gene",
    )
    assert "age" in out["metrics"]
    assert "rmse" in out["metrics"]["age"]


def test_late_fusion_shapes() -> None:
    rng = np.random.default_rng(1)
    a = rng.normal(size=(10, 3)).astype(np.float32)
    b = rng.normal(size=(10, 2)).astype(np.float32)
    scores = concatenate_score_blocks([a, b])
    assert scores.shape == (10, 5)
    age = rng.normal(size=10)
    tissue = rng.integers(0, 2, size=10)
    out = evaluate_late_fusion(
        scores_train=scores[:7],
        scores_test=scores[7:],
        age_train=age[:7],
        age_mask_train=np.ones(7, dtype=bool),
        tissue_train=tissue[:7],
        tissue_mask_train=np.ones(7, dtype=bool),
        sex_train=None,
        sex_mask_train=None,
        age_test=age[7:],
        age_mask_test=np.ones(3, dtype=bool),
        tissue_test=tissue[7:],
        tissue_mask_test=np.ones(3, dtype=bool),
        sex_test=None,
        sex_mask_test=None,
    )
    assert out["n_score_features"] == 5
    assert "age" in out["metrics"]


def test_inject_fold_disables_auto_and_reuse() -> None:
    fold = {
        "train_studies": ["A"],
        "validation_studies": ["B"],
        "external_test_studies": ["C"],
        "split_id": "hub-ats-7e-3fold-v1/fold-0",
        "outer_fold": 0,
    }
    cfg = inject_fold_into_config(
        {"pilot": {"auto_split": True, "reuse_flat_split": True}, "experiment": {}},
        fold,
        seed=99,
    )
    assert cfg["pilot"]["auto_split"] is False
    assert cfg["pilot"]["reuse_flat_split"] is False
    assert cfg["pilot"]["train_studies"] == ["A"]
    assert cfg["experiment"]["seed"] == 99
    assert arm_run_id("N-flat", 1, 0, tag="x") == "stage0-7e-N-flat-f1-r0-x"
