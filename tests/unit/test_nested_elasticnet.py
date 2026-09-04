"""Unit tests for nested elastic-net (train-fold scaler + inner-val α/l1)."""

from __future__ import annotations

import numpy as np

from mbs.training.transparent_baselines import (
    NESTED_ENET_ALPHA_GRID,
    run_nested_elasticnet_multitask,
)


def test_nested_enet_scaler_fit_on_train_only() -> None:
    """StandardScaler mean must equal train mean, not train+test."""
    rng = np.random.default_rng(0)
    n_train, n_test, n_feat = 60, 20, 8
    x_train = rng.normal(loc=0.0, scale=1.0, size=(n_train, n_feat))
    # Test shifted — leaked scaler would pull mean away from train.
    x_test = rng.normal(loc=50.0, scale=1.0, size=(n_test, n_feat))
    studies = np.asarray([f"s{i % 6}" for i in range(n_train)], dtype=object)
    age_train = 40.0 + x_train[:, 0] * 5.0 + rng.normal(0, 0.5, size=n_train)
    age_test = 40.0 + x_test[:, 0] * 5.0 + rng.normal(0, 0.5, size=n_test)
    tissue_train = (x_train[:, 1] > 0).astype(np.int64)
    tissue_test = (x_test[:, 1] > 0).astype(np.int64)
    # Ensure both classes in train.
    tissue_train[:2] = [0, 1]
    sex_train = (x_train[:, 2] > 0).astype(np.int64)
    sex_test = (x_test[:, 2] > 0).astype(np.int64)
    sex_train[:2] = [0, 1]

    out = run_nested_elasticnet_multitask(
        x_train=x_train,
        x_test=x_test,
        age_train=age_train,
        age_mask_train=np.ones(n_train, dtype=bool),
        tissue_train=tissue_train,
        tissue_mask_train=np.ones(n_train, dtype=bool),
        sex_train=sex_train,
        sex_mask_train=np.ones(n_train, dtype=bool),
        age_test=age_test,
        age_mask_test=np.ones(n_test, dtype=bool),
        tissue_test=tissue_test,
        tissue_mask_test=np.ones(n_test, dtype=bool),
        sex_test=sex_test,
        sex_mask_test=np.ones(n_test, dtype=bool),
        study_ids_train=studies,
        study_ids_test=np.asarray([f"t{i}" for i in range(n_test)], dtype=object),
        alpha_grid=(1e-3, 0.1),
        l1_grid=(0.5,),
        seed=7,
    )
    assert out["kind"] == "elasticnet_nested"
    assert out["scaler"]["fit_on"] == "outer_train_only"
    recorded = np.asarray(out["scaler"]["train_feature_mean"], dtype=np.float64)
    expected = x_train.mean(axis=0)
    leaked = np.concatenate([x_train, x_test], axis=0).mean(axis=0)
    np.testing.assert_allclose(recorded, expected, rtol=1e-6, atol=1e-6)
    # Must not match the train+test mean (test is heavily shifted).
    assert float(np.linalg.norm(recorded - leaked)) > 1.0
    age_sel = out["selected"]["age"]
    assert age_sel["alpha"] in (1e-3, 0.1)
    assert age_sel["l1_ratio"] == 0.5
    assert age_sel["n_inner_train"] >= 1
    assert age_sel["n_inner_val"] >= 1


def test_nested_enet_selected_alpha_from_grid() -> None:
    rng = np.random.default_rng(1)
    n, d = 40, 4
    x = rng.normal(size=(n, d))
    age = x[:, 0] * 3.0
    studies = np.asarray([f"g{i % 4}" for i in range(n)], dtype=object)
    out = run_nested_elasticnet_multitask(
        x_train=x[:30],
        x_test=x[30:],
        age_train=age[:30],
        age_mask_train=np.ones(30, dtype=bool),
        tissue_train=None,
        tissue_mask_train=None,
        sex_train=None,
        sex_mask_train=None,
        age_test=age[30:],
        age_mask_test=np.ones(10, dtype=bool),
        tissue_test=None,
        tissue_mask_test=None,
        sex_test=None,
        sex_mask_test=None,
        study_ids_train=studies[:30],
        study_ids_test=studies[30:],
        seed=3,
    )
    assert out["selected"]["age"]["alpha"] in NESTED_ENET_ALPHA_GRID
    assert out["selected"]["age"]["l1_ratio"] in (0.25, 0.5, 0.75)
    assert "age" in (out.get("metrics") or {})
