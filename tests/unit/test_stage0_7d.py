"""Milestone 7D Level-1 fold-fitted MAD robust-z."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pytest

from mbs.evaluation.metrics import binary_auroc_auprc, expected_calibration_error
from mbs.paths import DataPaths
from mbs.training.controls import apply_feature_control
from mbs.training.direct_cpg import (
    direct_cpg_design_matrix,
    fit_direct_elasticnet,
)
from mbs.training.features import beta_to_m_value, cpg_input_dim
from mbs.training.level1_norm import (
    MAD_SCALE,
    apply_level1_robust_z,
    fit_level1_robust_z,
    load_level1,
    persist_level1,
    resolve_level1_config,
)
from mbs.training.loop import train_flat_baseline


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _point_env(monkeypatch: pytest.MonkeyPatch, workspace: Path) -> None:
    monkeypatch.setenv("MBS_ROOT", str(workspace))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(workspace / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(workspace / "cache"))


@pytest.fixture
def isolated_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    scratch_base = _repo_root() / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"7d-{uuid4().hex}"
    workspace.mkdir()
    _point_env(monkeypatch, workspace)
    return workspace


def _base_cfg(*, robust_deviation: bool = False, seed: int = 0) -> dict[str, Any]:
    return {
        "experiment": {"name": "unit_7d", "stage": 0, "seed": seed},
        "features": {
            "methylation": {
                "beta": True,
                "m_value": True,
                "robust_deviation": robust_deviation,
                "epsilon": 0.001,
                "sigma_min": 1e-6,
            }
        },
        "model": {
            "phi_layers": 2,
            "phi_hidden_dimension": 16,
            "rho_layers": 2,
            "rho_hidden_dimension": 8,
            "pooling": "max",
            "neutral_score": 0.5,
            "dropout": 0.0,
        },
        "training": {
            "optimizer": "adam",
            "learning_rate": 0.05,
            "weight_decay": 0.0,
            "max_epochs": 2,
            "early_stopping_patience": 5,
            "gradient_clip_norm": 2.0,
            "precision": "fp32",
            "require_cuda": False,
        },
        "heads": {"tissue": {"enabled": True}},
        "logging": {"tensorboard": False, "auto_tensorboard": False},
    }


def test_level1_two_point_formula() -> None:
    # Two finite values -> median midpoints; MAD = |x-med|; sigma = 1.4826*MAD.
    m = np.array([[1.0, 10.0], [3.0, 10.0]], dtype=np.float64)
    params = fit_level1_robust_z(m, sigma_min=1e-6)
    assert params.estimated.tolist() == [True, True]
    assert abs(params.mu[0] - 2.0) < 1e-12
    assert abs(params.mu[1] - 10.0) < 1e-12
    mad0 = 1.0  # median(|1-2|, |3-2|)
    assert abs(params.sigma[0] - MAD_SCALE * mad0) < 1e-12
    # Constant column -> MAD=0 -> sigma_min floor
    assert abs(params.sigma[1] - 1e-6) < 1e-15
    z, present = apply_level1_robust_z(np.array([1.0, 10.0]), params)
    assert present.tolist() == [True, True]
    assert abs(float(z[0]) - (1.0 - 2.0) / params.sigma[0]) < 1e-6
    assert abs(float(z[1])) < 1e-6


def test_novel_locus_zero_z_kept() -> None:
    m = np.array([[1.0, np.nan], [2.0, np.nan]], dtype=np.float64)
    params = fit_level1_robust_z(m)
    assert params.estimated.tolist() == [True, False]
    z, present = apply_level1_robust_z(np.array([1.5, 9.0]), params)
    assert present.tolist() == [True, False]
    assert float(z[1]) == 0.0
    # Row kept (length 2)
    assert z.shape == (2,)


def test_no_leakage_val_only_locus() -> None:
    # Locus 1 only observed on "val" row → must not be estimated from train.
    train_m = np.array([[1.0, np.nan], [3.0, np.nan]], dtype=np.float64)
    params = fit_level1_robust_z(train_m)
    assert params.estimated[0]
    assert not params.estimated[1]
    z_val, present = apply_level1_robust_z(np.array([2.0, 100.0]), params)
    assert present[0]
    assert not present[1]
    assert float(z_val[1]) == 0.0


def test_persist_roundtrip_and_tamper(tmp_path: Path) -> None:
    m = np.array([[0.0, 1.0], [2.0, 3.0]], dtype=np.float64)
    params = fit_level1_robust_z(m, fold_id="f0", run_id="r0")
    manifest = persist_level1(tmp_path, params)
    loaded, loaded_manifest = load_level1(tmp_path)
    assert np.allclose(loaded.mu, params.mu)
    assert loaded_manifest["mu_sha256"] == manifest["mu_sha256"]
    # Tamper mu.npy
    np.save(tmp_path / "fold_norm" / "mu.npy", params.mu + 1.0)
    with pytest.raises(ValueError, match="sha256 mismatch"):
        load_level1(tmp_path)


def test_channel_a_vs_b_identical_holdout_folds(isolated_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    seed = 7
    cfg_a = _base_cfg(robust_deviation=False, seed=seed)
    cfg_b = _base_cfg(robust_deviation=True, seed=seed)
    result_a = train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg_a,
        run_id="7d-channel-a",
        device_str="cpu",
        study_holdout_fixture=True,
        max_epochs=1,
    )
    result_b = train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg_b,
        run_id="7d-channel-b",
        device_str="cpu",
        study_holdout_fixture=True,
        max_epochs=1,
    )
    split_a = json.loads((result_a.run_dir / "split.json").read_text(encoding="utf-8"))
    split_b = json.loads((result_b.run_dir / "split.json").read_text(encoding="utf-8"))
    assert split_a["train_sample_ids"] == split_b["train_sample_ids"]
    assert split_a["validation_sample_ids"] == split_b["validation_sample_ids"]
    dim_a = cpg_input_dim(4, include_m_value=True, include_robust_z=False)
    dim_b = cpg_input_dim(4, include_m_value=True, include_robust_z=True)
    assert dim_b == dim_a + 2
    resolved_b = (result_b.run_dir / "resolved_config.yaml").read_text(encoding="utf-8")
    assert f"input_dim: {dim_b}" in resolved_b or f"input_dim: {dim_b}\n" in resolved_b
    # Confirm via metrics level1 block and fold_norm artifacts
    assert result_b.metrics["level1_normalization"]["enabled"] is True
    assert result_a.metrics["level1_normalization"]["enabled"] is False
    assert (result_b.run_dir / "fold_norm" / "manifest.json").is_file()
    assert not (result_a.run_dir / "fold_norm").exists()
    # Never write under canonical matrices
    canonical = paths.data_root / "canonical" / "matrices"
    if canonical.exists():
        assert list(canonical.rglob("*")) == []
    else:
        # Isolated workspace should not create Hub matrices either
        assert not canonical.exists()


def test_controls_with_robust_z_channels() -> None:
    # beta, M, z, static(2), static_present, norm_present → dim 7
    feats = np.ones((3, 7), dtype=np.float32)
    static_only = apply_feature_control(
        feats, mode="static_only", include_m_value=True, include_robust_z=True
    )
    assert np.allclose(static_only[:, :3], 0.0)  # beta, M, z zeroed
    coverage = apply_feature_control(
        feats, mode="coverage_only", include_m_value=True, include_robust_z=True
    )
    assert np.allclose(coverage[:, :-2], 0.0)
    assert np.allclose(coverage[:, -2:], 1.0)  # present flags kept


def test_direct_elasticnet_on_level1_z() -> None:
    rng = np.random.default_rng(0)
    betas = rng.uniform(0.1, 0.9, size=(30, 6))
    m = beta_to_m_value(betas)
    obs = np.ones_like(m, dtype=bool)
    z, params = direct_cpg_design_matrix(m, obs, use_level1=True, sigma_min=1e-6)
    assert params is not None
    assert params.n_estimated == 6
    y = z[:, 0] + 0.01 * rng.normal(size=30)
    studies = np.array(["A"] * 15 + ["B"] * 15)
    fitted = fit_direct_elasticnet(z, obs, y, studies, min_studies=2)
    assert fitted["n_loci"] == 6
    assert fitted["weights"].shape[0] == 6


def test_config_robust_requires_m_value() -> None:
    with pytest.raises(ValueError, match="m_value"):
        resolve_level1_config(
            {"features": {"methylation": {"robust_deviation": True, "m_value": False}}}
        )


def test_config_level2_level3_not_implemented() -> None:
    with pytest.raises(NotImplementedError, match="Level-2"):
        resolve_level1_config({"features": {"methylation": {"level2_probe_adapter": True}}})
    with pytest.raises(NotImplementedError, match="Level-3"):
        resolve_level1_config({"features": {"methylation": {"level3_masked_ae": True}}})


def test_binary_metrics_helpers_available() -> None:
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    out = binary_auroc_auprc(y, s)
    assert "auroc" in out and "auprc" in out
    ece = expected_calibration_error(y, s)
    assert "ece" in ece
