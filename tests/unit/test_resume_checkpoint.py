"""§3.6 atomic checkpoint writes + fail-closed / faithful resume (CPU)."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from uuid import uuid4

import pytest
import torch

from mbs.paths import DataPaths
from mbs.training.loop import train_flat_baseline
from mbs.training.run_artifacts import (
    ResumeError,
    load_resume_checkpoint,
    save_checkpoint,
    sha256_file,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _point_env(monkeypatch: pytest.MonkeyPatch, scratch: Path) -> None:
    repo = _repo_root()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(scratch / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(scratch / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(scratch / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(scratch / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(scratch / "docker"))


@pytest.fixture
def isolated_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    scratch_base = _repo_root() / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"resume-{uuid4().hex}"
    workspace.mkdir()
    _point_env(monkeypatch, workspace)
    return workspace


def _tiny_payload(**overrides: Any) -> dict[str, Any]:
    base = {
        "model_state": {"w": torch.tensor([1.0])},
        "head_state": {"b": torch.tensor([0.0])},
        "optimizer_state": {"state": {}, "param_groups": []},
        "epoch": 1,
        "metrics": {"loss": 1.0},
        "config_hash": "abc",
        "resume_state": {
            "best_val": 1.0,
            "best_rank": None,
            "best_epoch": 1,
            "stale": 0,
            "history": [],
            "val_history": [],
            "n_samples_seen": 0,
            "n_optimizer_steps": 0,
            "batch_size": 1,
            "batch_token_budget": None,
            "torch_rng_state": torch.get_rng_state(),
            "numpy_rng_state": __import__("numpy").random.get_state(),
            "cuda_rng_state_all": None,
        },
    }
    base.update(overrides)
    return base


def test_atomic_checkpoint_preserves_previous_on_failed_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "last.pt"
    save_checkpoint(
        path,
        model_state={"w": torch.tensor([1.0])},
        head_state={"b": torch.tensor([0.0])},
        optimizer_state={"state": {}, "param_groups": []},
        epoch=1,
        metrics={"loss": 1.0},
        config_hash="abc",
        resume_state={"best_epoch": 1, "stale": 0},
    )
    digest = sha256_file(path)

    def _boom(*_args: Any, **_kwargs: Any) -> None:
        raise OSError("simulated disk full")

    monkeypatch.setattr(torch, "save", _boom)
    with pytest.raises(OSError, match="disk full"):
        save_checkpoint(
            path,
            model_state={"w": torch.tensor([2.0])},
            head_state={"b": torch.tensor([1.0])},
            optimizer_state={"state": {}, "param_groups": []},
            epoch=2,
            metrics={"loss": 0.5},
            config_hash="abc",
            resume_state={"best_epoch": 2, "stale": 0},
        )
    assert sha256_file(path) == digest
    assert not any(tmp_path.glob(".last.pt.*.tmp"))


def test_load_resume_rejects_truncated_checkpoint(tmp_path: Path) -> None:
    path = tmp_path / "last.pt"
    path.write_bytes(b"not-a-real-torch-checkpoint")
    with pytest.raises(ResumeError, match="cannot load"):
        load_resume_checkpoint(path, expected_config_hash="abc", max_epochs=10)


def test_load_resume_rejects_config_hash_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "last.pt"
    payload = _tiny_payload(config_hash="expected")
    torch.save(payload, path)
    with pytest.raises(ResumeError, match="config_hash mismatch"):
        load_resume_checkpoint(path, expected_config_hash="other", max_epochs=10)


def test_load_resume_rejects_missing_resume_state(tmp_path: Path) -> None:
    path = tmp_path / "last.pt"
    payload = _tiny_payload()
    del payload["resume_state"]
    torch.save(payload, path)
    with pytest.raises(ResumeError, match="no resume_state"):
        load_resume_checkpoint(path, expected_config_hash="abc", max_epochs=10)


def test_load_resume_rejects_finished_run(tmp_path: Path) -> None:
    path = tmp_path / "last.pt"
    torch.save(_tiny_payload(epoch=4), path)
    with pytest.raises(ResumeError, match="nothing to resume"):
        load_resume_checkpoint(path, expected_config_hash="abc", max_epochs=4)


def _overfit_cfg(*, max_epochs: int) -> dict[str, Any]:
    return {
        "experiment": {"name": "unit_resume", "stage": 0, "seed": 0},
        "model": {
            "phi_layers": 2,
            "phi_hidden_dimension": 32,
            "rho_layers": 2,
            "rho_hidden_dimension": 16,
            "pooling": "max",
            "neutral_score": 0.5,
            "dropout": 0.0,
        },
        "training": {
            "optimizer": "adam",
            "learning_rate": 1e-3,
            # Low LR so the overfit early-exit (acc≥0.999) does not fire mid-test.
            "weight_decay": 0.0,
            "max_epochs": max_epochs,
            "early_stopping_patience": 100,
            "gradient_clip_norm": 2.0,
            "precision": "fp32",
            "require_cuda": False,
            "resume": "auto",
            "batch_size": 4,
        },
        "heads": {"tissue": {"enabled": True}},
        "logging": {"tensorboard": False, "auto_tensorboard": False},
    }


def test_resume_matches_uninterrupted_cpu_run(isolated_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    n_total = 4
    n_partial = 2
    cfg = _overfit_cfg(max_epochs=n_total)

    straight = train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        run_id="resume-straight",
        device_str="cpu",
        overfit_fixture=True,
        max_epochs=n_total,
    )
    straight_last = torch.load(
        straight.checkpoint_dir / "last.pt",
        map_location="cpu",
        weights_only=False,
    )

    train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        run_id="resume-split",
        device_str="cpu",
        overfit_fixture=True,
        max_epochs=n_partial,
    )
    resumed = train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        run_id="resume-split",
        device_str="cpu",
        overfit_fixture=True,
        max_epochs=n_total,
    )
    resumed_last = torch.load(
        resumed.checkpoint_dir / "last.pt",
        map_location="cpu",
        weights_only=False,
    )

    assert straight_last["epoch"] == resumed_last["epoch"] == n_total
    assert "resume_state" in resumed_last
    for key in ("model_state", "head_state"):
        for t_name, t_a in straight_last[key].items():
            t_b = resumed_last[key][t_name]
            assert torch.equal(t_a, t_b), f"{key}.{t_name} diverged"
    assert straight_last["resume_state"]["history"] == resumed_last["resume_state"]["history"]
    assert straight_last["resume_state"]["best_epoch"] == resumed_last["resume_state"]["best_epoch"]
    assert straight_last["resume_state"]["stale"] == resumed_last["resume_state"]["stale"]
