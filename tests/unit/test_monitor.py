"""Unit tests for the lightweight train-run monitor."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import pytest
import yaml
from typer.testing import CliRunner

from mbs.cli import app
from mbs.training.monitor import (
    EpochMetrics,
    collect_snapshot,
    estimate_eta,
    parse_metrics_row,
    read_metrics_jsonl,
    render_snapshot,
    resolve_max_epochs,
    validate_run_id,
)

runner = CliRunner()


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
def art_root(monkeypatch: pytest.MonkeyPatch) -> Path:
    # Keep under /data via project scratch (path policy).
    scratch_base = _repo_root() / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"monitor-{uuid4().hex}"
    workspace.mkdir(parents=True)
    _point_env(monkeypatch, workspace)
    root = workspace / "artifacts"
    (root / "runs").mkdir(parents=True)
    (root / "checkpoints").mkdir(parents=True)
    return root


def test_validate_run_id_rejects_path_escape() -> None:
    with pytest.raises(ValueError):
        validate_run_id("../evil")
    assert validate_run_id("stage0-flat-multitask-age-tissue-v1")


def test_parse_and_read_jsonl(tmp_path: Path) -> None:
    path = tmp_path / "metrics.jsonl"
    path.write_text(
        "\n".join(
            [
                json.dumps(
                    {
                        "epoch": 1,
                        "train_loss": 1.2,
                        "val_loss": 0.9,
                        "train_mae": 0.5,
                        "val_mae": 0.4,
                        "train_accuracy": 0.1,
                        "val_accuracy": 0.0,
                        "task": "multitask",
                    }
                ),
                "{not-json",
                json.dumps(
                    {
                        "epoch": 2,
                        "train_loss": 0.8,
                        "val_loss": 0.7,
                        "train_mae": 0.3,
                        "val_mae": 0.2,
                        "train_accuracy": 0.5,
                        "val_accuracy": 0.0,
                        "val_macro_f1": 0.25,
                        "task": "multitask",
                    }
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    rows = read_metrics_jsonl(path)
    assert len(rows) == 2
    assert rows[-1].epoch == 2
    assert rows[-1].val_macro_f1 == 0.25
    assert parse_metrics_row({"no_epoch": 1}) is None


def test_resolve_max_epochs_from_yaml(tmp_path: Path) -> None:
    cfg = tmp_path / "exp.yaml"
    cfg.write_text(yaml.safe_dump({"training": {"max_epochs": 40}}), encoding="utf-8")
    assert resolve_max_epochs(run_root=tmp_path, config_path=cfg, max_epochs_override=None) == 40
    assert resolve_max_epochs(run_root=tmp_path, config_path=cfg, max_epochs_override=12) == 12


def test_estimate_eta() -> None:
    history = [
        EpochMetrics(
            epoch=2,
            train_loss=1.0,
            val_loss=1.0,
            train_mae=None,
            val_mae=None,
            train_accuracy=None,
            val_accuracy=None,
            train_macro_f1=None,
            val_macro_f1=None,
            learning_rate=None,
            task="multitask",
        )
    ]
    eta, spe = estimate_eta(
        history=history,
        max_epochs=10,
        epoch_timestamps={1: 100.0, 2: 130.0},
    )
    assert spe == pytest.approx(30.0)
    assert eta == pytest.approx(8 * 30.0)


def test_collect_and_render(art_root: Path) -> None:
    run_id = "stage0-monitor-fixture-v1"
    run_root = art_root / "runs" / run_id
    ckpt_root = art_root / "checkpoints" / run_id
    run_root.mkdir(parents=True)
    ckpt_root.mkdir(parents=True)
    (run_root / "metrics.jsonl").write_text(
        json.dumps(
            {
                "epoch": 3,
                "train_loss": 0.5,
                "val_loss": 0.6,
                "train_mae": 1.1,
                "val_mae": 1.2,
                "train_accuracy": 0.8,
                "val_accuracy": 0.0,
                "learning_rate": 0.001,
                "task": "multitask",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (ckpt_root / "best.pt").write_bytes(b"x")
    (ckpt_root / "last.pt").write_bytes(b"y")
    cfg = run_root / "resolved_config.yaml"
    cfg.write_text(yaml.safe_dump({"training": {"max_epochs": 40}}), encoding="utf-8")

    snap = collect_snapshot(run_id=run_id, artifact_root=art_root)
    assert snap.latest is not None
    assert snap.latest.epoch == 3
    assert snap.max_epochs == 40
    assert snap.best_ckpt is not None
    rendered = render_snapshot(snap)
    assert rendered is not None


def test_monitor_cli_once(art_root: Path) -> None:
    run_id = "stage0-monitor-cli-v1"
    run_root = art_root / "runs" / run_id
    run_root.mkdir(parents=True)
    (run_root / "metrics.jsonl").write_text(
        json.dumps(
            {
                "epoch": 1,
                "train_loss": 1.0,
                "val_loss": 1.0,
                "train_mae": 0.0,
                "val_mae": 0.0,
                "train_accuracy": 0.0,
                "val_accuracy": 0.0,
                "task": "regression",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (run_root / "metrics.json").write_text("{}\n", encoding="utf-8")
    result = runner.invoke(
        app,
        ["monitor", "--run-id", run_id, "--once", "--max-epochs", "5"],
    )
    assert result.exit_code == 0, result.output
    assert (
        "Run finished" in result.output or "FINISHED" in result.output or "epoch" in result.output
    )
