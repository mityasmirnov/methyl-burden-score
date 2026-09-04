from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import pytest
from typer.main import get_command
from typer.testing import CliRunner

from mbs.cli import app
from mbs.inspect_source import inventory_source, write_inspection_report

runner = CliRunner()

EXPECTED_TABLES = 31
EXPECTED_VIEWS = 23
POLICY_EXIT_CODE = 2


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _point_env_at_scratch(monkeypatch: pytest.MonkeyPatch, scratch: Path) -> None:
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
    """Create a disposable workspace under the real project scratch (/data)."""
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"cli-{uuid4().hex}"
    workspace.mkdir()
    _point_env_at_scratch(monkeypatch, workspace)
    return workspace


def test_convert_pack_help_lists_bmi_ancestry() -> None:
    convert = get_command(app).commands["matrix"].commands["convert-pack"]
    family_help = next(p.help or "" for p in convert.params if p.name == "phenotype_family")
    assert "bmi" in family_help
    assert "ancestry" in family_help


def test_doctor_ok(isolated_workspace: Path) -> None:
    result = runner.invoke(app, ["doctor", "--create-directories"])
    assert result.exit_code == 0, result.output
    assert "MBS environment" in result.output
    assert (isolated_workspace / "data" / "raw" / "cpgcorpus").is_dir()


def test_doctor_rejects_home_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MBS_DATA_ROOT", "/home/smirnov/mbs-data")
    result = runner.invoke(app, ["doctor"])
    assert result.exit_code == POLICY_EXIT_CODE
    assert "Path policy failure" in result.output


def test_catalog_init_cli(isolated_workspace: Path) -> None:
    result = runner.invoke(app, ["catalog", "init"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["tables"] == EXPECTED_TABLES
    assert payload["views"] == EXPECTED_VIEWS
    assert Path(payload["database"]).exists()
    assert Path(payload["database"]).is_relative_to(isolated_workspace / "data")


def test_inspect_source_empty_dir(isolated_workspace: Path) -> None:
    raw = isolated_workspace / "data" / "raw" / "cpgcorpus"
    raw.mkdir(parents=True, exist_ok=True)
    report = isolated_workspace / "reports" / "inspection" / "cpgcorpus"

    result = runner.invoke(
        app,
        [
            "inspect",
            "source",
            "--source-id",
            "cpgcorpus",
            "--report-dir",
            str(report),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["file_count"] == 0
    assert (report / "summary.json").exists()
    assert (report / "summary.md").exists()


def test_inventory_source_counts_files(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    nested = raw / "GSE1" / "GPL1"
    nested.mkdir(parents=True)
    (nested / "metadata.arrow").write_bytes(b"abc")
    (nested / "note.txt").write_text("hi", encoding="utf-8")

    inventory = inventory_source(raw, source_id="demo", max_entries=10)
    assert inventory["file_count"] == 2
    assert inventory["total_bytes"] == 5
    assert inventory["suffix_counts"][".arrow"] == 1
    assert inventory["suffix_counts"][".txt"] == 1

    report_dir = write_inspection_report(inventory, tmp_path / "report")
    summary = json.loads((report_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["source_id"] == "demo"
