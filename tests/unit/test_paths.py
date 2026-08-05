from __future__ import annotations

from pathlib import Path

import pytest

from mbs.paths import DataPaths, PathPolicyError


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def test_default_paths_are_project_local(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in (
        "MBS_ROOT",
        "MBS_PROJECT_ROOT",
        "MBS_DATA_ROOT",
        "MBS_SCRATCH_ROOT",
        "MBS_CACHE_ROOT",
        "MBS_ARTIFACT_ROOT",
        "MBS_DOCKER_ROOT",
    ):
        monkeypatch.delenv(name, raising=False)

    paths = DataPaths.from_environment()
    root = _repo_root()

    assert paths.project_root == root
    assert paths.data_root == root / "data"
    assert paths.scratch_root == root / "scratch"
    assert paths.cache_root == root / "cache"
    assert paths.artifact_root == root / "artifacts"
    assert paths.docker_root == root / "docker"
    for value in paths.as_dict().values():
        assert Path(value).is_relative_to(Path("/data"))
        assert Path(value).is_relative_to(root) or Path(value) == root


def test_home_path_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MBS_DATA_ROOT", "/home/smirnov/data")

    with pytest.raises(PathPolicyError, match="outside /data"):
        DataPaths.from_environment()


def test_relative_path_is_rejected() -> None:
    root = _repo_root()
    paths = DataPaths(
        project_root=Path("relative/project"),
        data_root=root / "data",
        scratch_root=root / "scratch",
        cache_root=root / "cache",
        artifact_root=root / "artifacts",
        docker_root=root / "docker",
    )

    with pytest.raises(PathPolicyError, match="not absolute"):
        paths.validate()


def test_tmp_home_cache_aliases_are_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MBS_CACHE_ROOT", "/tmp/mbs-cache")

    with pytest.raises(PathPolicyError, match="outside /data"):
        DataPaths.from_environment()
