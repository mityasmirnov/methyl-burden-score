from __future__ import annotations

from pathlib import Path

import pytest

from mbs.paths import DataPaths, PathPolicyError


def test_default_paths_are_under_data(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in (
        "MBS_PROJECT_ROOT",
        "MBS_DATA_ROOT",
        "MBS_SCRATCH_ROOT",
        "MBS_CACHE_ROOT",
        "MBS_ARTIFACT_ROOT",
        "MBS_DOCKER_ROOT",
    ):
        monkeypatch.delenv(name, raising=False)

    paths = DataPaths.from_environment()

    for value in paths.as_dict().values():
        assert Path(value).is_relative_to(Path("/data"))


def test_home_path_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MBS_DATA_ROOT", "/home/smirnov/data")

    with pytest.raises(PathPolicyError, match="outside /data"):
        DataPaths.from_environment()


def test_relative_path_is_rejected() -> None:
    paths = DataPaths(
        project_root=Path("relative/project"),
        data_root=Path("/data/datasets/methyl-burden-score"),
        scratch_root=Path("/data/scratch/methyl-burden-score"),
        cache_root=Path("/data/cache/methyl-burden-score"),
        artifact_root=Path("/data/artifacts/methyl-burden-score"),
        docker_root=Path("/data/docker"),
    )

    with pytest.raises(PathPolicyError, match="not absolute"):
        paths.validate()
