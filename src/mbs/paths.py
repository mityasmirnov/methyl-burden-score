"""Filesystem policy for the MBS project."""

from __future__ import annotations

import os
from dataclasses import dataclass, fields
from pathlib import Path


class PathPolicyError(ValueError):
    """Raised when a configured project path is outside ``/data``."""


def _from_environment(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser().absolute()


def _is_under_data(path: Path) -> bool:
    data_root = Path("/data")
    try:
        path.relative_to(data_root)
    except ValueError:
        return False
    return True


@dataclass(frozen=True, slots=True)
class DataPaths:
    """Canonical project paths, all constrained to the server's ``/data`` volume."""

    project_root: Path
    data_root: Path
    scratch_root: Path
    cache_root: Path
    artifact_root: Path
    docker_root: Path

    @classmethod
    def from_environment(cls) -> DataPaths:
        """Construct paths from environment variables and validate them."""
        paths = cls(
            project_root=_from_environment(
                "MBS_PROJECT_ROOT",
                "/data/projects/methyl-burden-score",
            ),
            data_root=_from_environment(
                "MBS_DATA_ROOT",
                "/data/datasets/methyl-burden-score",
            ),
            scratch_root=_from_environment(
                "MBS_SCRATCH_ROOT",
                "/data/scratch/methyl-burden-score",
            ),
            cache_root=_from_environment(
                "MBS_CACHE_ROOT",
                "/data/cache/methyl-burden-score",
            ),
            artifact_root=_from_environment(
                "MBS_ARTIFACT_ROOT",
                "/data/artifacts/methyl-burden-score",
            ),
            docker_root=_from_environment("MBS_DOCKER_ROOT", "/data/docker"),
        )
        paths.validate()
        return paths

    def validate(self) -> None:
        """Require every configured path to be absolute and underneath ``/data``."""
        failures: list[str] = []
        for field in fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, Path):
                failures.append(f"{field.name} is not a pathlib.Path")
                continue
            if not value.is_absolute():
                failures.append(f"{field.name} is not absolute: {value}")
            elif not _is_under_data(value):
                failures.append(f"{field.name} is outside /data: {value}")

        if failures:
            raise PathPolicyError("; ".join(failures))

    def required_directories(self) -> tuple[Path, ...]:
        """Return directories required by the Stage 0 workspace."""
        return (
            self.project_root,
            self.data_root / "raw",
            self.data_root / "staging",
            self.data_root / "canonical" / "catalog" / "tables",
            self.data_root / "canonical" / "matrices",
            self.data_root / "canonical" / "annotations",
            self.data_root / "canonical" / "graphs",
            self.data_root / "canonical" / "static_features",
            self.scratch_root / "tmp",
            self.cache_root,
            self.artifact_root / "runs",
            self.artifact_root / "checkpoints",
            self.artifact_root / "scores",
            self.artifact_root / "reports",
            self.docker_root,
        )

    def ensure_directories(self) -> None:
        """Create required directories without touching locations outside ``/data``."""
        self.validate()
        for directory in self.required_directories():
            directory.mkdir(parents=True, exist_ok=True)

    def as_dict(self) -> dict[str, str]:
        """Return a JSON-serializable representation."""
        return {
            field.name: str(getattr(self, field.name))
            for field in fields(self)
        }
