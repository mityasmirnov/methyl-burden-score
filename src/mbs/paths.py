"""Filesystem policy for the MBS project."""

from __future__ import annotations

import os
from dataclasses import dataclass, fields
from pathlib import Path


class PathPolicyError(ValueError):
    """Raised when a configured project path is outside ``/data``."""


def get_env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser().resolve()


def _default_repo_root() -> Path:
    # src/mbs/paths.py -> repo root
    return Path(__file__).resolve().parents[2]


def _is_under_data(path: Path) -> bool:
    try:
        path.relative_to(Path("/data"))
    except ValueError:
        return False
    return True


@dataclass(frozen=True, slots=True)
class DataPaths:
    """Canonical project paths, constrained to the server's ``/data`` volume."""

    project_root: Path
    data_root: Path
    scratch_root: Path
    cache_root: Path
    artifact_root: Path
    docker_root: Path

    @classmethod
    def from_environment(cls) -> DataPaths:
        """Construct paths from environment variables and validate them."""
        root = get_env_path(
            "MBS_ROOT",
            os.environ.get("MBS_PROJECT_ROOT", str(_default_repo_root())),
        )
        paths = cls(
            project_root=root,
            data_root=get_env_path("MBS_DATA_ROOT", str(root / "data")),
            scratch_root=get_env_path("MBS_SCRATCH_ROOT", str(root / "scratch")),
            cache_root=get_env_path("MBS_CACHE_ROOT", str(root / "cache")),
            artifact_root=get_env_path("MBS_ARTIFACT_ROOT", str(root / "artifacts")),
            docker_root=get_env_path("MBS_DOCKER_ROOT", str(root / "docker")),
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
            self.data_root / "raw" / "cpgcorpus",
            self.data_root / "raw" / "cpgcorpus" / "_partial_fullsync",
            self.data_root / "raw" / "ewas_datahub",
            self.data_root / "raw" / "ewas_datahub" / "download",
            self.data_root / "raw" / "ewas_datahub" / "EWAS_db",
            self.data_root / "raw" / "ewas_atlas",
            self.data_root / "raw" / "manifests",
            self.data_root / "raw" / "manifests" / "epicv2",
            self.data_root / "raw" / "gencode",
            self.data_root / "raw" / "annotations",
            self.data_root / "staging",
            self.data_root / "staging" / "infinium_export",
            self.scratch_root / "downloads",
            self.artifact_root / "logs" / "downloads",
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
        return {field.name: str(getattr(self, field.name)) for field in fields(self)}


# Module-level convenience aliases (resolved at import from the current env).
MBS_ROOT = get_env_path(
    "MBS_ROOT",
    os.environ.get("MBS_PROJECT_ROOT", str(_default_repo_root())),
)
MBS_DATA_ROOT = get_env_path("MBS_DATA_ROOT", str(MBS_ROOT / "data"))
MBS_SCRATCH_ROOT = get_env_path("MBS_SCRATCH_ROOT", str(MBS_ROOT / "scratch"))
MBS_CACHE_ROOT = get_env_path("MBS_CACHE_ROOT", str(MBS_ROOT / "cache"))
MBS_ARTIFACT_ROOT = get_env_path("MBS_ARTIFACT_ROOT", str(MBS_ROOT / "artifacts"))
