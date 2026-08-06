"""Phenotype / source dataset registry (Milestone 5b)."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import yaml

SourceSystem = Literal["ewas_datahub", "ewas_atlas", "cpgcorpus"]
PhenotypeFamily = Literal[
    "age",
    "tissue",
    "disease",
    "cancer",
    "blood",
    "brain",
    "sex",
    "ancestry",
    "bmi",
    "cell_type",
]
LabelType = Literal[
    "regression",
    "binary",
    "multiclass",
    "pack_profile",
    "sample_info",
    "none",
]
SplitRole = Literal[
    "train",
    "validation",
    "external_test",
    "pilot",
    "registered",
    "secondary",
    "benchmark",
]

_SOURCE_SYSTEMS = frozenset({"ewas_datahub", "ewas_atlas", "cpgcorpus"})
_FAMILIES = frozenset(
    {
        "age",
        "tissue",
        "disease",
        "cancer",
        "blood",
        "brain",
        "sex",
        "ancestry",
        "bmi",
        "cell_type",
    }
)
_LABEL_TYPES = frozenset(
    {"regression", "binary", "multiclass", "pack_profile", "sample_info", "none"}
)
_SPLIT_ROLES = frozenset(
    {
        "train",
        "validation",
        "external_test",
        "pilot",
        "registered",
        "secondary",
        "benchmark",
    }
)


@dataclass(frozen=True, slots=True)
class RegistryEntry:
    entry_id: str
    source_system: SourceSystem
    phenotype_family: PhenotypeFamily
    label_type: LabelType
    split_role: SplitRole
    download_path: str
    study_id: str | None = None
    study_ids: tuple[str, ...] | None = None
    platform: str | None = None
    sample_count: int | None = None
    matrix_path: str | None = None
    checksum: str | None = None
    notes: str | None = None

    def resolve_download_path(self, data_root: Path) -> Path:
        rel = Path(self.download_path)
        return rel if rel.is_absolute() else (data_root / rel).resolve()

    def resolve_matrix_path(self, data_root: Path) -> Path | None:
        if self.matrix_path is None:
            return None
        rel = Path(self.matrix_path)
        return rel if rel.is_absolute() else (data_root / rel).resolve()


@dataclass(frozen=True, slots=True)
class PhenotypeRegistry:
    registry_version: str
    entries: tuple[RegistryEntry, ...]
    description: str | None = None

    def by_id(self) -> dict[str, RegistryEntry]:
        return {entry.entry_id: entry for entry in self.entries}

    def filter(
        self,
        *,
        phenotype_family: str | None = None,
        split_role: str | None = None,
        source_system: str | None = None,
    ) -> list[RegistryEntry]:
        out: list[RegistryEntry] = []
        for entry in self.entries:
            if phenotype_family is not None and entry.phenotype_family != phenotype_family:
                continue
            if split_role is not None and entry.split_role != split_role:
                continue
            if source_system is not None and entry.source_system != source_system:
                continue
            out.append(entry)
        return out


def _require_str(payload: dict[str, Any], key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"registry entry missing non-empty string {key!r}")
    return value


def validate_phenotype_registry(payload: dict[str, Any]) -> None:
    """Validate against ``schemas/phenotype_registry.schema.json`` rules."""
    if "registry_version" not in payload:
        raise ValueError("phenotype registry missing registry_version")
    if not isinstance(payload["registry_version"], str) or not payload["registry_version"]:
        raise ValueError("registry_version must be a non-empty string")
    entries = payload.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ValueError("phenotype registry entries must be a non-empty list")
    seen: set[str] = set()
    for raw in entries:
        if not isinstance(raw, dict):
            raise TypeError("each registry entry must be a mapping")
        entry_id = _require_str(raw, "entry_id")
        if entry_id in seen:
            raise ValueError(f"duplicate registry entry_id: {entry_id}")
        seen.add(entry_id)
        source = _require_str(raw, "source_system")
        if source not in _SOURCE_SYSTEMS:
            raise ValueError(f"invalid source_system: {source}")
        family = _require_str(raw, "phenotype_family")
        if family not in _FAMILIES:
            raise ValueError(f"invalid phenotype_family: {family}")
        label_type = _require_str(raw, "label_type")
        if label_type not in _LABEL_TYPES:
            raise ValueError(f"invalid label_type: {label_type}")
        split_role = _require_str(raw, "split_role")
        if split_role not in _SPLIT_ROLES:
            raise ValueError(f"invalid split_role: {split_role}")
        _require_str(raw, "download_path")
        checksum = raw.get("checksum")
        if checksum is not None and not (
            isinstance(checksum, str)
            and len(checksum) == 64
            and set(checksum) <= set("0123456789abcdef")
        ):
            raise ValueError(f"checksum must be 64 lowercase hex or null: {entry_id}")
        sample_count = raw.get("sample_count")
        if sample_count is not None and (not isinstance(sample_count, int) or sample_count < 0):
            raise ValueError(f"sample_count must be non-negative int or null: {entry_id}")


def _entry_from_dict(raw: dict[str, Any]) -> RegistryEntry:
    study_ids_raw = raw.get("study_ids")
    study_ids: tuple[str, ...] | None
    if study_ids_raw is None:
        study_ids = None
    elif isinstance(study_ids_raw, list) and all(isinstance(x, str) for x in study_ids_raw):
        study_ids = tuple(study_ids_raw)
    else:
        raise TypeError(f"study_ids must be a list of strings or null: {raw.get('entry_id')}")
    return RegistryEntry(
        entry_id=_require_str(raw, "entry_id"),
        source_system=_require_str(raw, "source_system"),  # type: ignore[arg-type]
        phenotype_family=_require_str(raw, "phenotype_family"),  # type: ignore[arg-type]
        label_type=_require_str(raw, "label_type"),  # type: ignore[arg-type]
        split_role=_require_str(raw, "split_role"),  # type: ignore[arg-type]
        download_path=_require_str(raw, "download_path"),
        study_id=raw.get("study_id"),
        study_ids=study_ids,
        platform=raw.get("platform"),
        sample_count=raw.get("sample_count"),
        matrix_path=raw.get("matrix_path"),
        checksum=raw.get("checksum"),
        notes=raw.get("notes"),
    )


def load_phenotype_registry(path: Path) -> PhenotypeRegistry:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"phenotype registry not found: {path}")
    with path.open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"phenotype registry must be a mapping: {path}")
    validate_phenotype_registry(payload)
    entries = tuple(_entry_from_dict(raw) for raw in payload["entries"])
    description = payload.get("description")
    return PhenotypeRegistry(
        registry_version=str(payload["registry_version"]),
        entries=entries,
        description=str(description) if description is not None else None,
    )


def default_registry_path(project_root: Path) -> Path:
    return project_root / "configs" / "data" / "phenotype_registry.yaml"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_download_checksums(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> Path:
    """Write / merge download checksum sidecar Parquet under canonical/registries."""
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    required = {"entry_id", "download_path", "sha256", "bytes"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"checksum rows missing columns: {sorted(missing)}")
    if output_path.is_file():
        existing = pd.read_parquet(output_path)
        merged = pd.concat([existing, frame], ignore_index=True)
        merged = merged.drop_duplicates(subset=["entry_id"], keep="last")
    else:
        merged = frame
    merged.to_parquet(output_path, index=False)
    return output_path


def export_registry_parquet(registry: PhenotypeRegistry, output_path: Path) -> Path:
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "entry_id": e.entry_id,
            "source_system": e.source_system,
            "phenotype_family": e.phenotype_family,
            "study_id": e.study_id,
            "study_ids": None if e.study_ids is None else ",".join(e.study_ids),
            "platform": e.platform,
            "sample_count": e.sample_count,
            "label_type": e.label_type,
            "split_role": e.split_role,
            "download_path": e.download_path,
            "matrix_path": e.matrix_path,
            "checksum": e.checksum,
            "notes": e.notes,
        }
        for e in registry.entries
    ]
    pd.DataFrame(rows).to_parquet(output_path, index=False)
    return output_path
