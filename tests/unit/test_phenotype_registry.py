"""Unit tests for the phenotype / source dataset registry."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from mbs.registry import (
    load_phenotype_registry,
    validate_phenotype_registry,
    write_download_checksums,
)


def test_seed_registry_loads() -> None:
    root = Path(__file__).resolve().parents[2]
    path = root / "configs" / "data" / "phenotype_registry.yaml"
    registry = load_phenotype_registry(path)
    assert registry.registry_version == "phenotype-registry-v1"
    assert registry.by_id()["gse35069-ewasdb-celltype-pilot"].study_id == "GSE35069"
    age = registry.filter(phenotype_family="age")
    assert len(age) >= 2
    assert any(e.label_type == "pack_profile" for e in age)


def test_validate_rejects_bad_source() -> None:
    payload = {
        "registry_version": "t",
        "entries": [
            {
                "entry_id": "x",
                "source_system": "not_a_source",
                "phenotype_family": "age",
                "label_type": "regression",
                "split_role": "train",
                "download_path": "raw/x",
            }
        ],
    }
    with pytest.raises(ValueError, match="source_system"):
        validate_phenotype_registry(payload)


def test_write_download_checksums_merges(tmp_path: Path) -> None:
    out = tmp_path / "download_checksums.parquet"
    write_download_checksums(
        [{"entry_id": "a", "download_path": "p/a", "sha256": "a" * 64, "bytes": 1}],
        out,
    )
    write_download_checksums(
        [{"entry_id": "a", "download_path": "p/a", "sha256": "b" * 64, "bytes": 2}],
        out,
    )
    frame = pd.read_parquet(out)
    assert len(frame) == 1
    assert frame.iloc[0]["sha256"] == "b" * 64


def test_fixture_yaml_roundtrip(tmp_path: Path) -> None:
    payload = {
        "registry_version": "fixture-v1",
        "entries": [
            {
                "entry_id": "fixture-age",
                "source_system": "ewas_datahub",
                "phenotype_family": "age",
                "study_id": "GSE1",
                "platform": "HM450",
                "sample_count": 3,
                "label_type": "regression",
                "split_role": "train",
                "download_path": "raw/fake.zip",
                "checksum": None,
                "notes": "unit fixture",
            }
        ],
    }
    path = tmp_path / "reg.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    registry = load_phenotype_registry(path)
    assert registry.entries[0].sample_count == 3
