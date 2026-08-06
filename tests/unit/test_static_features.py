"""Unit tests for static locus feature export helpers (no CpGPT runtime)."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import torch

from mbs.static_features.coordinates import mbs_locus_to_cpgpt_location
from mbs.static_features.cpgpt_adapter import SequenceAdapterMLP
from mbs.static_features.manifest import (
    ARTIFACT_VERSION,
    SOURCE_REPOSITORY,
    validate_static_feature_manifest,
)
from mbs.static_features.store import (
    open_embeddings_zarr,
    read_loci_index,
    static_feature_store_paths,
    write_artifact,
    write_embeddings_zarr,
    write_loci_index,
)
from mbs.static_features.validate_export import (
    embedding_summary_stats,
    validate_embeddings_array,
    validate_loci_frame,
)


def test_mbs_locus_to_cpgpt_location_primary() -> None:
    assert mbs_locus_to_cpgpt_location("chr1", 10848) == "1:10847"
    assert mbs_locus_to_cpgpt_location("chrM", 100) == "MT:99"
    assert mbs_locus_to_cpgpt_location("1", 10) == "1:9"


def test_mbs_locus_to_cpgpt_location_rejects_alt() -> None:
    assert mbs_locus_to_cpgpt_location("chr14_GL000009v2_random", 100) is None
    assert mbs_locus_to_cpgpt_location("chr1", 0) is None
    assert mbs_locus_to_cpgpt_location("", 10) is None


def _valid_manifest(tmp_path: Path) -> dict[str, Any]:
    embedding = tmp_path / "embeddings.zarr"
    loci = tmp_path / "loci.parquet"
    # Paths must be under /data for schema validation; use absolute /data-looking
    # strings that do not need to exist for validate_static_feature_manifest.
    return {
        "artifact_version": ARTIFACT_VERSION,
        "feature_set_id": "cpgpt2m_adapter_128_v1",
        "source_model": "nucleotide-transformer-v2-500m-multi-species",
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": "a" * 40,
        "checkpoint_sha256": "b" * 64,
        "configuration_sha256": "c" * 64,
        "vocabulary_sha256": None,
        "context_length": 2001,
        "genome_build": "GRCh38",
        "input_dimension": 1024,
        "output_dimension": 128,
        "storage_dtype": "float16",
        "normalization": None,
        "n_loci": 2,
        "locus_table_sha256": "d" * 64,
        "embedding_path": (
            f"/data/projects/methyl-burden-score/data/canonical/static_features/x/{embedding.name}"
        ),
        "locus_index_path": (
            f"/data/projects/methyl-burden-score/data/canonical/static_features/x/{loci.name}"
        ),
        "export_command": "uv run --extra cpgpt mbs features export-cpgpt",
        "created_at": "2026-08-06T00:00:00Z",
        "notes": "unit test",
    }


def test_validate_static_feature_manifest_accepts_valid(tmp_path: Path) -> None:
    validate_static_feature_manifest(_valid_manifest(tmp_path))


def test_validate_static_feature_manifest_rejects_bad_genome(tmp_path: Path) -> None:
    manifest = _valid_manifest(tmp_path)
    manifest["genome_build"] = "hg19"
    with pytest.raises(ValueError, match="GRCh38"):
        validate_static_feature_manifest(manifest)


def test_validate_static_feature_manifest_rejects_non_data_path(tmp_path: Path) -> None:
    manifest = _valid_manifest(tmp_path)
    manifest["embedding_path"] = "/tmp/embeddings.zarr"
    with pytest.raises(ValueError, match="/data"):
        validate_static_feature_manifest(manifest)


def test_store_round_trip(tmp_path: Path) -> None:
    store_root = tmp_path / "cpgpt2m_adapter_128_v1"
    paths = static_feature_store_paths(store_root)
    embeddings = np.arange(8, dtype=np.float32).reshape(2, 4).astype(np.float16)
    write_embeddings_zarr(paths.embeddings_path, embeddings)
    loci = pd.DataFrame(
        {
            "embedding_row": pd.array([0, 1, None], dtype="Int64"),
            "locus_id": np.array([1, 2, 3], dtype=np.uint64),
            "canonical_key": ["GRCh38:chr1:10", "GRCh38:chr1:20", "GRCh38:chr1:30"],
            "source_location_key": ["1:9", "1:19", None],
            "source_embedding_row": pd.array([5, 6, None], dtype="Int64"),
            "mapping_status": ["mapped", "mapped", "missing"],
        }
    )
    write_loci_index(paths.loci_path, loci)

    loaded = np.asarray(open_embeddings_zarr(paths.embeddings_path)[:])
    assert loaded.shape == (2, 4)
    assert loaded.dtype == np.float16
    reloaded = read_loci_index(paths.loci_path)
    assert list(reloaded["mapping_status"]) == ["mapped", "mapped", "missing"]
    validate_embeddings_array(loaded, output_dimension=4)
    coverage = validate_loci_frame(reloaded, n_mapped=2)
    assert coverage == {"n_loci": 3, "n_mapped": 2, "n_missing": 1}
    stats = embedding_summary_stats(loaded)
    assert stats["n_rows"] == 2
    assert stats["n_dims"] == 4


def test_sequence_adapter_mlp_shapes() -> None:
    adapter = SequenceAdapterMLP(
        1024,
        128,
        128,
        dropout=0.0,
        n_blocks=3,
        pre_norm=False,
        post_norm=False,
    )
    adapter.eval()
    with torch.inference_mode():
        out = adapter(torch.randn(2, 1, 1024))
    assert out.shape == (2, 1, 128)


def test_write_artifact_under_data_root() -> None:
    data_root = Path("/data/projects/methyl-burden-score/data")
    if not data_root.is_dir():
        pytest.skip("/data project root unavailable")
    out = data_root / "scratch" / "tmp" / "test_static_features_artifact"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    paths = static_feature_store_paths(out)
    embeddings = np.zeros((1, 128), dtype=np.float16)
    write_embeddings_zarr(paths.embeddings_path, embeddings)
    loci = pd.DataFrame(
        {
            "embedding_row": pd.array([0], dtype="Int64"),
            "locus_id": np.array([1], dtype=np.uint64),
            "canonical_key": ["GRCh38:chr1:10"],
            "source_location_key": ["1:9"],
            "source_embedding_row": pd.array([0], dtype="Int64"),
            "mapping_status": ["mapped"],
        }
    )
    write_loci_index(paths.loci_path, loci)
    manifest = _valid_manifest(out)
    manifest["embedding_path"] = str(paths.embeddings_path)
    manifest["locus_index_path"] = str(paths.loci_path)
    manifest["n_loci"] = 1
    write_artifact(paths.artifact_path, manifest)
    assert paths.artifact_path.is_file()
