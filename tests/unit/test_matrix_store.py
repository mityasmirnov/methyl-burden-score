"""Unit tests for EWAS Data Hub → canonical matrix conversion."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from mbs.cli import app
from mbs.matrix.convert import convert_ewas_db_study
from mbs.matrix.ewas_db import list_ewas_db_sample_files, read_ewas_db_sample
from mbs.matrix.roundtrip import verify_roundtrip
from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
    validate_matrix_manifest,
)

runner = CliRunner()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def isolated_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"matrix-{uuid4().hex}"
    workspace.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(workspace / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(workspace / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(workspace / "docker"))
    return workspace


def _write_hub_fixture(source_dir: Path) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)
    # Shared probe vocabulary; sample B misses cg00000003 (→ NaN in matrix).
    (source_dir / "GSM000001.txt").write_text(
        "cg00000001\t0.10\ncg00000002\t0.20\ncg00000003\t0.30\ncg99999999\t0.50\n",
        encoding="utf-8",
    )
    (source_dir / "GSM000002.txt").write_text(
        "cg00000001\t0.11\ncg00000002\t1.50\ncg00000003\t0.33\n",
        encoding="utf-8",
    )
    (source_dir / "GSM000003.txt").write_text(
        "cg00000001\t0.12\ncg00000002\t0.22\ncg00000003\t\n",
        encoding="utf-8",
    )


def _write_mini_annotations(annotations_dir: Path) -> None:
    annotations_dir.mkdir(parents=True, exist_ok=True)
    loci = pd.DataFrame(
        {
            "locus_id": np.array([101, 102, 103], dtype=np.uint64),
            "genome_build": ["GRCh38", "GRCh38", "GRCh38"],
            "chromosome": ["chr1", "chr1", "chr2"],
            "position": [100, 200, 300],
            "canonical_key": ["GRCh38:chr1:100", "GRCh38:chr1:200", "GRCh38:chr2:300"],
            "mapping_status": ["mapped", "mapped", "mapped"],
            "cpg_context": ["island", "shore", "open_sea"],
        }
    )
    edges = pd.DataFrame(
        {
            "probe_id": ["cg00000001", "cg00000002", "cg00000003"],
            "platform_id": ["HM450", "HM450", "HM450"],
            "locus_id": np.array([101, 102, 103], dtype=np.uint64),
            "mapping_source": ["test", "test", "test"],
            "mapping_confidence": [1.0, 1.0, 1.0],
            "is_primary": [True, True, True],
        }
    )
    probes = pd.DataFrame(
        {
            "probe_id": ["cg00000001", "cg00000002", "cg00000003"],
            "platform_id": ["HM450", "HM450", "HM450"],
            "probe_design": ["2", "2", "2"],
            "core_probe_id": ["cg00000001", "cg00000002", "cg00000003"],
            "M_mapping": [False, False, False],
            "M_nonuniq": [False, False, False],
            "M_general": [False, False, False],
            "mapQ": [60, 60, 60],
            "strand": ["+", "+", "-"],
            "mapping_status": ["mapped", "mapped", "mapped"],
        }
    )
    loci.to_parquet(annotations_dir / "loci.parquet", index=False)
    edges.to_parquet(annotations_dir / "probe_locus_edges.parquet", index=False)
    probes.to_parquet(annotations_dir / "probes.parquet", index=False)


def test_list_and_read_ewas_db_sample(tmp_path: Path) -> None:
    # Use /data scratch via repo path for policy-safe ops where needed; local tmp ok for reader
    source = tmp_path / "GSE_TEST"
    _write_hub_fixture(source)
    files = list_ewas_db_sample_files(source)
    assert [f.sample_id for f in files] == ["GSM000001", "GSM000002", "GSM000003"]
    table = read_ewas_db_sample(files[0].path)
    assert list(table.probe_ids) == ["cg00000001", "cg00000002", "cg00000003", "cg99999999"]
    assert table.betas[0] == pytest.approx(0.10)


def test_convert_roundtrip_fixture(isolated_workspace: Path) -> None:
    repo = _repo_root()
    source = isolated_workspace / "data" / "raw" / "ewas_datahub" / "EWAS_db" / "GSE_TEST"
    annotations = isolated_workspace / "data" / "canonical" / "annotations"
    output = isolated_workspace / "data" / "canonical" / "matrices" / "matrix-test-v1"
    report = isolated_workspace / "reports" / "inspection" / "GSE_TEST_ewas_db"
    _write_hub_fixture(source)
    _write_mini_annotations(annotations)

    result = convert_ewas_db_study(
        project_root=repo,
        source_dir=source,
        annotations_dir=annotations,
        output_dir=output,
        study_id="GSE_TEST",
        matrix_id="matrix-test-v1",
        platform_id="HM450",
        processing_level="gmqn",
        report_dir=report,
        verify=True,
        verify_max_probes=10,
    )
    assert result.roundtrip is not None
    assert result.roundtrip.ok
    assert result.stats["n_samples"] == 3
    assert result.stats["n_study_loci"] == 3
    assert result.stats["n_unmapped_probes"] == 1
    assert result.stats["n_out_of_range"] == 1  # 1.50 in GSM000002
    assert (report / "summary.md").is_file()

    paths = matrix_store_paths(output)
    manifest = json.loads(paths.manifest_path.read_text(encoding="utf-8"))
    validate_matrix_manifest(manifest)
    assert manifest["shape"] == [3, 3]
    assert manifest["dtype"] == "float32"
    assert len(manifest["source_files"]) == 3

    sample_index = read_sample_index(paths.sample_index_path)
    locus_index = read_locus_index(paths.locus_index_path)
    betas = open_betas_zarr(paths.betas_path)
    assert list(sample_index["sample_id"]) == ["GSM000001", "GSM000002", "GSM000003"]
    assert list(locus_index["probe_id"]) == ["cg00000001", "cg00000002", "cg00000003"]
    assert np.asarray(betas[0, 0]) == pytest.approx(np.float32(0.10))
    assert np.asarray(betas[1, 1]) == pytest.approx(np.float32(1.50))  # not clipped
    assert np.isnan(np.asarray(betas[2, 2]))

    rt = verify_roundtrip(source, output, max_probes=None)
    assert rt.ok
    assert rt.n_mismatch == 0


def test_matrix_convert_cli_help() -> None:
    result = runner.invoke(app, ["matrix", "convert", "--help"])
    assert result.exit_code == 0, result.output
    assert "EWAS Data Hub" in result.output or "canonical" in result.output.lower()


def test_matrix_convert_cli_fixture(isolated_workspace: Path) -> None:
    source = isolated_workspace / "data" / "raw" / "ewas_datahub" / "EWAS_db" / "GSE_TEST"
    annotations = isolated_workspace / "data" / "canonical" / "annotations"
    output = isolated_workspace / "data" / "canonical" / "matrices" / "matrix-cli-v1"
    report = isolated_workspace / "reports" / "inspection" / "GSE_TEST_ewas_db"
    _write_hub_fixture(source)
    _write_mini_annotations(annotations)

    result = runner.invoke(
        app,
        [
            "matrix",
            "convert",
            "--study-id",
            "GSE_TEST",
            "--source-dir",
            str(source),
            "--annotations-dir",
            str(annotations),
            "--output-dir",
            str(output),
            "--matrix-id",
            "matrix-cli-v1",
            "--platform-id",
            "HM450",
            "--report-dir",
            str(report),
            "--verify",
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_samples"] == 3
    assert payload["roundtrip_ok"] is True
    assert (output / "matrix_manifest.json").is_file()
