"""Unit tests for EWAS Data Hub baseline pack → canonical matrix conversion."""

from __future__ import annotations

import zipfile
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest

from mbs.matrix.hub_pack import convert_hub_pack_subset, stream_pack_betas
from mbs.matrix.store import open_betas_zarr, read_sample_index


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def isolated_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"hubpack-{uuid4().hex}"
    workspace.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(workspace / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(workspace / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(workspace / "docker"))
    return workspace


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
            "is_primary": [True, True, True],
        }
    )
    loci.to_parquet(annotations_dir / "loci.parquet", index=False)
    edges.to_parquet(annotations_dir / "probe_locus_edges.parquet", index=False)


def _write_tiny_age_pack(data_root: Path) -> tuple[Path, Path]:
    """Write a miniature Hub age pack zip + sample-info parquet."""
    download = data_root / "raw" / "ewas_datahub" / "download"
    download.mkdir(parents=True, exist_ok=True)
    phenotypes = data_root / "canonical" / "phenotypes"
    phenotypes.mkdir(parents=True, exist_ok=True)

    # Header + metadata + three probes; sample order GSM_A, GSM_B, GSM_C
    lines = [
        "sample_id\tGSM_A\tGSM_B\tGSM_C\n",
        "age\t40\t50\t60\n",
        "tissue\tblood\tblood\tblood\n",
        "cg00000001\t0.10\t0.20\t0.30\n",
        "cg00000002\t0.11\tNA\t0.31\n",
        "cg00000003\t0.12\t0.22\t0.32\n",
        "not_a_probe\t1\t2\t3\n",
    ]
    zip_path = download / "age_methylation_v1.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("age_methylation_v1.txt", "".join(lines))

    info = pd.DataFrame(
        {
            "sample_id": ["GSM_A", "GSM_B", "GSM_C"],
            "study_id": ["GSE_TRAIN", "GSE_TRAIN", "GSE_VAL"],
            "platform": ["450K", "450K", "450K"],
            "sample_type": ["normal", "normal", "normal"],
            "phenotype_value": ["40", "50", "60"],
            "phenotype_value_numeric": [40.0, 50.0, 60.0],
            "phenotype_family": ["age", "age", "age"],
        }
    )
    info_path = phenotypes / "age_sample_info.parquet"
    info.to_parquet(info_path, index=False)
    return zip_path, info_path


def test_stream_pack_betas_column_order(isolated_workspace: Path) -> None:
    data_root = isolated_workspace / "data"
    zip_path, _ = _write_tiny_age_pack(data_root)
    # Request samples in reverse pack order to exercise reorder path.
    betas, probes, meta = stream_pack_betas(
        zip_path=zip_path,
        family="age",
        sample_ids=["GSM_C", "GSM_A"],
    )
    assert list(probes) == ["cg00000001", "cg00000002", "cg00000003"]
    assert meta["n_selected_samples"] == 2
    np.testing.assert_allclose(betas[0], [0.30, 0.31, 0.32], rtol=0, atol=1e-6)
    np.testing.assert_allclose(betas[1], [0.10, 0.11, 0.12], rtol=0, atol=1e-6)


def test_convert_hub_pack_subset(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    _write_tiny_age_pack(data_root)
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations)
    out = data_root / "canonical" / "matrices" / "matrix-hub-age-tiny-v1"
    result = convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="age",
        study_ids=["GSE_TRAIN", "GSE_VAL"],
        matrix_id="matrix-hub-age-tiny-v1",
        output_dir=out,
        max_per_study=2,
    )
    assert result.stats["n_samples"] == 3
    assert (out / "sample_phenotypes.parquet").is_file()
    assert (out / "study_subset.json").is_file()
    sample_index = read_sample_index(out / "sample_index.parquet")
    assert sample_index["sample_id"].tolist() == ["GSM_A", "GSM_B", "GSM_C"]
    betas = np.asarray(open_betas_zarr(out / "betas.zarr")[:])
    assert betas.shape == (3, 3)
    assert np.isnan(betas[1, 1])
