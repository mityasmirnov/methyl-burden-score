"""Unit tests for EWAS Data Hub baseline pack → canonical matrix conversion."""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest

from mbs.annotation.manifest import sha256_file
from mbs.matrix.hub_pack import (
    SUPPORTED_PACK_FAMILIES,
    convert_hub_pack_subset,
    pack_txt_name,
    pack_zip_path,
    select_samples_for_studies,
    stream_pack_betas,
)
from mbs.matrix.hub_pack_index import (
    build_hub_pack_matrix_index,
    check_overlapping_gsm_betas,
)
from mbs.matrix.locus_map import (
    COLLAPSE_MEAN,
    COLLAPSE_MEDIAN,
    build_probe_locus_map,
)
from mbs.matrix.store import open_betas_zarr, read_locus_index, read_sample_index


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


def _write_mini_annotations(annotations_dir: Path, *, collapse: bool = False) -> None:
    annotations_dir.mkdir(parents=True, exist_ok=True)
    if collapse:
        loci = pd.DataFrame(
            {
                "locus_id": np.array([101, 102], dtype=np.uint64),
                "genome_build": ["GRCh38", "GRCh38"],
                "chromosome": ["chr1", "chr1"],
                "position": [100, 200],
                "canonical_key": ["GRCh38:chr1:100", "GRCh38:chr1:200"],
                "mapping_status": ["mapped", "mapped"],
                "cpg_context": ["island", "shore"],
            }
        )
        edges = pd.DataFrame(
            {
                "probe_id": ["cg00000001", "cg00000002", "cg00000003"],
                "platform_id": ["HM450", "HM450", "HM450"],
                "locus_id": np.array([101, 101, 102], dtype=np.uint64),
                "is_primary": [True, True, True],
            }
        )
    else:
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


def _write_tiny_pack(
    data_root: Path,
    *,
    family: str = "age",
    lines: list[str] | None = None,
    info: pd.DataFrame | None = None,
) -> tuple[Path, Path]:
    download = data_root / "raw" / "ewas_datahub" / "download"
    download.mkdir(parents=True, exist_ok=True)
    phenotypes = data_root / "canonical" / "phenotypes"
    phenotypes.mkdir(parents=True, exist_ok=True)

    if lines is None:
        lines = [
            "sample_id\tGSM_A\tGSM_B\tGSM_C\n",
            "age\t40\t50\t60\n",
            "tissue\tblood\tblood\tblood\n",
            "cg00000001\t0.10\t0.20\t0.30\n",
            "cg00000002\t0.11\tNA\t0.31\n",
            "cg00000003\t0.12\t0.22\t0.32\n",
            "not_a_probe\t1\t2\t3\n",
        ]
    zip_path = pack_zip_path(data_root, family)
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(pack_txt_name(family), "".join(lines))

    if info is None:
        info = pd.DataFrame(
            {
                "sample_id": ["GSM_A", "GSM_B", "GSM_C"],
                "study_id": ["GSE_TRAIN", "GSE_TRAIN", "GSE_VAL"],
                "platform": ["450K", "450K", "450K"],
                "sample_type": ["normal", "normal", "normal"],
                "phenotype_value": ["40", "50", "60"],
                "phenotype_value_numeric": [40.0, 50.0, 60.0],
                "phenotype_family": [family, family, family],
            }
        )
    info_path = phenotypes / f"{family}_sample_info.parquet"
    info.to_parquet(info_path, index=False)
    return zip_path, info_path


def test_bmi_ancestry_in_pack_maps(isolated_workspace: Path) -> None:
    data_root = isolated_workspace / "data"
    assert "bmi" in SUPPORTED_PACK_FAMILIES
    assert "ancestry" in SUPPORTED_PACK_FAMILIES
    assert pack_zip_path(data_root, "bmi").name == "bmi_methylation_v1.zip"
    assert pack_zip_path(data_root, "ancestry").name == "ancestry_category_methylation_v1.zip"
    with pytest.raises(ValueError, match="unsupported"):
        pack_zip_path(data_root, "not_a_family")


def test_stream_pack_betas_column_order(isolated_workspace: Path) -> None:
    data_root = isolated_workspace / "data"
    zip_path, _ = _write_tiny_pack(data_root)
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
    _write_tiny_pack(data_root)
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
    locus = read_locus_index(out / "locus_index.parquet")
    assert "contributing_probe_ids" in locus.columns
    assert "collapse_method" in locus.columns


def test_content_checksum_not_name_size(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    zip_path, _ = _write_tiny_pack(data_root)
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations)
    out = data_root / "canonical" / "matrices" / "matrix-hub-age-hash-v1"
    convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="age",
        study_ids=["GSE_TRAIN", "GSE_VAL"],
        matrix_id="matrix-hub-age-hash-v1",
        output_dir=out,
    )
    manifest = json.loads((out / "matrix_manifest.json").read_text())
    entry = manifest["source_files"][0]
    assert entry["sha256"] == sha256_file(zip_path)
    fake = hashlib.sha256(
        f"{zip_path.name}:{zip_path.stat().st_size}:age_methylation_v1.txt".encode()
    ).hexdigest()
    assert entry["sha256"] != fake


def test_duplicate_gsm_long_form_no_overwrite(isolated_workspace: Path) -> None:
    info = pd.DataFrame(
        {
            "sample_id": ["GSM_A", "GSM_A", "GSM_B"],
            "study_id": ["GSE1", "GSE1", "GSE1"],
            "platform": ["450K", "450K", "450K"],
            "phenotype_value": ["diabetes", "obesity", "control"],
            "sample_type": ["case", "case", "control"],
            "phenotype_family": ["disease", "disease", "disease"],
        }
    )
    unique, long_form = select_samples_for_studies(info, ["GSE1"])
    assert len(unique) == 2
    assert len(long_form) == 3
    conflict = info.copy()
    conflict.loc[1, "platform"] = "EPIC"
    with pytest.raises(ValueError, match="conflicting platform"):
        select_samples_for_studies(conflict, ["GSE1"])
    assert long_form["phenotype_value"].tolist() == ["diabetes", "obesity", "control"]

    repo = _repo_root()
    data_root = isolated_workspace / "data"
    lines = [
        "sample_id\tGSM_A\tGSM_B\n",
        "disease\tdiabetes\tcontrol\n",
        "cg00000001\t0.10\t0.20\n",
        "cg00000002\t0.11\t0.21\n",
        "cg00000003\t0.12\t0.22\n",
    ]
    _write_tiny_pack(data_root, family="disease", lines=lines, info=info)
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations)
    out = data_root / "canonical" / "matrices" / "matrix-hub-disease-dup-v1"
    result = convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="disease",
        study_ids=["GSE1"],
        matrix_id="matrix-hub-disease-dup-v1",
        output_dir=out,
    )
    assert result.stats["n_samples"] == 2
    assert result.stats["n_phenotype_rows"] == 3
    pheno = pd.read_parquet(out / "sample_phenotypes.parquet")
    assert len(pheno) == 3
    assert set(pheno.loc[pheno["sample_id"] == "GSM_A", "phenotype_value"]) == {
        "diabetes",
        "obesity",
    }


def test_probe_collapse_mean_and_median() -> None:
    edges = pd.DataFrame(
        {
            "probe_id": ["cg001", "cg002", "cg003", "cg004"],
            "locus_id": [10, 10, 20, 20],
            "canonical_key": ["chr1:10", "chr1:10", "chr1:20", "chr1:20"],
            "genome_build": ["GRCh38"] * 4,
            "platform_id": ["HM450"] * 4,
            "is_primary": [True] * 4,
        }
    )
    # two probes → mean; three at locus 20 needs a third
    edges_three = pd.concat(
        [
            edges,
            pd.DataFrame(
                {
                    "probe_id": ["cg005"],
                    "locus_id": [20],
                    "canonical_key": ["chr1:20"],
                    "genome_build": ["GRCh38"],
                    "platform_id": ["HM450"],
                    "is_primary": [True],
                }
            ),
        ],
        ignore_index=True,
    )
    m2 = build_probe_locus_map(["cg001", "cg002"], edges, platform_id="HM450")
    assert m2.n_collapsed_probes == 1
    assert m2.collapse_method[0] == COLLAPSE_MEAN
    assert m2.contributing_probe_ids[0] == ("cg001", "cg002")

    m3 = build_probe_locus_map(["cg003", "cg004", "cg005"], edges_three, platform_id="HM450")
    assert m3.collapse_method[0] == COLLAPSE_MEDIAN
    assert m3.contributing_probe_ids[0] == ("cg003", "cg004", "cg005")


def test_chunked_convert_matches_dense_oracle(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    _write_tiny_pack(data_root)
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations)
    out = data_root / "canonical" / "matrices" / "matrix-hub-age-oracle-v1"
    convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="age",
        study_ids=["GSE_TRAIN", "GSE_VAL"],
        matrix_id="matrix-hub-age-oracle-v1",
        output_dir=out,
    )
    zip_path = pack_zip_path(data_root, "age")
    dense, probes, _ = stream_pack_betas(
        zip_path=zip_path, family="age", sample_ids=["GSM_A", "GSM_B", "GSM_C"]
    )
    # Map dense probe matrix through identity locus map (1:1 in mini annotations).
    got = np.asarray(open_betas_zarr(out / "betas.zarr")[:])
    # Columns follow locus_id order matching probe order in this fixture.
    np.testing.assert_allclose(got, dense, rtol=0, atol=1e-6, equal_nan=True)
    assert list(probes) == ["cg00000001", "cg00000002", "cg00000003"]


def test_collapse_mean_values_in_convert(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    lines = [
        "sample_id\tGSM_A\n",
        "age\t40\n",
        "cg00000001\t0.20\n",
        "cg00000002\t0.40\n",
        "cg00000003\t0.50\n",
    ]
    info = pd.DataFrame(
        {
            "sample_id": ["GSM_A"],
            "study_id": ["GSE1"],
            "platform": ["450K"],
            "phenotype_value": ["40"],
            "phenotype_value_numeric": [40.0],
            "phenotype_family": ["age"],
        }
    )
    _write_tiny_pack(data_root, lines=lines, info=info)
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations, collapse=True)
    out = data_root / "canonical" / "matrices" / "matrix-hub-age-collapse-v1"
    convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="age",
        study_ids=["GSE1"],
        matrix_id="matrix-hub-age-collapse-v1",
        output_dir=out,
    )
    betas = np.asarray(open_betas_zarr(out / "betas.zarr")[:])
    # locus 101 = mean(0.20, 0.40); locus 102 = 0.50
    np.testing.assert_allclose(betas[0], [0.30, 0.50], rtol=0, atol=1e-6)
    locus = read_locus_index(out / "locus_index.parquet")
    assert locus.loc[0, "contributing_probe_ids"] == "cg00000001|cg00000002"
    assert locus.loc[0, "collapse_method"] == "mean"


def test_collapse_median_values_in_convert(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    lines = [
        "sample_id\tGSM_A\n",
        "age\t40\n",
        "cg00000001\t0.10\n",
        "cg00000002\t0.20\n",
        "cg00000003\t0.90\n",
    ]
    info = pd.DataFrame(
        {
            "sample_id": ["GSM_A"],
            "study_id": ["GSE1"],
            "platform": ["450K"],
            "phenotype_value": ["40"],
            "phenotype_value_numeric": [40.0],
            "phenotype_family": ["age"],
        }
    )
    _write_tiny_pack(data_root, lines=lines, info=info)
    annotations = data_root / "canonical" / "annotations"
    annotations.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "locus_id": np.array([101], dtype=np.uint64),
            "genome_build": ["GRCh38"],
            "chromosome": ["chr1"],
            "position": [100],
            "canonical_key": ["GRCh38:chr1:100"],
            "mapping_status": ["mapped"],
            "cpg_context": ["island"],
        }
    ).to_parquet(annotations / "loci.parquet", index=False)
    pd.DataFrame(
        {
            "probe_id": ["cg00000001", "cg00000002", "cg00000003"],
            "platform_id": ["HM450", "HM450", "HM450"],
            "locus_id": np.array([101, 101, 101], dtype=np.uint64),
            "is_primary": [True, True, True],
        }
    ).to_parquet(annotations / "probe_locus_edges.parquet", index=False)
    out = data_root / "canonical" / "matrices" / "matrix-hub-age-median-v1"
    convert_hub_pack_subset(
        project_root=repo,
        data_root=data_root,
        annotations_dir=annotations,
        phenotype_family="age",
        study_ids=["GSE1"],
        matrix_id="matrix-hub-age-median-v1",
        output_dir=out,
    )
    betas = np.asarray(open_betas_zarr(out / "betas.zarr")[:])
    np.testing.assert_allclose(betas[0], [0.20], rtol=0, atol=1e-6)
    locus = read_locus_index(out / "locus_index.parquet")
    assert locus.loc[0, "collapse_method"] == "median"
    assert locus.loc[0, "contributing_probe_ids"] == "cg00000001|cg00000002|cg00000003"


def test_overlap_concordance_and_discordant(isolated_workspace: Path) -> None:
    repo = _repo_root()
    data_root = isolated_workspace / "data"
    annotations = data_root / "canonical" / "annotations"
    _write_mini_annotations(annotations)

    def _make(family: str, matrix_id: str, val_a: float) -> None:
        lines = [
            "sample_id\tGSM_A\tGSM_B\n",
            f"{family}\tx\ty\n",
            f"cg00000001\t{val_a}\t0.20\n",
            "cg00000002\t0.11\t0.21\n",
            "cg00000003\t0.12\t0.22\n",
        ]
        info = pd.DataFrame(
            {
                "sample_id": ["GSM_A", "GSM_B"],
                "study_id": ["GSE1", "GSE1"],
                "platform": ["450K", "450K"],
                "phenotype_value": ["x", "y"],
                "phenotype_family": [family, family],
            }
        )
        _write_tiny_pack(data_root, family=family, lines=lines, info=info)
        convert_hub_pack_subset(
            project_root=repo,
            data_root=data_root,
            annotations_dir=annotations,
            phenotype_family=family,
            study_ids=["GSE1"],
            matrix_id=matrix_id,
            output_dir=data_root / "canonical" / "matrices" / matrix_id,
        )

    _make("age", "matrix-hub-age-full-v1", 0.10)
    _make("tissue", "matrix-hub-tissue-full-v1", 0.10)
    index = build_hub_pack_matrix_index(data_root)
    assert set(index["family"]) == {"age", "tissue"}
    ok = check_overlapping_gsm_betas(data_root, index=index)
    assert ok.n_discordant == 0
    assert ok.report["merge_allowed"] is True

    # Overwrite tissue pack with discordant GSM_A
    _make("sex", "matrix-hub-sex-full-v1", 0.99)
    index2 = build_hub_pack_matrix_index(data_root)
    bad = check_overlapping_gsm_betas(data_root, index=index2, tolerance=1e-4)
    assert bad.n_discordant >= 1
    assert bad.report["merge_allowed"] is False
