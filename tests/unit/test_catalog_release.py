"""Unit tests for Milestone 7A catalog release refresh."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import duckdb
import pandas as pd
import pytest
from typer.testing import CliRunner

from mbs.cli import app
from mbs.paths import DataPaths
from mbs.release import (
    RELEASE_ID,
    compute_trait_eligibility,
    head_ranking_eligible,
    head_training_allowed,
    refresh_release,
    release_paths,
    scan_ewas_db_tree,
    validate_multitask_head_eligibility,
    validate_release,
    validate_release_manifest,
)

runner = CliRunner()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def release_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"release-{uuid4().hex}"
    workspace.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(workspace / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(workspace / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(workspace / "docker"))
    return workspace


def _write_gsm(path: Path, probe: str = "cg00000029", beta: str = "0.5") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{probe}\t{beta}\n", encoding="utf-8")


def _seed_hub_and_ewas(workspace: Path, *, with_third_study: bool = False) -> DataPaths:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    pheno = paths.data_root / "canonical" / "phenotypes"
    pheno.mkdir(parents=True, exist_ok=True)

    # Overlapping GSM across age + tissue → unique N < row sum
    age = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE_A",
                "platform": "450K",
                "sample_type": "sample",
                "sex": "Male",
                "phenotype_value": "40",
                "phenotype_value_numeric": 40.0,
                "phenotype_family": "age",
                "source_zip": "sample_age.txt",
            },
            {
                "sample_id": "GSM2",
                "study_id": "GSE_A",
                "platform": "450K",
                "sample_type": "sample",
                "sex": "Female",
                "phenotype_value": "50",
                "phenotype_value_numeric": 50.0,
                "phenotype_family": "age",
                "source_zip": "sample_age.txt",
            },
        ]
    )
    tissue = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE_A",
                "platform": "450K",
                "sample_type": "sample",
                "sex": "Male",
                "phenotype_value": "blood",
                "phenotype_family": "tissue",
                "source_zip": "sample_tissue.txt",
            },
            {
                "sample_id": "GSM3",
                "study_id": "GSE_B",
                "platform": "EPIC",
                "sample_type": "sample",
                "sex": "Female",
                "phenotype_value": "brain",
                "phenotype_family": "tissue",
                "source_zip": "sample_tissue.txt",
            },
        ]
    )
    disease = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE_A",
                "platform": "450K",
                "sample_type": "case",
                "phenotype_value": "asthma",
                "phenotype_family": "disease",
                "source_zip": "sample_disease.txt",
            },
            {
                "sample_id": "GSM4",
                "study_id": "GSE_C",
                "platform": "450K",
                "sample_type": "control",
                "phenotype_value": "control",
                "phenotype_family": "disease",
                "source_zip": "sample_disease.txt",
            },
            {
                "sample_id": "GSM5",
                "study_id": "GSE_C",
                "platform": "450K",
                "sample_type": "sample",
                "phenotype_value": None,
                "phenotype_family": "disease",
                "source_zip": "sample_disease.txt",
            },
        ]
    )
    age.to_parquet(pheno / "age_sample_info.parquet", index=False)
    tissue.to_parquet(pheno / "tissue_sample_info.parquet", index=False)
    disease.to_parquet(pheno / "disease_sample_info.parquet", index=False)

    ewas = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    _write_gsm(ewas / "GSE_A" / "GSM1.txt")
    _write_gsm(ewas / "GSE_A" / "GSM2.txt")
    _write_gsm(ewas / "GSE_B" / "GSM3.txt")
    if with_third_study:
        _write_gsm(ewas / "GSE_NEW" / "GSM9.txt")
    return paths


def test_unique_gsm_less_than_pack_row_sum(release_workspace: Path) -> None:
    paths = _seed_hub_and_ewas(release_workspace)
    result = refresh_release(paths=paths, report_dir=None)
    assert result.n_samples >= 3
    rp = release_paths(paths.data_root)
    connection = duckdb.connect(str(rp.catalog_db), read_only=True)
    try:
        unique_row = connection.execute(
            "SELECT count(DISTINCT sample_id) FROM sample"
        ).fetchone()
        pack_row = connection.execute(
            "SELECT count(*) FROM sample_source_membership"
        ).fetchone()
        assert unique_row is not None and pack_row is not None
        unique = unique_row[0]
        pack_sum = pack_row[0]
    finally:
        connection.close()
    assert int(unique) < int(pack_sum)


def test_incremental_ewas_db_third_study(release_workspace: Path) -> None:
    paths = _seed_hub_and_ewas(release_workspace, with_third_study=False)
    first = refresh_release(paths=paths, report_dir=None)
    n1 = first.ewas_db_n_local_studies
    _write_gsm(paths.data_root / "raw" / "ewas_datahub" / "EWAS_db" / "GSE_NEW" / "GSM9.txt")
    second = refresh_release(paths=paths, report_dir=None)
    assert second.ewas_db_n_local_studies == n1 + 1
    assert second.ewas_db_n_local_gsm == first.ewas_db_n_local_gsm + 1


def test_disease_unknown_not_control() -> None:
    phenotypes = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "phenotype_id": "disease",
                "numeric_value": None,
                "categorical_value": "asthma",
                "label_status": "case",
                "is_observed": True,
                "source_family": "disease",
                "source_record_id": "disease:0",
                "ontology_id": None,
            },
            {
                "sample_id": "GSM2",
                "phenotype_id": "disease",
                "numeric_value": None,
                "categorical_value": "control",
                "label_status": "control",
                "is_observed": True,
                "source_family": "disease",
                "source_record_id": "disease:1",
                "ontology_id": None,
            },
            {
                "sample_id": "GSM3",
                "phenotype_id": "disease",
                "numeric_value": None,
                "categorical_value": None,
                "label_status": "unknown",
                "is_observed": False,
                "source_family": "disease",
                "source_record_id": "disease:2",
                "ontology_id": None,
            },
        ]
    )
    samples = pd.DataFrame(
        [
            {"sample_id": "GSM1", "study_id": "G1", "platform_id": "HM450", "tissue_raw": None},
            {"sample_id": "GSM2", "study_id": "G2", "platform_id": "HM450", "tissue_raw": None},
            {"sample_id": "GSM3", "study_id": "G3", "platform_id": "HM450", "tissue_raw": None},
        ]
    )
    elig = compute_trait_eligibility(phenotypes, samples)
    row = elig[(elig["phenotype_id"] == "disease") & (elig["phenotype_family"] == "disease")].iloc[
        0
    ]
    assert int(row["n_cases"]) == 1
    assert int(row["n_controls"]) == 1
    assert int(row["n_unknown"]) == 1
    assert not bool(row["eligible_core_task"])


def _eligibility_frame(*rows: dict[str, object]) -> pd.DataFrame:
    defaults = {
        "phenotype_id": "disease",
        "phenotype_family": "disease",
        "eligible_core_task": False,
        "eligible_auxiliary_task": False,
    }
    return pd.DataFrame([{**defaults, **row} for row in rows])


def test_head_training_allowed_core_or_aux() -> None:
    elig = _eligibility_frame(
        {"eligible_core_task": True},
        {
            "phenotype_id": "cancer",
            "phenotype_family": "cancer",
            "eligible_auxiliary_task": True,
        },
    )
    assert head_training_allowed(elig, phenotype_id="disease", phenotype_family="disease")
    assert head_training_allowed(elig, phenotype_id="cancer", phenotype_family="cancer")
    assert not head_training_allowed(elig, phenotype_id="disease", phenotype_family="cancer")


def test_head_ranking_eligible_core_only() -> None:
    elig = _eligibility_frame({"eligible_auxiliary_task": True})
    assert head_training_allowed(elig, phenotype_id="disease", phenotype_family="disease")
    assert not head_ranking_eligible(elig, phenotype_id="disease", phenotype_family="disease")


def test_validate_multitask_head_eligibility_raises_when_both_false(
    release_workspace: Path,
) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    rp = release_paths(paths.data_root)
    rp.catalog_tables.mkdir(parents=True, exist_ok=True)
    frame = _eligibility_frame({})
    frame.to_parquet(rp.catalog_tables / "trait_eligibility.parquet", index=False)
    with pytest.raises(ValueError, match="disease head not eligible"):
        validate_multitask_head_eligibility(
            data_root=paths.data_root,
            disease_enabled=True,
            cancer_enabled=False,
        )


def test_validate_multitask_head_eligibility_allows_aux(release_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    rp = release_paths(paths.data_root)
    rp.catalog_tables.mkdir(parents=True, exist_ok=True)
    frame = _eligibility_frame({"eligible_auxiliary_task": True})
    frame.to_parquet(rp.catalog_tables / "trait_eligibility.parquet", index=False)
    meta = validate_multitask_head_eligibility(
        data_root=paths.data_root,
        disease_enabled=True,
        cancer_enabled=False,
    )
    assert meta["disease"]["training_allowed"] is True
    assert meta["disease"]["ranking_eligible"] is False


def test_validate_multitask_head_eligibility_skip_check() -> None:
    assert (
        validate_multitask_head_eligibility(
            data_root=Path("/nonexistent"),
            disease_enabled=True,
            cancer_enabled=True,
            skip_check=True,
        )
        == {}
    )


def test_validate_release_fails_without_manifest(release_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    with pytest.raises(FileNotFoundError):
        validate_release(data_root=paths.data_root, release_id=RELEASE_ID)


def test_validate_release_manifest_missing_key() -> None:
    with pytest.raises(ValueError, match="missing keys"):
        validate_release_manifest({"release_id": RELEASE_ID})


def test_scan_ewas_db_tree_lists_txt_only(release_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    root = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    _write_gsm(root / "GSE_X" / "GSM1.txt")
    (root / "GSE_X" / "notes.md").write_text("ignore", encoding="utf-8")
    (root / "GSE_EMPTY").mkdir(parents=True, exist_ok=True)
    frame = scan_ewas_db_tree(root)
    assert list(frame["sample_id"]) == ["GSM1"]
    assert frame.iloc[0]["sha256"]


def test_catalog_refresh_release_cli(release_workspace: Path) -> None:
    _seed_hub_and_ewas(release_workspace)
    report_dir = release_workspace / "reports" / "inspection" / "refresh-cli"
    result = runner.invoke(
        app,
        [
            "catalog",
            "refresh-release",
            "--no-fetch-remote-index",
            "--report-dir",
            str(report_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["release_id"] == RELEASE_ID
    assert payload["n_samples"] > 0
    assert payload["report_dir"] == str(report_dir)
    assert (report_dir / "census.json").is_file()
    assert (report_dir / "trait_eligibility.json").is_file()
    val = runner.invoke(app, ["catalog", "validate-release"])
    assert val.exit_code == 0, val.output


def test_catalog_census_and_eligibility_cli_temp_report_dir(release_workspace: Path) -> None:
    _seed_hub_and_ewas(release_workspace)
    report_dir = release_workspace / "reports" / "inspection" / "census-cli"
    refresh = runner.invoke(
        app,
        [
            "catalog",
            "refresh-release",
            "--no-fetch-remote-index",
            "--report-dir",
            str(report_dir),
        ],
    )
    assert refresh.exit_code == 0, refresh.output
    census_dir = release_workspace / "reports" / "inspection" / "census-only"
    census = runner.invoke(
        app,
        ["catalog", "phenotype-census", "--report-dir", str(census_dir)],
    )
    assert census.exit_code == 0, census.output
    census_payload = json.loads((census_dir / "census.json").read_text(encoding="utf-8"))
    assert "age_distribution_by_study" in census_payload
    assert "bmi_distribution_by_study" in census_payload
    assert census_payload["donor_replicate"].get("n_with_donor", 0) == 0
    elig_dir = release_workspace / "reports" / "inspection" / "elig-only"
    elig = runner.invoke(
        app,
        ["catalog", "trait-eligibility", "--report-dir", str(elig_dir)],
    )
    assert elig.exit_code == 0, elig.output
    assert (elig_dir / "trait_eligibility.json").is_file()


def test_platform_alias_450k_to_hm450_on_refresh(release_workspace: Path) -> None:
    _seed_hub_and_ewas(release_workspace)
    report_dir = release_workspace / "reports" / "inspection" / "platform-alias"
    result = runner.invoke(
        app,
        [
            "catalog",
            "refresh-release",
            "--no-fetch-remote-index",
            "--report-dir",
            str(report_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    paths = DataPaths.from_environment()
    rp = release_paths(paths.data_root, RELEASE_ID)
    con = duckdb.connect(str(rp.catalog_db), read_only=True)
    platforms = {
        row[0]
        for row in con.execute("SELECT DISTINCT platform_id FROM study").fetchall()
    }
    con.close()
    assert "450K" not in platforms
    assert "HM450" in platforms


def test_geo_backfill_merges_ewas_only_not_hub(release_workspace: Path) -> None:
    paths = _seed_hub_and_ewas(release_workspace)
    ewas = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    _write_gsm(ewas / "GSE_ONLY" / "GSM_ONLY.txt")
    geo_frame = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE_A",
                "source_name": "should-not-apply",
                "platform_id": "GPL21145",
                "catalog_platform_id": "EPIC",
                "pubmed_ids": "[]",
                "characteristics_raw": "{}",
                "age": 99.0,
                "sex": "Female",
                "sex_label_status": "observed",
                "tissue": "ignored",
                "disease": None,
                "disease_label_status": None,
                "cancer": None,
                "cancer_label_status": None,
                "fetched_at": "2026-01-01T00:00:00Z",
                "soft_sha256": "abc",
            },
            {
                "sample_id": "GSM_ONLY",
                "study_id": "GSE_ONLY",
                "source_name": "ewas-only blood",
                "platform_id": "GPL6883",
                "catalog_platform_id": "HM450",
                "pubmed_ids": "[\"999\"]",
                "characteristics_raw": "{\"tissue\": \"blood\"}",
                "age": 33.0,
                "sex": "Male",
                "sex_label_status": "observed",
                "tissue": "blood",
                "disease": None,
                "disease_label_status": None,
                "cancer": None,
                "cancer_label_status": None,
                "fetched_at": "2026-01-01T00:00:00Z",
                "soft_sha256": "def",
            },
        ]
    )
    geo_frame.to_parquet(
        paths.data_root / "canonical" / "phenotypes" / "geo_sample_metadata.parquet",
        index=False,
    )
    refresh_release(paths=paths, report_dir=None)
    rp = release_paths(paths.data_root)
    con = duckdb.connect(str(rp.catalog_db), read_only=True)
    try:
        hub_geo = con.execute(
            """
            SELECT count(*) FROM sample_phenotype
            WHERE sample_id = 'GSM1' AND source_family = 'geo_metadata_backfill'
            """
        ).fetchone()
        ewas_geo = con.execute(
            """
            SELECT count(*) FROM sample_phenotype
            WHERE sample_id = 'GSM_ONLY' AND source_family = 'geo_metadata_backfill'
            """
        ).fetchone()
    finally:
        con.close()
    assert hub_geo is not None and int(hub_geo[0]) == 0
    assert ewas_geo is not None and int(ewas_geo[0]) >= 1