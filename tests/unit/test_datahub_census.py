"""Unit tests for DataHub census + Hub/EWAS_db lane flags."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mbs.datahub_census import (
    build_sample_lane_flags,
    build_study_lane_flags,
    catalog_platform_from_datahub,
    fetch_datahub_repository_census,
    merge_datahub_census,
    stamp_lane_flags_into_metadata,
    write_datahub_census_report,
    write_lane_flags_report,
)
from mbs.platform_id import normalize_platform


def test_catalog_platform_from_datahub() -> None:
    assert catalog_platform_from_datahub("450K") == "HM450"
    assert catalog_platform_from_datahub("850K") == "EPIC"
    assert catalog_platform_from_datahub("935K") == "EPICv2"
    assert normalize_platform("850K") == "EPIC"
    assert normalize_platform("935K") == "EPICv2"


def test_build_sample_and_study_lane_flags() -> None:
    samples = pd.DataFrame(
        [
            {"sample_id": "GSM1", "study_id": "GSE1"},
            {"sample_id": "GSM2", "study_id": "GSE1"},
            {"sample_id": "GSM3", "study_id": "GSE2"},
        ]
    )
    membership = pd.DataFrame(
        [
            {"sample_id": "GSM1", "phenotype_family": "age"},
            {"sample_id": "GSM1", "phenotype_family": "tissue"},
        ]
    )
    ewas_files = pd.DataFrame(
        [
            {"sample_id": "GSM2", "study_id": "GSE1"},
            {"sample_id": "GSM3", "study_id": "GSE2"},
        ]
    )
    flags = build_sample_lane_flags(
        samples=samples, membership=membership, ewas_files=ewas_files
    )
    by = flags.set_index("sample_id")
    assert bool(by.loc["GSM1", "in_hub_baseline"]) is True
    assert bool(by.loc["GSM1", "in_ewas_db"]) is False
    assert json.loads(str(by.loc["GSM1", "hub_families"])) == ["age", "tissue"]
    assert bool(by.loc["GSM2", "in_hub_baseline"]) is False
    assert bool(by.loc["GSM2", "in_ewas_db"]) is True
    assert bool(by.loc["GSM3", "in_ewas_db"]) is True

    studies = pd.DataFrame([{"study_id": "GSE1"}, {"study_id": "GSE2"}])
    study_flags = build_study_lane_flags(studies=studies, sample_flags=flags)
    sby = study_flags.set_index("study_id")
    assert bool(sby.loc["GSE1", "in_hub_baseline"]) is True
    assert bool(sby.loc["GSE1", "in_ewas_db"]) is True
    assert int(sby.loc["GSE1", "n_hub_samples"]) == 1
    assert int(sby.loc["GSE1", "n_ewas_db_samples"]) == 1
    assert bool(sby.loc["GSE2", "in_hub_baseline"]) is False


def test_stamp_lane_flags_into_metadata() -> None:
    samples = pd.DataFrame(
        [{"sample_id": "GSM1", "study_id": "GSE1", "metadata_json": None}]
    )
    studies = pd.DataFrame(
        [{"study_id": "GSE1", "metadata_json": json.dumps({"lanes": ["other"]})}]
    )
    sample_flags = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE1",
                "in_hub_baseline": True,
                "in_ewas_db": True,
                "hub_families": "[]",
            }
        ]
    )
    study_flags = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "in_hub_baseline": True,
                "in_ewas_db": True,
                "n_hub_samples": 1,
                "n_ewas_db_samples": 1,
                "n_samples": 1,
            }
        ]
    )
    out_s, out_st = stamp_lane_flags_into_metadata(
        samples=samples,
        studies=studies,
        sample_flags=sample_flags,
        study_flags=study_flags,
    )
    s_meta = json.loads(str(out_s.iloc[0]["metadata_json"]))
    assert s_meta["lanes"] == ["ewas_datahub_baseline", "ewas_datahub_db"]
    st_meta = json.loads(str(out_st.iloc[0]["metadata_json"]))
    assert "ewas_datahub_baseline" in st_meta["lanes"]
    assert "ewas_datahub_db" in st_meta["lanes"]
    assert "other" in st_meta["lanes"]


def test_merge_datahub_census_fills_platform_and_hub_wins(tmp_path: Path) -> None:
    samples = pd.DataFrame(
        [
            {
                "sample_id": "GSM_HUB",
                "study_id": "GSE1",
                "tissue_raw": "blood",
                "case_control": None,
                "metadata_json": None,
            },
            {
                "sample_id": "GSM_DB",
                "study_id": "GSE2",
                "tissue_raw": None,
                "case_control": None,
                "metadata_json": None,
            },
        ]
    )
    phenotypes = pd.DataFrame(
        columns=[
            "sample_id",
            "phenotype_id",
            "numeric_value",
            "categorical_value",
            "label_status",
            "is_observed",
            "source_family",
            "source_record_id",
            "ontology_id",
        ]
    )
    studies = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "platform_id": None,
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_baseline"]}),
            },
            {
                "study_id": "GSE2",
                "platform_id": None,
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_db"]}),
            },
        ]
    )
    census = pd.DataFrame(
        [
            {
                "sample_id": "GSM_HUB",
                "project_id": "GSE1",
                "platform": "450K",
                "catalog_platform_id": "HM450",
                "tissue": "liver",
                "sample_type": "control",
                "disease": "should_not_add",
            },
            {
                "sample_id": "GSM_DB",
                "project_id": "GSE2",
                "platform": "850K",
                "catalog_platform_id": "EPIC",
                "tissue": "brain",
                "sample_type": "disease tissue",
                "disease": "Alzheimer",
                "sex": "female",
            },
        ]
    )
    out_s, out_p, out_st, stats = merge_datahub_census(
        samples=samples,
        phenotypes=phenotypes,
        studies=studies,
        census=census,
        hub_sample_ids={"GSM_HUB"},
    )
    assert stats["n_samples_matched"] == 2
    hub = out_s.loc[out_s["sample_id"] == "GSM_HUB"].iloc[0]
    assert hub["tissue_raw"] == "blood"  # Hub tissue kept
    db = out_s.loc[out_s["sample_id"] == "GSM_DB"].iloc[0]
    assert db["tissue_raw"] == "brain"
    assert db["case_control"] == "case"
    meta = json.loads(str(db["metadata_json"]))
    assert meta["datahub"]["platform"] == "850K"
    # Hub sample must not get repository phenotype rows
    assert not ((out_p["sample_id"] == "GSM_HUB") & (out_p["source_family"] == "ewas_datahub_repository")).any()
    assert (out_p["sample_id"] == "GSM_DB").any()
    gse2 = out_st.loc[out_st["study_id"] == "GSE2"].iloc[0]
    assert gse2["platform_id"] == "EPIC"


def test_fetch_census_from_cache_only(tmp_path: Path) -> None:
    cache = tmp_path / "cache" / "450K"
    cache.mkdir(parents=True)
    page = {
        "total": "2",
        "content": [
            {
                "sample id": "GSM_A",
                "project id": "GSE1",
                "platform": "450K",
                "tissue": "blood",
                "sample type": "control",
                "disease": "none",
            },
            {
                "sample id": "GSM_B",
                "project id": "GSE2",
                "platform": "450K",
                "tissue": "brain",
                "sample type": "disease tissue",
                "disease": "ALS",
            },
        ],
    }
    (cache / "page_000000.json").write_text(json.dumps(page), encoding="utf-8")
    out = tmp_path / "phenotypes"
    manifest = fetch_datahub_repository_census(
        cache_root=tmp_path / "cache",
        output_dir=out,
        platforms=("450K",),
        from_cache_only=True,
    )
    assert manifest["n_samples"] == 2
    assert (out / "ewas_datahub_sample_census.parquet").is_file()
    frame = pd.read_parquet(out / "ewas_datahub_sample_census.parquet")
    assert set(frame["sample_id"]) == {"GSM_A", "GSM_B"}
    assert set(frame["catalog_platform_id"]) == {"HM450"}


def test_write_reports(tmp_path: Path) -> None:
    sample_flags = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE1",
                "in_hub_baseline": True,
                "in_ewas_db": False,
                "hub_families": "[]",
            },
            {
                "sample_id": "GSM2",
                "study_id": "GSE2",
                "in_hub_baseline": False,
                "in_ewas_db": True,
                "hub_families": "[]",
            },
        ]
    )
    study_flags = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "in_hub_baseline": True,
                "in_ewas_db": False,
                "n_hub_samples": 1,
                "n_ewas_db_samples": 0,
                "n_samples": 1,
            },
            {
                "study_id": "GSE2",
                "in_hub_baseline": False,
                "in_ewas_db": True,
                "n_hub_samples": 0,
                "n_ewas_db_samples": 1,
                "n_samples": 1,
            },
        ]
    )
    assert write_lane_flags_report(
        sample_flags=sample_flags, study_flags=study_flags, report_dir=tmp_path
    ).is_file()
    census = pd.DataFrame(
        [
            {"sample_id": "GSM2", "platform": "450K", "tissue": "blood", "disease": "x"},
        ]
    )
    assert write_datahub_census_report(
        census=census,
        sample_flags=sample_flags,
        merge_stats={"enabled": False},
        report_dir=tmp_path,
        manifest={"n_fields": 10},
    ).is_file()
