"""Unit tests for GEO sample metadata backfill."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mbs.geo_metadata import (
    GEO_SOURCE_FAMILY,
    catalog_platform_from_gpl,
    characteristics_to_phenotypes,
    consolidate_geo_sample_rows,
    family_soft_url,
    map_geo_tissue,
    merge_geo_sample_metadata,
    parse_family_soft,
)
from mbs.training.phenotype_table import TissueOntology


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _fixture_soft() -> str:
    return (_repo_root() / "tests/fixtures/geo/GSE_FIXTURE_family.soft").read_text(encoding="utf-8")


def _mini_ontology() -> TissueOntology:
    labels = ("colon", "whole blood", "peripheral blood mononuclear cell")
    return TissueOntology(
        labels=labels,
        min_n=1,
        label_to_id={lab: i for i, lab in enumerate(labels)},
    )


def test_family_soft_url() -> None:
    url = family_soft_url("GSE197678")
    assert url.endswith("GSE197678/soft/GSE197678_family.soft.gz")
    assert "GSE197nnn" in url


def test_catalog_platform_from_gpl() -> None:
    assert catalog_platform_from_gpl("GPL13534") == "HM450"
    assert catalog_platform_from_gpl("GPL21145") == "EPIC"
    assert catalog_platform_from_gpl("GPL23976") == "EPICv2"
    assert catalog_platform_from_gpl("GPL6883") is None
    assert catalog_platform_from_gpl("GPL29753") is None


def test_parse_age_units() -> None:
    rows_years = characteristics_to_phenotypes({"age": "45"})
    assert rows_years[0]["numeric_value"] == 45.0
    assert rows_years[0]["age_raw"] == "45"
    assert rows_years[0]["age_unit"] == "years"

    rows_mo = characteristics_to_phenotypes({"age": "6 months"})
    assert abs(rows_mo[0]["numeric_value"] - 0.5) < 1e-9
    assert rows_mo[0]["age_raw"] == "6 months"
    assert rows_mo[0]["age_unit"] == "months"

    rows_d = characteristics_to_phenotypes({"age": "120 days"})
    assert abs(rows_d[0]["numeric_value"] - (120.0 / 365.25)) < 1e-9
    assert rows_d[0]["age_unit"] == "days"

    # Must not treat "6 months" as 6 years.
    assert characteristics_to_phenotypes({"age": "6 months"})[0]["numeric_value"] != 6.0
    assert characteristics_to_phenotypes({"age": "approx six"}) == []


def test_map_geo_tissue_aliases() -> None:
    ont = _mini_ontology()
    aliases = {
        "blood": "whole blood",
        "peripheral blood": "whole blood",
        "whole blood": "whole blood",
    }
    label, oid, status = map_geo_tissue("blood", ontology=ont, aliases=aliases)
    assert status == "mapped"
    assert label == "whole blood"
    assert oid == "1"
    label2, _, status2 = map_geo_tissue("peripheral blood", ontology=ont, aliases=aliases)
    assert status2 == "mapped" and label2 == "whole blood"
    _, _, status3 = map_geo_tissue("mystery organ", ontology=ont, aliases=aliases)
    assert status3 == "unmapped"


def test_parse_family_soft_samples() -> None:
    series, samples = parse_family_soft(_fixture_soft())
    assert series["study_id"] == "GSE_FIXTURE"
    assert series["pubmed_ids"] == ["12345678"]
    by_id = {s["sample_id"]: s for s in samples}
    assert set(by_id) == {"GSM_HUB", "GSM001", "GSM002"}
    assert by_id["GSM001"]["age"] == 45.0
    assert by_id["GSM001"]["age_raw"] == "45"
    assert by_id["GSM001"]["sex"] == "Female"
    assert by_id["GSM001"]["tissue"] == "whole blood"
    assert by_id["GSM001"]["disease_label_status"] == "control"
    assert "disease" not in by_id["GSM002"] or by_id["GSM002"].get("disease") is None
    assert by_id["GSM002"].get("age") is None


def test_characteristics_to_phenotypes_eligibility() -> None:
    rows = characteristics_to_phenotypes(
        {
            "age": "10002",
            "sex": "ambiguous",
            "tissue": "colon",
            "disease status": "Crohn's disease",
            "batch": "3",
        }
    )
    pheno_ids = {r["phenotype_id"] for r in rows}
    assert "age" not in pheno_ids
    assert "disease" not in pheno_ids
    assert "tissue" in pheno_ids
    sex_row = next(r for r in rows if r["phenotype_id"] == "sex")
    assert sex_row["label_status"] == "unknown"
    assert sex_row["is_observed"] is False


def test_characteristics_control_disease() -> None:
    rows = characteristics_to_phenotypes({"disease status": "control"})
    disease = next(r for r in rows if r["phenotype_id"] == "disease")
    assert disease["label_status"] == "control"
    assert disease["is_observed"] is True


def test_consolidate_geo_agrees_and_conflicts() -> None:
    agree = pd.DataFrame(
        [
            {
                "sample_id": "GSM_X",
                "study_id": "GSE1",
                "age": 40.0,
                "sex": "Female",
                "tissue": "whole blood",
                "platform_id": "GPL21145",
            },
            {
                "sample_id": "GSM_X",
                "study_id": "GSE2",
                "age": 40.0,
                "sex": "Female",
                "tissue": "whole blood",
                "platform_id": "GPL21145",
            },
        ]
    )
    out, stats = consolidate_geo_sample_rows(agree)
    assert len(out) == 1
    assert stats["n_conflict_samples"] == 0
    assert stats["n_multi_study_gsm"] == 1
    assert json.loads(out.iloc[0]["study_ids"]) == ["GSE1", "GSE2"]

    conflict = pd.DataFrame(
        [
            {
                "sample_id": "GSM_Y",
                "study_id": "GSE1",
                "age": 40.0,
                "sex": "Female",
                "tissue": "whole blood",
                "platform_id": "GPL21145",
            },
            {
                "sample_id": "GSM_Y",
                "study_id": "GSE2",
                "age": 41.0,
                "sex": "Male",
                "tissue": "whole blood",
                "platform_id": "GPL21145",
            },
        ]
    )
    out2, stats2 = consolidate_geo_sample_rows(conflict)
    assert len(out2) == 1
    assert stats2["n_conflict_samples"] == 1
    row = out2.iloc[0]
    assert row["age"] is None or (isinstance(row["age"], float) and pd.isna(row["age"]))
    assert row["sex"] is None or (isinstance(row["sex"], float) and pd.isna(row["sex"]))
    assert row["tissue"] == "whole blood"


def test_merge_geo_skips_hub_gsm() -> None:
    samples = pd.DataFrame(
        [
            {
                "sample_id": "GSM_HUB",
                "study_id": "GSE_FIXTURE",
                "source_sample_id": "GSM_HUB",
                "donor_id": None,
                "replicate_group": None,
                "age": 40.0,
                "sex": "Male",
                "tissue_raw": "blood",
                "tissue_ontology_id": None,
                "case_control": None,
                "metadata_json": None,
            },
            {
                "sample_id": "GSM001",
                "study_id": "GSE_FIXTURE",
                "source_sample_id": "GSM001",
                "donor_id": None,
                "replicate_group": None,
                "age": None,
                "sex": None,
                "tissue_raw": None,
                "tissue_ontology_id": None,
                "case_control": None,
                "metadata_json": json.dumps({"source": "ewas_db"}),
            },
        ]
    )
    studies = pd.DataFrame(
        [
            {
                "study_id": "GSE_FIXTURE",
                "source_release_id": "ewas-datahub-db-v1",
                "gse_id": "GSE_FIXTURE",
                "cohort_id": None,
                "platform_id": None,
                "processing_level": "raw_beta_txt",
                "genome_build": "GRCh38",
                "retrieved_at": "2026-01-01T00:00:00Z",
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_db"]}),
            }
        ]
    )
    _, geo_samples = parse_family_soft(_fixture_soft())
    geo_frame = pd.DataFrame(
        [
            {
                "sample_id": s["sample_id"],
                "study_id": s["study_id"],
                "study_ids": json.dumps([s["study_id"]]),
                "source_name": s.get("source_name"),
                "platform_id": s.get("platform_id"),
                "catalog_platform_id": s.get("catalog_platform_id"),
                "pubmed_ids": json.dumps(s.get("pubmed_ids") or []),
                "characteristics_raw": json.dumps(s.get("characteristics_raw") or {}),
                "age": s.get("age"),
                "age_raw": s.get("age_raw"),
                "age_unit": s.get("age_unit"),
                "sex": s.get("sex"),
                "sex_label_status": s.get("sex_label_status"),
                "tissue_raw": s.get("tissue_raw") or s.get("tissue"),
                "tissue": s.get("tissue"),
                "tissue_ontology_id": None,
                "tissue_map_status": None,
                "disease": s.get("disease"),
                "disease_label_status": s.get("disease_label_status"),
                "cancer": s.get("cancer"),
                "cancer_label_status": s.get("cancer_label_status"),
                "fetched_at": "2026-01-01T00:00:00Z",
                "soft_sha256": "abc",
            }
            for s in geo_samples
        ]
    )
    phenotypes = pd.DataFrame()
    out_samples, out_pheno, out_studies, stats = merge_geo_sample_metadata(
        samples=samples,
        phenotypes=phenotypes,
        studies=studies,
        geo_frame=geo_frame,
        ontology=_mini_ontology(),
        aliases={"blood": "whole blood", "whole blood": "whole blood"},
    )
    assert stats["n_samples_skipped_hub"] == 1
    assert stats["n_samples_touched"] == 1
    assert stats["n_geo_phenotype_rows_before_merge"] == 0
    assert set(out_pheno["sample_id"].tolist()) == {"GSM001"}
    assert (out_pheno["source_family"] == GEO_SOURCE_FAMILY).all()
    tissue_row = out_pheno.loc[out_pheno["phenotype_id"] == "tissue"].iloc[0]
    assert tissue_row["categorical_value"] == "whole blood"
    assert tissue_row["ontology_id"] == "1"
    hub_meta = out_samples.loc[out_samples["sample_id"] == "GSM_HUB", "metadata_json"].iloc[0]
    assert hub_meta is None or (isinstance(hub_meta, float) and pd.isna(hub_meta))
    ewas_meta = json.loads(
        out_samples.loc[out_samples["sample_id"] == "GSM001", "metadata_json"].iloc[0]
    )
    assert "geo" in ewas_meta
    assert ewas_meta["geo"]["source_name"] == "peripheral blood"
    assert ewas_meta["geo"]["age_raw"] == "45"
    assert ewas_meta["geo"]["tissue_raw"] == "whole blood"
    study_meta = json.loads(out_studies.iloc[0]["metadata_json"])
    assert study_meta["lanes"] == ["ewas_datahub_db"]
    assert study_meta["geo"]["pubmed_ids"] == ["12345678"]
    assert out_studies.iloc[0]["platform_id"] == "EPIC"


def test_merge_leaves_mixed_gpl_study_platform_null() -> None:
    samples = pd.DataFrame(
        [
            {
                "sample_id": "GSM_A",
                "study_id": "GSE_MIX",
                "source_sample_id": "GSM_A",
                "donor_id": None,
                "replicate_group": None,
                "age": None,
                "sex": None,
                "tissue_raw": None,
                "tissue_ontology_id": None,
                "case_control": None,
                "metadata_json": json.dumps({"source": "ewas_db"}),
            },
            {
                "sample_id": "GSM_B",
                "study_id": "GSE_MIX",
                "source_sample_id": "GSM_B",
                "donor_id": None,
                "replicate_group": None,
                "age": None,
                "sex": None,
                "tissue_raw": None,
                "tissue_ontology_id": None,
                "case_control": None,
                "metadata_json": json.dumps({"source": "ewas_db"}),
            },
        ]
    )
    studies = pd.DataFrame(
        [
            {
                "study_id": "GSE_MIX",
                "source_release_id": "ewas-datahub-db-v1",
                "gse_id": "GSE_MIX",
                "cohort_id": None,
                "platform_id": None,
                "processing_level": "raw_beta_txt",
                "genome_build": "GRCh38",
                "retrieved_at": "2026-01-01T00:00:00Z",
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_db"]}),
            }
        ]
    )
    geo_frame = pd.DataFrame(
        [
            {
                "sample_id": "GSM_A",
                "study_id": "GSE_MIX",
                "source_name": "a",
                "platform_id": "GPL13534",
                "catalog_platform_id": None,
                "pubmed_ids": "[]",
                "characteristics_raw": "{}",
                "age": None,
                "sex": "Female",
                "sex_label_status": "observed",
                "tissue": None,
                "disease": None,
                "disease_label_status": None,
                "cancer": None,
                "cancer_label_status": None,
                "fetched_at": "2026-01-01T00:00:00Z",
                "soft_sha256": "x",
            },
            {
                "sample_id": "GSM_B",
                "study_id": "GSE_MIX",
                "source_name": "b",
                "platform_id": "GPL21145",
                "catalog_platform_id": None,
                "pubmed_ids": "[]",
                "characteristics_raw": "{}",
                "age": None,
                "sex": "Male",
                "sex_label_status": "observed",
                "tissue": None,
                "disease": None,
                "disease_label_status": None,
                "cancer": None,
                "cancer_label_status": None,
                "fetched_at": "2026-01-01T00:00:00Z",
                "soft_sha256": "x",
            },
        ]
    )
    _, _, out_studies, _ = merge_geo_sample_metadata(
        samples=samples,
        phenotypes=pd.DataFrame(),
        studies=studies,
        geo_frame=geo_frame,
    )
    plat = out_studies.iloc[0]["platform_id"]
    assert plat is None or (isinstance(plat, float) and pd.isna(plat))
