"""Unit tests for sample overview platform / species enrichment."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mbs.sample_overview_enrich import (
    fill_null_study_platforms_from_ewas_db,
    fill_species_for_rows,
    infer_study_platform_from_assay_dir,
    platform_from_n_probes,
    platform_from_sample_id_hint,
    species_from_series_taxon,
)


def test_platform_from_n_probes_bands() -> None:
    assert platform_from_n_probes(485512) == "HM450"
    assert platform_from_n_probes(866836) == "EPIC"
    assert platform_from_n_probes(937690) == "EPICv2"
    assert platform_from_n_probes(17000) is None


def test_fill_null_study_platforms_from_ewas_db(tmp_path: Path) -> None:
    study_dir = tmp_path / "GSE_TEST"
    study_dir.mkdir()
    # 485511 newline-terminated rows → HM450 band (compact payload)
    (study_dir / "GSM1.txt").write_bytes(b"x\n" * 485511)
    studies = pd.DataFrame(
        [
            {"study_id": "GSE_TEST", "platform_id": None},
            {"study_id": "GSE_KEEP", "platform_id": "EPIC"},
            {"study_id": "GSE_EMPTY", "platform_id": None},
        ]
    )
    (tmp_path / "GSE_EMPTY").mkdir()
    out, stats = fill_null_study_platforms_from_ewas_db(studies, ewas_db_root=tmp_path)
    by_id = out.set_index("study_id")["platform_id"].to_dict()
    assert by_id["GSE_TEST"] == "HM450"
    assert by_id["GSE_KEEP"] == "EPIC"
    assert by_id["GSE_EMPTY"] is None or (isinstance(by_id["GSE_EMPTY"], float) and pd.isna(by_id["GSE_EMPTY"]))
    assert stats["n_studies_platform_set"] == 1
    assert infer_study_platform_from_assay_dir(study_dir) == "HM450"


def test_platform_from_sample_id_935k() -> None:
    assert platform_from_sample_id_hint("abc_935k") == "EPICv2"
    assert platform_from_sample_id_hint("GSM1") is None


def test_species_from_series_taxon_single_only() -> None:
    status, org = species_from_series_taxon("Homo sapiens")
    assert status == "human"
    assert org == "Homo sapiens"
    status2, _ = species_from_series_taxon("Macaca mulatta; Homo sapiens")
    assert status2 is None
    status3, _ = species_from_series_taxon("Mus musculus")
    assert status3 == "non_human"


def test_fill_species_covers_all_rows() -> None:
    frame = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c", "d"],
            "study_id": ["GSE1", "TCGA-BRCA", "GSE2", "GSE3"],
            "in_hub_baseline": [True, False, False, False],
            "species_status": ["human", None, None, None],
            "series_taxon": [None, None, "Homo sapiens", None],
            "taxon_id": [9606, None, None, None],
            "organism": ["Homo sapiens", None, None, None],
        }
    )
    out = fill_species_for_rows(frame)
    assert out["species_status"].notna().all()
    assert set(out["species_status"]) <= {"human", "non_human"}
    assert out.loc[0, "species_source"] == "geo_soft"
    assert out.loc[1, "species_source"] == "project_namespace_human"
    assert out.loc[2, "species_source"] == "geo_series_taxon"
    assert out.loc[3, "species_source"] == "ewas_datahub_assumed_human"
