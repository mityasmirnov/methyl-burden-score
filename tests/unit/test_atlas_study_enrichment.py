"""Unit tests for EWAS Atlas study-level enrichment."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

from mbs.atlas_study_enrichment import (
    build_study_atlas_enrichment,
    gds_uid_for_gse,
    load_atlas_gse_map,
    merge_atlas_enrichment_into_studies,
    resolve_atlas_study_ids,
    write_study_atlas_enrichment_report,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _fixture_atlas(tmp_path: Path) -> Path:
    src = _repo_root() / "tests/fixtures/ewas_metadata/atlas"
    dest = tmp_path / "ewas_atlas"
    shutil.copytree(src, dest)
    return dest


def test_gds_uid_for_gse() -> None:
    assert gds_uid_for_gse("GSE100197") == "200100197"
    assert gds_uid_for_gse("GSE80970") == "200080970"


def test_load_atlas_gse_map_ignores_comments(tmp_path: Path) -> None:
    path = tmp_path / "map.tsv"
    path.write_text(
        "\n".join(
            [
                "# comment",
                "gse_id\tatlas_study_id\tpmid\tsource",
                "GSE1\tES1\t123\tmanual",
            ]
        ),
        encoding="utf-8",
    )
    frame = load_atlas_gse_map(path)
    assert list(frame["gse_id"]) == ["GSE1"]


def test_resolve_atlas_study_ids_gse_es_map() -> None:
    gse_map = pd.DataFrame(
        [{"gse_id": "GSE9", "atlas_study_id": "ES1", "pmid": "123", "source": "manual"}]
    )
    atlas_ids, method, pmid = resolve_atlas_study_ids(
        catalog_study_id="GSE9",
        gse_map=gse_map,
        pmid_to_study={},
    )
    assert atlas_ids == ["ES1"]
    assert method == "gse_es_map"
    assert pmid == "123"


def test_resolve_atlas_study_ids_pmid_bridge() -> None:
    gse_map = pd.DataFrame([{"gse_id": "GSE9", "atlas_study_id": "", "pmid": "456", "source": "geo"}])
    atlas_ids, method, pmid = resolve_atlas_study_ids(
        catalog_study_id="GSE9",
        gse_map=gse_map,
        pmid_to_study={"456": ["ES2"]},
    )
    assert atlas_ids == ["ES2"]
    assert method == "pmid"
    assert pmid == "456"


def test_build_study_atlas_enrichment_summarizes_cohorts(tmp_path: Path) -> None:
    atlas_root = _fixture_atlas(tmp_path)
    map_path = tmp_path / "map.tsv"
    map_path.write_text(
        "gse_id\tatlas_study_id\tpmid\tsource\nGSE_HUB\tES1\t123\tmanual\n",
        encoding="utf-8",
    )
    frame = build_study_atlas_enrichment(
        catalog_study_ids=["GSE_HUB", "GSE_MISSING"],
        atlas_root=atlas_root,
        gse_map_path=map_path,
    )
    matched = frame.loc[frame["study_id"] == "GSE_HUB"].iloc[0]
    assert matched["join_method"] == "gse_es_map"
    assert json.loads(str(matched["atlas_study_ids"])) == ["ES1"]
    assert matched["n_atlas_cohorts"] == 1
    assert matched["total_sample_size"] == 10
    assert "blood" in json.loads(str(matched["tissues"]))
    missing = frame.loc[frame["study_id"] == "GSE_MISSING"].iloc[0]
    assert missing["join_method"] == "none"


def test_merge_atlas_enrichment_into_studies() -> None:
    studies = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_db"]}),
            }
        ]
    )
    enrichment = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "join_method": "gse_es_map",
                "atlas_study_ids": json.dumps(["ES1"]),
                "pmid": "123",
                "n_atlas_cohorts": 1,
                "total_sample_size": 10,
                "tissues": json.dumps(["blood"]),
                "cohort_descriptions": json.dumps(["desc"]),
                "platforms": json.dumps(["450K"]),
                "ancestries": json.dumps(["European"]),
                "atlas_traits": json.dumps(["age"]),
            }
        ]
    )
    merged = merge_atlas_enrichment_into_studies(studies, enrichment)
    meta = json.loads(str(merged.iloc[0]["metadata_json"]))
    assert meta["lanes"] == ["ewas_datahub_db"]
    assert meta["atlas_enrichment"]["join_method"] == "gse_es_map"
    assert meta["atlas_enrichment"]["tissues"] == ["blood"]


def test_write_study_atlas_enrichment_report(tmp_path: Path) -> None:
    enrichment = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "join_method": "gse_es_map",
                "atlas_study_ids": '["ES1"]',
                "pmid": "123",
                "n_atlas_cohorts": 1,
                "total_sample_size": 10,
                "tissues": '["blood"]',
                "cohort_descriptions": '["desc"]',
                "platforms": '["450K"]',
                "ancestries": '["European"]',
                "atlas_traits": '["age"]',
            }
        ]
    )
    out = write_study_atlas_enrichment_report(enrichment=enrichment, report_dir=tmp_path)
    assert out.is_file()
    assert (tmp_path / "study_atlas_enrichment.md").is_file()
