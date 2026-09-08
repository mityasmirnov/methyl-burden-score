"""GEO series metadata parsing and study merge."""

from __future__ import annotations

import json

import pandas as pd

from mbs.geo_series_metadata import (
    merge_geo_series_into_studies,
    parse_series_soft_header,
)


def test_parse_series_soft_header_stops_at_sample() -> None:
    soft = (
        "^SERIES = GSE1\n"
        "!Series_geo_accession = GSE1\n"
        "!Series_title = Example methylation study\n"
        "!Series_summary = Short summary.\n"
        "!Series_overall_design = Case vs control blood.\n"
        "!Series_type = Methylation profiling by array\n"
        "!Series_pubmed_id = 111\n"
        "!Series_platform_id = GPL13534\n"
        "^SAMPLE = GSM1\n"
        "!Sample_geo_accession = GSM1\n"
        "!Sample_characteristics_ch1 = age: 40\n"
    )
    series = parse_series_soft_header(soft)
    assert series["study_id"] == "GSE1"
    assert series["title"] == "Example methylation study"
    assert series["summary"] == "Short summary."
    assert series["overall_design"] == "Case vs control blood."
    assert series["series_type"] == "Methylation profiling by array"
    assert series["pubmed_ids"] == ["111"]
    assert series["platform_ids"] == ["GPL13534"]


def test_merge_geo_series_into_studies() -> None:
    studies = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "metadata_json": json.dumps({"source": "ewas_db", "geo": {"pubmed_ids": ["999"]}}),
            },
            {"study_id": "GSE2", "metadata_json": None},
        ]
    )
    series = pd.DataFrame(
        [
            {
                "study_id": "GSE1",
                "title": "Title One",
                "summary": "Sum",
                "overall_design": "Design",
                "series_type": "Methylation",
                "pubmed_ids": json.dumps(["111"]),
                "platform_ids": json.dumps(["GPL13534"]),
                "n_samples_ncbi": 10,
                "taxon": "Homo sapiens",
                "platform_title": "Illumina 450K",
                "source": "ncbi_esummary+soft_header",
                "fetched_at": "2026-09-07T00:00:00Z",
            }
        ]
    )
    out = merge_geo_series_into_studies(studies, series)
    meta = json.loads(out.loc[out.study_id == "GSE1", "metadata_json"].iloc[0])
    assert meta["geo"]["title"] == "Title One"
    assert meta["geo"]["overall_design"] == "Design"
    assert meta["geo"]["pubmed_ids"] == ["111", "999"]
    # Untouched study stays null / unchanged.
    assert out.loc[out.study_id == "GSE2", "metadata_json"].iloc[0] is None or pd.isna(
        out.loc[out.study_id == "GSE2", "metadata_json"].iloc[0]
    )
