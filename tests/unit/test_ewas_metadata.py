from __future__ import annotations

from pathlib import Path

import pandas as pd

from mbs.inspect_ewas_metadata import (
    inspect_ewas_metadata,
    profile_atlas_table,
    profile_sample_pack,
    write_ewas_metadata_report,
)
from mbs.registry.sample_info import (
    export_family_from_data_root,
    export_sample_info_parquet,
    read_r_style_table,
    resolve_sample_info_txt,
)


def _fixtures() -> Path:
    return Path(__file__).resolve().parents[1] / "fixtures" / "ewas_metadata"


def test_read_r_style_table_and_export_from_txt(tmp_path: Path) -> None:
    txt = _fixtures() / "samples" / "sample_age_methylation_v1" / "sample_age.txt"
    frame = read_r_style_table(txt)
    assert list(frame.columns)[:3] == ["sample_id", "age", "project_id"]
    assert len(frame) == 2

    out = tmp_path / "age_sample_info.parquet"
    export_sample_info_parquet(txt_path=txt, family="age", output_path=out)
    exported = pd.read_parquet(out)
    assert "phenotype_value" in exported.columns
    assert "phenotype_value_numeric" in exported.columns
    assert set(exported["sample_id"]) == {"GSM1", "GSM2"}


def test_resolve_prefers_unpacked_over_missing_zip(tmp_path: Path) -> None:
    project = tmp_path / "proj"
    data_root = project / "data"
    samples = (
        project
        / "reports"
        / "inspection"
        / "ewas_datahub_samples"
        / "sample_age_methylation_v1"
    )
    samples.mkdir(parents=True)
    src = _fixtures() / "samples" / "sample_age_methylation_v1" / "sample_age.txt"
    (samples / "sample_age.txt").write_text(src.read_text(encoding="utf-8"), encoding="utf-8")

    path, source = resolve_sample_info_txt(
        data_root=data_root,
        family="age",
        project_root=project,
    )
    assert source == "unpacked"
    assert path.name == "sample_age.txt"

    out = export_family_from_data_root(data_root, "age", project_root=project)
    assert out.is_file()
    assert pd.read_parquet(out)["phenotype_family"].iloc[0] == "age"


def test_inspect_ewas_metadata_on_fixtures(tmp_path: Path) -> None:
    project = tmp_path / "proj"
    data_root = project / "data"
    atlas = data_root / "raw" / "ewas_atlas"
    atlas.mkdir(parents=True)
    fx = _fixtures()
    for name in (
        "EWAS_Atlas_studies.tsv",
        "EWAS_Atlas_cohorts.tsv",
        "EWAS_trait_trait_logP.txt",
    ):
        (atlas / name).write_text((fx / "atlas" / name).read_text(encoding="utf-8"), encoding="utf-8")

    samples = (
        project
        / "reports"
        / "inspection"
        / "ewas_datahub_samples"
        / "sample_age_methylation_v1"
    )
    samples.mkdir(parents=True)
    (samples / "sample_age.txt").write_text(
        (fx / "samples" / "sample_age_methylation_v1" / "sample_age.txt").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )

    report = inspect_ewas_metadata(
        data_root=data_root,
        project_root=project,
        sample_families=("age",),
    )
    assert sum(1 for t in report["atlas_tables"] if t["exists"]) == 3
    age = next(p for p in report["sample_packs"] if p["family"] == "age")
    assert age["exists"] is True
    assert age["n_rows"] == 2
    assert age["primary_phenotype_column"] == "age"

    written = write_ewas_metadata_report(report, tmp_path / "report")
    assert (written / "summary.md").is_file()
    assert (written / "summary.json").is_file()
    assert "EWAS metadata structure" in (written / "summary.md").read_text(encoding="utf-8")


def test_profile_helpers_missing_paths(tmp_path: Path) -> None:
    missing = profile_atlas_table(tmp_path / "nope.tsv", table_id="studies", fmt="tsv")
    assert missing["exists"] is False
    missing_pack = profile_sample_pack(tmp_path / "nope.txt", family="age")
    assert missing_pack["exists"] is False
