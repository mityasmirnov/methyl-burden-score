"""EWAS_db download failure log parsing."""

from __future__ import annotations

from pathlib import Path

from mbs.ewas_download import (
    ewas_db_failures_still_missing,
    is_retryable_ewas_db_filename,
    parse_ewas_db_download_log,
)


def test_parse_ewas_db_download_log(tmp_path: Path) -> None:
    log = tmp_path / "ewas.log"
    log.write_text(
        "\n".join(
            [
                "[1/1989] GSE100",
                "WARN: failed GSE100/GSM1.txt",
                "WARN: failed GSE100/GSM2.txt",
                "[2/1989] GSE200",
                "WARN: failed GSE200/(.+?)",
            ]
        ),
        encoding="utf-8",
    )
    by_study, progress = parse_ewas_db_download_log(log)
    assert progress["last_study_index"] == 2
    assert progress["advertised_studies"] == 1989
    assert by_study["GSE100"] == ["GSM1.txt", "GSM2.txt"]
    assert by_study["GSE200"] == ["(.+?)"]


def test_is_retryable_ewas_db_filename() -> None:
    assert is_retryable_ewas_db_filename("GSM123456.txt")
    assert not is_retryable_ewas_db_filename("(.+?)")
    assert not is_retryable_ewas_db_filename("(.+)")
    assert not is_retryable_ewas_db_filename("present.txt")
    assert not is_retryable_ewas_db_filename("index.html")


def test_ewas_db_failures_still_missing(tmp_path: Path) -> None:
    root = tmp_path / "EWAS_db"
    study = root / "GSE1"
    study.mkdir(parents=True)
    (study / "GSM1.txt").write_text("ok\n", encoding="utf-8")
    (study / "GSM2.txt").write_text("", encoding="utf-8")
    by_study = {
        "GSE1": ["GSM1.txt", "GSM2.txt", "GSM3.txt"],
        "GSE2": ["(.+?)", "GSM9.txt"],
    }
    missing = ewas_db_failures_still_missing(root, by_study)
    assert missing["GSE1"] == ["GSM2.txt", "GSM3.txt"]
    assert missing["GSE2"] == ["GSM9.txt"]


def test_ewas_db_failures_still_missing_skips_artifacts_only(tmp_path: Path) -> None:
    root = tmp_path / "EWAS_db"
    (root / "GSE200").mkdir(parents=True)
    by_study = {"GSE200": ["(.+?)"]}
    assert ewas_db_failures_still_missing(root, by_study) == {}
