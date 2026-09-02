"""EWAS_db download log parsing and retry manifest helpers."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

WARN_RE = re.compile(r"^WARN: failed (?P<study>[^/]+)/(?P<file>.+)$")
STUDY_PROGRESS_RE = re.compile(r"^\[(?P<current>\d+)/(?P<total>\d+)\] (?P<study>\S+)$")
PARSE_ARTIFACT_NAMES = frozenset({"(.+?)", "(.+)"})


def parse_ewas_db_download_log(log_path: Path) -> tuple[dict[str, list[str]], dict[str, int | None]]:
    """Parse study→failed files and last progress line from an EWAS_db wget log."""
    by_study: dict[str, list[str]] = defaultdict(list)
    progress: dict[str, int | None] = {"last_study_index": None, "advertised_studies": None}
    if not log_path.is_file():
        return dict(by_study), progress
    with log_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            prog = STUDY_PROGRESS_RE.match(line)
            if prog:
                progress["last_study_index"] = int(prog.group("current"))
                progress["advertised_studies"] = int(prog.group("total"))
                continue
            match = WARN_RE.match(line)
            if not match:
                continue
            by_study[match.group("study")].append(match.group("file"))
    return dict(by_study), progress


def ewas_db_failures_still_missing(
    ewas_db_root: Path,
    by_study: dict[str, list[str]],
) -> dict[str, list[str]]:
    """Return failed files that are still absent or zero-length under EWAS_db."""
    still_missing: dict[str, list[str]] = {}
    for study, files in sorted(by_study.items()):
        missing: list[str] = []
        study_dir = ewas_db_root / study
        for name in files:
            path = study_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(name)
        if missing:
            still_missing[study] = missing
    return still_missing
