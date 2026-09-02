#!/usr/bin/env python3
"""Summarize EWAS_db wget failures from download log for retry + inspection."""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

from mbs.ewas_download import (
    PARSE_ARTIFACT_NAMES,
    ewas_db_failures_still_missing,
    parse_ewas_db_download_log,
)
from mbs.paths import DataPaths


def write_manifest(manifest_path: Path, still_missing: dict[str, list[str]]) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["study_id\tfilename"]
    for study, files in sorted(still_missing.items()):
        for name in sorted(set(files)):
            lines.append(f"{study}\t{name}")
    manifest_path.write_text("\n".join(lines) + ("\n" if len(lines) > 1 else ""), encoding="utf-8")


def write_markdown(
    report_path: Path,
    *,
    log_path: Path,
    by_study: dict[str, list[str]],
    still_missing: dict[str, list[str]],
    progress: dict[str, int | None],
    ewas_db_root: Path,
) -> None:
    n_failures = sum(len(v) for v in by_study.values())
    n_studies = len(by_study)
    n_missing = sum(len(v) for v in still_missing.values())
    parse_artifacts = sum(
        1 for files in by_study.values() for f in files if f in PARSE_ARTIFACT_NAMES
    )
    top = sorted(
        ((study, len(files)) for study, files in by_study.items()),
        key=lambda item: (-item[1], item[0]),
    )[:25]
    generated = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "# EWAS_db download failures",
        "",
        f"- Generated: `{generated}`",
        f"- Log: `{log_path}`",
        f"- EWAS_db root: `{ewas_db_root}`",
        f"- Studies with ≥1 logged failure: **{n_studies}**",
        f"- Total `WARN: failed` lines: **{n_failures}**",
        f"- Still missing or empty on disk: **{n_missing}**",
        f"- HTML-parse artifact filenames (`(.+?)`): **{parse_artifacts}**",
    ]
    if progress.get("last_study_index") is not None:
        lines.append(
            f"- Last study progress in log: **{progress['last_study_index']}**"
            f" / **{progress['advertised_studies']}**"
        )
    lines.extend(
        [
            "",
            "## Retry",
            "",
            "Manifest: `artifacts/logs/downloads/ewas_db_retry_manifest.tsv`",
            "",
            "```bash",
            "bash scripts/retry_ewas_db_download_failures.sh",
            "```",
            "",
            "## Top studies by failure count",
            "",
            "| study_id | failures | still_missing |",
            "| --- | ---: | ---: |",
        ]
    )
    for study, count in top:
        miss = len(still_missing.get(study, []))
        lines.append(f"| `{study}` | {count} | {miss} |")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Re-run `bash scripts/download_ewas_datahub.sh EWAS_db` to resume; "
            "successful files are skipped via `wget -c`.",
            "- Post-download hook runs `mbs catalog refresh-release` automatically "
            "(disable with `EWAS_DATAHUB_SKIP_POST_HOOK=1`).",
            "",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    paths = DataPaths.from_environment()
    default_log = paths.artifact_root / "logs" / "downloads" / "ewas_datahub_EWAS_db.log"
    default_report = paths.project_root / "reports" / "inspection" / "deepmat_data_v1"
    default_manifest = paths.artifact_root / "logs" / "downloads" / "ewas_db_retry_manifest.tsv"
    parser.add_argument("--log", type=Path, default=default_log)
    parser.add_argument(
        "--ewas-db-root",
        type=Path,
        default=paths.data_root / "raw" / "ewas_datahub" / "EWAS_db",
    )
    parser.add_argument("--report-dir", type=Path, default=default_report)
    parser.add_argument("--manifest-out", type=Path, default=default_manifest)
    args = parser.parse_args()

    by_study, progress = parse_ewas_db_download_log(args.log)
    still_missing = ewas_db_failures_still_missing(args.ewas_db_root, by_study)
    write_manifest(args.manifest_out, still_missing)

    summary = {
        "generated_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "log_path": str(args.log),
        "ewas_db_root": str(args.ewas_db_root),
        "n_studies_with_failures": len(by_study),
        "n_logged_failures": sum(len(v) for v in by_study.values()),
        "n_still_missing": sum(len(v) for v in still_missing.values()),
        "n_parse_artifact_filenames": sum(
            1 for files in by_study.values() for f in files if f in PARSE_ARTIFACT_NAMES
        ),
        "progress": progress,
        "manifest_path": str(args.manifest_out),
        "top_studies": [
            {"study_id": study, "n_failures": count}
            for study, count in sorted(
                ((s, len(f)) for s, f in by_study.items()),
                key=lambda item: (-item[1], item[0]),
            )[:50]
        ],
        "by_study": {study: sorted(set(files)) for study, files in sorted(by_study.items())},
        "still_missing": {
            study: sorted(set(files)) for study, files in sorted(still_missing.items())
        },
    }
    json_path = args.report_dir / "ewas_db_download_failures.json"
    md_path = args.report_dir / "ewas_db_download_failures.md"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_markdown(
        md_path,
        log_path=args.log,
        by_study=by_study,
        still_missing=still_missing,
        progress=progress,
        ewas_db_root=args.ewas_db_root,
    )
    print(f"wrote {json_path}", flush=True)
    print(f"wrote {md_path}", flush=True)
    print(f"wrote {args.manifest_out}", flush=True)
    print(
        f"failures={summary['n_logged_failures']} "
        f"still_missing={summary['n_still_missing']} "
        f"studies={summary['n_studies_with_failures']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
