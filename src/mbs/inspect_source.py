"""Shallow source-directory inventory for Stage 0 inspection reports."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TypedDict


class FileEntry(TypedDict):
    relative_path: str
    size_bytes: int
    suffix: str
    modified_at: str


class SourceInventory(TypedDict):
    source_id: str
    raw_root: str
    exists: bool
    file_count: int
    total_bytes: int
    truncated: bool
    max_entries: int
    generated_at: str
    files: list[FileEntry]
    suffix_counts: dict[str, int]


def inventory_source(
    raw_root: Path,
    *,
    source_id: str,
    max_entries: int = 10_000,
) -> SourceInventory:
    """Walk ``raw_root`` for file metadata only (no binary matrix reads)."""
    raw_root = raw_root.absolute()
    generated_at = datetime.now(UTC).isoformat()
    files: list[FileEntry] = []
    suffix_counts: dict[str, int] = {}
    total_bytes = 0
    truncated = False

    if not raw_root.exists():
        return {
            "source_id": source_id,
            "raw_root": str(raw_root),
            "exists": False,
            "file_count": 0,
            "total_bytes": 0,
            "truncated": False,
            "max_entries": max_entries,
            "generated_at": generated_at,
            "files": [],
            "suffix_counts": {},
        }

    for path in sorted(raw_root.rglob("*")):
        if not path.is_file():
            continue
        if len(files) >= max_entries:
            truncated = True
            break
        try:
            stat = path.stat()
        except OSError:
            continue
        relative = path.relative_to(raw_root).as_posix()
        suffix = path.suffix.lower() or "(none)"
        suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1
        total_bytes += int(stat.st_size)
        files.append(
            {
                "relative_path": relative,
                "size_bytes": int(stat.st_size),
                "suffix": suffix,
                "modified_at": datetime.fromtimestamp(stat.st_mtime, tz=UTC).isoformat(),
            }
        )

    return {
        "source_id": source_id,
        "raw_root": str(raw_root),
        "exists": True,
        "file_count": len(files),
        "total_bytes": total_bytes,
        "truncated": truncated,
        "max_entries": max_entries,
        "generated_at": generated_at,
        "files": files,
        "suffix_counts": dict(sorted(suffix_counts.items())),
    }


def _summary_markdown(inventory: SourceInventory) -> str:
    lines = [
        f"# Source inspection: `{inventory['source_id']}`",
        "",
        f"- Generated at: `{inventory['generated_at']}`",
        f"- Raw root: `{inventory['raw_root']}`",
        f"- Exists: `{inventory['exists']}`",
        f"- File count: `{inventory['file_count']}`",
        f"- Total bytes: `{inventory['total_bytes']}`",
        f"- Truncated: `{inventory['truncated']}` (max entries `{inventory['max_entries']}`)",
        "",
        "## Suffix counts",
        "",
    ]
    if not inventory["suffix_counts"]:
        lines.append("_No files inventoried._")
    else:
        lines.append("| Suffix | Count |")
        lines.append("| --- | ---: |")
        for suffix, count in inventory["suffix_counts"].items():
            lines.append(f"| `{suffix}` | {count} |")

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This report is a shallow file inventory only.",
            "- Binary matrices are not opened; checksums are not computed for large files.",
            "- Re-run after downloads before converting a source into a canonical release.",
            "",
        ]
    )
    return "\n".join(lines)


def write_inspection_report(
    inventory: SourceInventory,
    report_dir: Path,
) -> Path:
    """Write ``summary.json`` and ``summary.md`` under ``report_dir``."""
    report_dir = report_dir.absolute()
    report_dir.mkdir(parents=True, exist_ok=True)
    summary_json = report_dir / "summary.json"
    summary_md = report_dir / "summary.md"
    payload: dict[str, Any] = dict(inventory)
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary_md.write_text(_summary_markdown(inventory), encoding="utf-8")
    return report_dir
