#!/usr/bin/env python3
"""Write Milestone 7B convert progress JSON/MD and refresh the plan Progress block.

Safe to run while ``convert_hub_full_packs.sh`` is mid-flight. Does not walk
``$MBS_DATA_ROOT`` recursively — only known ``matrix-hub-*-full-v1`` paths.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from mbs.paths import DataPaths

SEVEN_B = ("ancestry", "bmi", "brain", "blood", "cancer", "disease")
FROZEN = ("age", "tissue", "sex")
PROGRESS_MARK_START = "<!-- 7B-PROGRESS-START -->"
PROGRESS_MARK_END = "<!-- 7B-PROGRESS-END -->"


def _utc_now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _dir_bytes(path: Path) -> int | None:
    if not path.exists():
        return None
    if path.is_file():
        return int(path.stat().st_size)
    try:
        out = subprocess.check_output(["du", "-sb", str(path)], text=True)
        return int(out.split()[0])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError):
        return None


def _infer_phase(root: Path, *, scratch: Path, betas: Path) -> str:
    """Best-effort convert phase from on-disk artifacts (no process attach)."""
    if (root / "conversion_stats.json").is_file():
        return "done"
    if scratch.is_file():
        return "stream_scratch"
    if (root / "sample_index.parquet").is_file() and betas.exists():
        # Indices written after stream; QC + zip sha256 remain.
        if not (root / "matrix_manifest.json").is_file():
            return "qc_or_checksum"
        return "finishing"
    if betas.exists():
        return "write_zarr"
    if root.is_dir():
        return "started"
    return "pending"


def _matrix_status(data_root: Path, family: str) -> dict[str, Any]:
    matrix_id = f"matrix-hub-{family}-full-v1"
    root = data_root / "canonical" / "matrices" / matrix_id
    stats_path = root / "conversion_stats.json"
    man_path = root / "matrix_manifest.json"
    scratch = root / ".betas_scratch.f32"
    betas = root / "betas.zarr"
    betas_bytes = _dir_bytes(betas) if betas.exists() else None
    scratch_bytes = int(scratch.stat().st_size) if scratch.is_file() else None
    if stats_path.is_file():
        stats = _load_json(stats_path) or {}
        man = _load_json(man_path) or {}
        return {
            "family": family,
            "matrix_id": matrix_id,
            "status": "done",
            "phase": "done",
            "n_samples": stats.get("n_samples"),
            "n_phenotype_rows": stats.get("n_phenotype_rows"),
            "n_loci": stats.get("n_study_loci"),
            "platform_id": man.get("platform_id"),
            "has_scratch": scratch.is_file(),
            "has_betas": betas.exists(),
            "betas_bytes": betas_bytes,
            "scratch_bytes": scratch_bytes,
        }
    if root.is_dir() and (
        betas.exists() or scratch.is_file() or (root / "sample_index.parquet").is_file()
    ):
        return {
            "family": family,
            "matrix_id": matrix_id,
            "status": "in_progress",
            "phase": _infer_phase(root, scratch=scratch, betas=betas),
            "n_samples": None,
            "n_phenotype_rows": None,
            "n_loci": None,
            "platform_id": None,
            "has_scratch": scratch.is_file(),
            "has_betas": betas.exists(),
            "betas_bytes": betas_bytes,
            "scratch_bytes": scratch_bytes,
        }
    return {
        "family": family,
        "matrix_id": matrix_id,
        "status": "pending",
        "phase": "pending",
        "n_samples": None,
        "n_phenotype_rows": None,
        "n_loci": None,
        "platform_id": None,
        "has_scratch": False,
        "has_betas": False,
        "betas_bytes": None,
        "scratch_bytes": None,
    }


def _running_convert_families() -> list[str]:
    try:
        out = subprocess.check_output(
            ["pgrep", "-af", "mbs matrix convert-pack"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    found: list[str] = []
    for line in out.splitlines():
        parts = line.split()
        for i, tok in enumerate(parts):
            if tok == "--phenotype-family" and i + 1 < len(parts):
                fam = parts[i + 1]
                if fam not in found:
                    found.append(fam)
                break
    return found


def build_progress(data_root: Path) -> dict[str, Any]:
    seven_b = [_matrix_status(data_root, f) for f in SEVEN_B]
    frozen = [_matrix_status(data_root, f) for f in FROZEN]
    n_done = sum(1 for r in seven_b if r["status"] == "done")
    n_prog = sum(1 for r in seven_b if r["status"] == "in_progress")
    running = _running_convert_families()
    index_path = data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
    report_dir = (
        Path(__file__).resolve().parents[1] / "reports" / "inspection" / "stage0_7b_hub_matrices"
    )
    summary = report_dir / "summary.json"
    return {
        "milestone": "7B",
        "updated_at": _utc_now(),
        "pid": os.getpid(),
        "seven_b_done": n_done,
        "seven_b_total": len(SEVEN_B),
        "seven_b_in_progress": n_prog,
        "all_seven_b_done": n_done == len(SEVEN_B),
        "running_convert_families": running,
        "virtual_index_present": index_path.is_file(),
        "inspection_summary_present": summary.is_file(),
        "seven_b_matrices": seven_b,
        "frozen_matrices": frozen,
        "next_pending": next((r["family"] for r in seven_b if r["status"] == "pending"), None),
        "remaining": _remaining(n_done, seven_b),
    }


def _remaining(n_done: int, seven_b: list[dict[str, Any]]) -> list[str]:
    pending = [r["family"] for r in seven_b if r["status"] != "done"]
    if n_done < len(SEVEN_B):
        packs = ", ".join(pending) if pending else "remaining packs"
        return [
            f"Finish convert for: {packs} (skip-if-exists resume)",
            "mbs matrix index-hub-packs --check-overlap",
            "scripts/write_stage0_7b_report.py",
            "mbs catalog refresh-release (7A pointers only)",
            "Required checks; mark TODO_PIPELINE 7B done with evidence",
        ]
    return [
        "Confirm index + overlap report",
        "scripts/write_stage0_7b_report.py",
        "mbs catalog refresh-release",
        "Required checks; mark TODO_PIPELINE 7B done",
    ]


def _fmt_bytes(n: int | None) -> str:
    if n is None:
        return "—"
    gib = n / (1024**3)
    if gib >= 0.1:
        return f"{gib:.1f} GiB"
    mib = n / (1024**2)
    return f"{mib:.0f} MiB"


def _markdown(progress: dict[str, Any]) -> str:
    lines = [
        "# Milestone 7B convert progress",
        "",
        f"Updated: `{progress['updated_at']}`",
        "",
        f"**7B packs:** {progress['seven_b_done']} / {progress['seven_b_total']} done"
        + (
            f"; in progress: {progress['seven_b_in_progress']}"
            if progress["seven_b_in_progress"]
            else ""
        ),
        "",
        "| Family | Status | Phase | Samples | Phenotype rows | Loci | betas.zarr |",
        "|--------|--------|-------|--------:|---------------:|-----:|-----------:|",
    ]
    for row in progress["seven_b_matrices"]:
        samples = row["n_samples"] if row["n_samples"] is not None else "—"
        pheno = row["n_phenotype_rows"] if row["n_phenotype_rows"] is not None else "—"
        loci = row["n_loci"] if row["n_loci"] is not None else "—"
        phase = row.get("phase") or row["status"]
        lines.append(
            f"| `{row['family']}` | `{row['status']}` | `{phase}` | {samples} | "
            f"{pheno} | {loci} | {_fmt_bytes(row.get('betas_bytes'))} |"
        )
    running = progress.get("running_convert_families") or []
    lines.extend(
        [
            "",
            f"- Running convert-pack: {', '.join(f'`{x}`' for x in running) or 'none'}",
            f"- Virtual index present: `{progress['virtual_index_present']}`",
            f"- Inspection summary present: `{progress['inspection_summary_present']}`",
            "",
            "Phases: `pending` → `stream_scratch` → `write_zarr` → "
            "`qc_or_checksum` → `done`.",
            "",
            "## Remaining",
            "",
        ]
    )
    for item in progress["remaining"]:
        lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines)


def _plan_progress_block(progress: dict[str, Any]) -> str:
    rows = []
    for row in progress["seven_b_matrices"]:
        mark = {"done": "done", "in_progress": "in progress", "pending": "pending"}[
            row["status"]
        ]
        detail = ""
        if row["status"] == "done":
            detail = f" `{row['n_samples']}×{row['n_loci']}`"
        elif row["status"] == "in_progress":
            phase = row.get("phase") or "in_progress"
            detail = f" (`{phase}`; {_fmt_bytes(row.get('betas_bytes'))})"
        rows.append(f"| `{row['family']}` | {mark}{detail} |")
    running = progress.get("running_convert_families") or []
    return "\n".join(
        [
            PROGRESS_MARK_START,
            "",
            f"_Auto-updated `{progress['updated_at']}` by `scripts/update_7b_convert_progress.py`._",
            "",
            f"**{progress['seven_b_done']}/{progress['seven_b_total']}** 7B packs done."
            + (f" Active: `{', '.join(running)}`." if running else ""),
            "",
            "| Family | Status |",
            "|--------|--------|",
            *rows,
            "",
            "Track live: `reports/inspection/stage0_7b_hub_matrices/progress.md`",
            "",
            PROGRESS_MARK_END,
        ]
    )


def _update_plan(plan_path: Path, progress: dict[str, Any]) -> None:
    text = plan_path.read_text(encoding="utf-8")
    block = _plan_progress_block(progress)
    if PROGRESS_MARK_START in text and PROGRESS_MARK_END in text:
        before, rest = text.split(PROGRESS_MARK_START, 1)
        _, after = rest.split(PROGRESS_MARK_END, 1)
        text = before + block + after
    else:
        # Insert after Status paragraph.
        needle = "Status: **in_progress**"
        idx = text.find(needle)
        if idx < 0:
            text = text.rstrip() + "\n\n## Progress\n\n" + block + "\n"
        else:
            # Find end of Status sentence block (blank line after).
            insert_at = text.find("\n\n", idx)
            if insert_at < 0:
                text = text + "\n\n## Progress\n\n" + block + "\n"
            else:
                text = (
                    text[: insert_at + 2]
                    + "## Progress\n\n"
                    + block
                    + "\n"
                    + text[insert_at + 2 :]
                )
    plan_path.write_text(text, encoding="utf-8")


def main() -> None:
    paths = DataPaths.from_environment()
    progress = build_progress(paths.data_root)
    report_dir = (
        Path(__file__).resolve().parents[1] / "reports" / "inspection" / "stage0_7b_hub_matrices"
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "progress.json").write_text(
        json.dumps(progress, indent=2) + "\n", encoding="utf-8"
    )
    (report_dir / "progress.md").write_text(_markdown(progress), encoding="utf-8")
    plan = Path(__file__).resolve().parents[1] / "docs" / "plans" / "milestone-7b-complete-hub-matrices.md"
    if plan.is_file():
        _update_plan(plan, progress)
    print(json.dumps({"updated_at": progress["updated_at"], "seven_b_done": progress["seven_b_done"]}))


if __name__ == "__main__":
    main()
