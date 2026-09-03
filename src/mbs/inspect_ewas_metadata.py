"""Profile Cursor-visible EWAS Atlas small tables and DataHub sample-info .txt packs."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.registry.sample_info import (
    FAMILY_VALUE_COLUMN,
    read_r_style_table,
    sample_txt_filename,
    unpacked_sample_info_dir,
)

ATLAS_SMALL_FILES: tuple[tuple[str, str, str], ...] = (
    ("studies", "EWAS_Atlas_studies.tsv", "tsv"),
    ("cohorts", "EWAS_Atlas_cohorts.tsv", "tsv"),
    ("trait_trait_logP", "EWAS_trait_trait_logP.txt", "tsv_matrix"),
)

SAMPLE_FAMILIES_PRESENT: tuple[str, ...] = (
    "age",
    "ancestry",
    "blood",
    "bmi",
    "brain",
    "cancer",
    "disease",
    "sex",
    "tissue",
)


def _non_null_rate(series: pd.Series) -> float:
    if len(series) == 0:
        return 0.0
    return float(series.notna().mean())


def _infer_kind(series: pd.Series) -> str:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().mean() >= 0.8 and series.notna().any():
        return "numeric"
    nunique = int(series.nunique(dropna=True))
    if nunique <= max(20, int(0.05 * max(len(series), 1))):
        return "categorical"
    return "string"


def _column_profiles(frame: pd.DataFrame, *, max_top: int = 8) -> list[dict[str, Any]]:
    profiles: list[dict[str, Any]] = []
    for name in frame.columns:
        series = frame[name]
        kind = _infer_kind(series)
        entry: dict[str, Any] = {
            "name": str(name),
            "kind": kind,
            "non_null_rate": round(_non_null_rate(series), 4),
            "n_unique": int(series.nunique(dropna=True)),
        }
        if kind == "numeric":
            numeric = pd.to_numeric(series, errors="coerce")
            entry["min"] = float(numeric.min()) if numeric.notna().any() else None
            entry["max"] = float(numeric.max()) if numeric.notna().any() else None
            entry["mean"] = float(numeric.mean()) if numeric.notna().any() else None
        elif kind == "categorical":
            counts = series.value_counts(dropna=True).head(max_top)
            entry["top_values"] = [
                {"value": str(idx), "count": int(count)} for idx, count in counts.items()
            ]
        profiles.append(entry)
    return profiles


def _example_rows(frame: pd.DataFrame, *, n: int = 2) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for _, row in frame.head(n).iterrows():
        rows.append({str(k): (None if pd.isna(v) else str(v)) for k, v in row.items()})
    return rows


def read_atlas_tsv(path: Path) -> pd.DataFrame:
    """Read Atlas TSV; tolerate rare malformed rows with extra tabs."""
    frame, _ = _read_atlas_tsv(path)
    return frame


def _read_atlas_tsv(path: Path) -> tuple[pd.DataFrame, int]:
    """Read Atlas TSV; tolerate rare malformed rows with extra tabs."""
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, encoding="latin-1")
        return frame, 0
    except pd.errors.ParserError:
        frame = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            encoding="latin-1",
            engine="python",
            on_bad_lines="warn",
        )
        # Count how many lines failed a strict field-count check.
        expected = None
        bad = 0
        with path.open(encoding="latin-1") as handle:
            for i, line in enumerate(handle):
                n_fields = line.rstrip("\n").count("\t") + 1
                if i == 0:
                    expected = n_fields
                    continue
                if expected is not None and n_fields != expected:
                    bad += 1
        return frame, bad


def profile_atlas_table(path: Path, *, table_id: str, fmt: str) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        return {
            "table_id": table_id,
            "path": str(path),
            "exists": False,
        }

    if fmt == "tsv_matrix":
        frame, bad_rows = _read_atlas_tsv(path)
        row_labels = frame.iloc[:, 0]
        value_cols = list(frame.columns[1:])
        # Sample a few numeric cells for range (avoid full melt on large grid).
        sample_vals = pd.to_numeric(frame.iloc[:50, 1:51].stack(), errors="coerce")
        return {
            "table_id": table_id,
            "path": str(path),
            "exists": True,
            "bytes": path.stat().st_size,
            "parse_recipe": "tab-separated latin-1; first column = trait name; remaining = square logP matrix",
            "n_rows": int(len(frame)),
            "n_cols": int(frame.shape[1]),
            "n_malformed_rows_skipped": bad_rows,
            "row_key": str(frame.columns[0]),
            "n_traits": int(len(value_cols)),
            "matrix_is_square": len(value_cols) == len(frame)
            and set(value_cols) == set(row_labels.astype(str)),
            "value_sample_min": float(sample_vals.min()) if sample_vals.notna().any() else None,
            "value_sample_max": float(sample_vals.max()) if sample_vals.notna().any() else None,
            "example_traits": [str(x) for x in row_labels.head(5).tolist()],
            "join_keys": ["trait (row/col labels; not Atlas study_ID)"],
        }

    frame, bad_rows = _read_atlas_tsv(path)
    return {
        "table_id": table_id,
        "path": str(path),
        "exists": True,
        "bytes": path.stat().st_size,
        "parse_recipe": "tab-separated TSV with header (latin-1); skip rare malformed rows",
        "n_rows": int(len(frame)),
        "n_cols": int(frame.shape[1]),
        "n_malformed_rows_skipped": bad_rows,
        "columns": _column_profiles(frame),
        "join_keys": [c for c in ("study_ID", "cohort_ID", "PMID") if c in frame.columns],
        "example_rows": _example_rows(frame),
    }


def profile_sample_pack(path: Path, *, family: str) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        return {
            "family": family,
            "path": str(path),
            "exists": False,
            "primary_phenotype_column": FAMILY_VALUE_COLUMN.get(family),
        }

    frame = read_r_style_table(path)
    value_col = FAMILY_VALUE_COLUMN.get(family)
    primary_stats: dict[str, Any] | None = None
    if value_col and value_col in frame.columns:
        series = frame[value_col]
        kind = _infer_kind(series)
        primary_stats = {
            "column": value_col,
            "kind": kind,
            "non_null_rate": round(_non_null_rate(series), 4),
            "n_unique": int(series.nunique(dropna=True)),
        }
        if kind == "numeric":
            numeric = pd.to_numeric(series, errors="coerce")
            primary_stats["min"] = float(numeric.min()) if numeric.notna().any() else None
            primary_stats["max"] = float(numeric.max()) if numeric.notna().any() else None
            primary_stats["mean"] = float(numeric.mean()) if numeric.notna().any() else None
        else:
            counts = series.value_counts(dropna=True).head(10)
            primary_stats["top_values"] = [
                {"value": str(idx), "count": int(count)} for idx, count in counts.items()
            ]

    return {
        "family": family,
        "path": str(path),
        "exists": True,
        "bytes": path.stat().st_size,
        "parse_recipe": "R write.table: space-separated, double-quoted fields, index_col=0 dropped",
        "n_rows": int(len(frame)),
        "n_cols": int(frame.shape[1]),
        "columns": [str(c) for c in frame.columns],
        "column_profiles": _column_profiles(frame),
        "join_keys": [c for c in ("sample_id", "project_id") if c in frame.columns],
        "primary_phenotype_column": value_col,
        "primary_phenotype_stats": primary_stats,
        "example_rows": _example_rows(frame.loc[:, [c for c in frame.columns if c in {
            "sample_id", "project_id", "platform", "sex", "tissue", "disease", "age",
            "cell_component", "sample_type", value_col or "",
        }]]),
        "sample_ids": set(frame["sample_id"].dropna().astype(str)) if "sample_id" in frame.columns else set(),
        "project_ids": set(frame["project_id"].dropna().astype(str)) if "project_id" in frame.columns else set(),
    }


def _cross_pack_analysis(sample_profiles: list[dict[str, Any]]) -> dict[str, Any]:
    present = [p for p in sample_profiles if p.get("exists")]
    if not present:
        return {"n_packs": 0}

    col_sets = [set(p.get("columns") or []) for p in present]
    shared = set.intersection(*col_sets) if col_sets else set()
    union = set.union(*col_sets) if col_sets else set()

    pairwise_sample: list[dict[str, Any]] = []
    pairwise_project: list[dict[str, Any]] = []
    for i, a in enumerate(present):
        for b in present[i + 1 :]:
            sa, sb = a.get("sample_ids") or set(), b.get("sample_ids") or set()
            pa, pb = a.get("project_ids") or set(), b.get("project_ids") or set()
            pairwise_sample.append(
                {
                    "a": a["family"],
                    "b": b["family"],
                    "n_shared_sample_id": len(sa & sb),
                    "n_a": len(sa),
                    "n_b": len(sb),
                }
            )
            pairwise_project.append(
                {
                    "a": a["family"],
                    "b": b["family"],
                    "n_shared_project_id": len(pa & pb),
                    "n_a": len(pa),
                    "n_b": len(pb),
                }
            )

    return {
        "n_packs": len(present),
        "n_shared_columns": len(shared),
        "shared_columns": sorted(shared),
        "n_union_columns": len(union),
        "pairwise_sample_id_overlap": pairwise_sample,
        "pairwise_project_id_overlap": pairwise_project,
    }


def _atlas_hub_project_overlap(
    studies_path: Path,
    sample_profiles: list[dict[str, Any]],
) -> dict[str, Any]:
    if not studies_path.is_file():
        return {"exists": False}
    studies = pd.read_csv(studies_path, sep="\t", dtype=str, encoding="latin-1")
    atlas_ids = set(studies["study_ID"].dropna().astype(str)) if "study_ID" in studies.columns else set()
    hub_projects: set[str] = set()
    for p in sample_profiles:
        hub_projects |= set(p.get("project_ids") or set())
    # Hub project_id is typically GSE*; Atlas study_ID is ES##### — rarely equal.
    gse_like = {x for x in hub_projects if str(x).upper().startswith("GSE")}
    return {
        "n_atlas_study_id": len(atlas_ids),
        "n_hub_project_id": len(hub_projects),
        "n_hub_gse_like": len(gse_like),
        "n_exact_study_id_equals_project_id": len(atlas_ids & hub_projects),
        "note": (
            "Atlas study_ID (ES*) and Hub project_id (usually GSE*) are different namespaces; "
            "do not join on raw equality. Use PMID / curated maps when needed."
        ),
    }


def inspect_ewas_metadata(
    *,
    data_root: Path,
    project_root: Path,
    sample_families: tuple[str, ...] = SAMPLE_FAMILIES_PRESENT,
) -> dict[str, Any]:
    """Build a structure profile for Atlas small tables + unpacked sample packs."""
    data_root = data_root.resolve()
    project_root = project_root.resolve()
    atlas_root = data_root / "raw" / "ewas_atlas"

    atlas_tables = [
        profile_atlas_table(atlas_root / name, table_id=table_id, fmt=fmt)
        for table_id, name, fmt in ATLAS_SMALL_FILES
    ]

    sample_packs: list[dict[str, Any]] = []
    for family in sample_families:
        path = unpacked_sample_info_dir(project_root, family) / sample_txt_filename(family)
        sample_packs.append(profile_sample_pack(path, family=family))

    cross = _cross_pack_analysis(sample_packs)
    atlas_hub = _atlas_hub_project_overlap(atlas_root / "EWAS_Atlas_studies.tsv", sample_packs)

    # Strip non-JSON-serializable sets before return
    serializable_packs = []
    for p in sample_packs:
        clean = {k: v for k, v in p.items() if k not in {"sample_ids", "project_ids"}}
        if p.get("exists"):
            clean["n_sample_id"] = len(p.get("sample_ids") or set())
            clean["n_project_id"] = len(p.get("project_ids") or set())
        serializable_packs.append(clean)

    return {
        "generated_at": datetime.now(UTC).isoformat(),
        "data_root": str(data_root),
        "project_root": str(project_root),
        "atlas_tables": atlas_tables,
        "sample_packs": serializable_packs,
        "cross_pack": cross,
        "atlas_hub_id_overlap": atlas_hub,
        "family_value_column": dict(FAMILY_VALUE_COLUMN),
    }


def _md_escape(text: str) -> str:
    return text.replace("|", "\\|")


def render_ewas_metadata_markdown(report: dict[str, Any]) -> str:
    lines: list[str] = [
        "# EWAS metadata structure",
        "",
        f"Generated at: `{report['generated_at']}`",
        "",
        "Scope: small EWAS Atlas tables + unpacked DataHub `sample_*.txt` packs.",
        "Large Atlas associations / probe annotations and matrix zips are out of scope.",
        "",
        "## Parse recipes",
        "",
        "| Source | Recipe |",
        "|--------|--------|",
        "| Atlas studies / cohorts | TSV, header row, tab-separated |",
        "| Atlas trait×trait | TSV matrix; col0 = trait; remaining = traits |",
        "| DataHub sample-info | R write.table (space + quotes); `read_r_style_table` |",
        "",
        "## Family → primary phenotype column",
        "",
        "| Family | Column |",
        "|--------|--------|",
    ]
    for family, col in sorted(report["family_value_column"].items()):
        lines.append(f"| `{family}` | `{col}` |")

    lines.extend(["", "## Atlas small tables", ""])
    for table in report["atlas_tables"]:
        lines.append(f"### `{table['table_id']}`")
        lines.append("")
        if not table.get("exists"):
            lines.append(f"Missing: `{table['path']}`")
            lines.append("")
            continue
        lines.append(f"- Path: `{table['path']}`")
        lines.append(f"- Bytes: {table.get('bytes')}")
        lines.append(f"- Shape: {table.get('n_rows')} × {table.get('n_cols')}")
        lines.append(f"- Parse: {table.get('parse_recipe')}")
        if table.get("n_malformed_rows_skipped"):
            lines.append(
                f"- Malformed rows skipped: {table.get('n_malformed_rows_skipped')}"
            )
        if table.get("join_keys"):
            lines.append(f"- Join keys: {', '.join(f'`{k}`' for k in table['join_keys'])}")
        if "columns" in table:
            lines.append("")
            lines.append("| Column | Kind | Non-null | N unique | Notes |")
            lines.append("|--------|------|---------:|---------:|-------|")
            for col in table["columns"]:
                notes = ""
                if col["kind"] == "numeric":
                    notes = f"min={col.get('min')}, max={col.get('max')}"
                elif col.get("top_values"):
                    notes = ", ".join(
                        f"{_md_escape(str(t['value']))} ({t['count']})"
                        for t in col["top_values"][:3]
                    )
                lines.append(
                    f"| `{col['name']}` | {col['kind']} | {col['non_null_rate']:.3f} | "
                    f"{col['n_unique']} | {notes} |"
                )
        if table.get("n_traits") is not None:
            lines.append(f"- Traits: {table['n_traits']}; square={table.get('matrix_is_square')}")
            lines.append(
                f"- Value sample range: [{table.get('value_sample_min')}, {table.get('value_sample_max')}]"
            )
        if table.get("example_rows"):
            lines.append("")
            lines.append("Example rows:")
            lines.append("")
            lines.append("```json")
            payload = json.dumps(table["example_rows"], indent=2)
            lines.append(payload if len(payload) < 2500 else payload[:2500] + "\n…")
            lines.append("```")
        lines.append("")

    lines.extend(["", "## DataHub sample-info packs", ""])
    for pack in report["sample_packs"]:
        lines.append(f"### Family `{pack['family']}`")
        lines.append("")
        if not pack.get("exists"):
            lines.append(f"Missing: `{pack['path']}`")
            lines.append("")
            continue
        lines.append(f"- Path: `{pack['path']}`")
        lines.append(f"- Bytes: {pack.get('bytes')}")
        lines.append(f"- Shape: {pack.get('n_rows')} × {pack.get('n_cols')}")
        lines.append(f"- Join keys: {', '.join(f'`{k}`' for k in pack.get('join_keys') or [])}")
        lines.append(f"- Primary phenotype: `{pack.get('primary_phenotype_column')}`")
        stats = pack.get("primary_phenotype_stats") or {}
        if stats:
            lines.append(
                f"- Primary non-null={stats.get('non_null_rate')}, "
                f"n_unique={stats.get('n_unique')}, kind={stats.get('kind')}"
            )
            if stats.get("top_values"):
                tops = ", ".join(f"{t['value']} ({t['count']})" for t in stats["top_values"][:5])
                lines.append(f"- Top values: {tops}")
            if stats.get("min") is not None:
                lines.append(
                    f"- Numeric range: [{stats.get('min')}, {stats.get('max')}] "
                    f"mean={stats.get('mean')}"
                )
        lines.append(
            f"- N sample_id / project_id: {pack.get('n_sample_id')} / {pack.get('n_project_id')}"
        )
        lines.append("")

    cross = report.get("cross_pack") or {}
    lines.extend(
        [
            "## Cross-pack column / ID overlap",
            "",
            f"- Packs profiled: {cross.get('n_packs')}",
            f"- Shared columns ({cross.get('n_shared_columns')}): "
            + ", ".join(f"`{c}`" for c in (cross.get("shared_columns") or [])),
            "",
            "Pairwise `sample_id` overlap:",
            "",
            "| A | B | Shared | N_A | N_B |",
            "|---|---|-------:|----:|----:|",
        ]
    )
    for row in cross.get("pairwise_sample_id_overlap") or []:
        lines.append(
            f"| `{row['a']}` | `{row['b']}` | {row['n_shared_sample_id']} | "
            f"{row['n_a']} | {row['n_b']} |"
        )

    lines.extend(
        [
            "",
            "Pairwise `project_id` overlap:",
            "",
            "| A | B | Shared | N_A | N_B |",
            "|---|---|-------:|----:|----:|",
        ]
    )
    for row in cross.get("pairwise_project_id_overlap") or []:
        lines.append(
            f"| `{row['a']}` | `{row['b']}` | {row['n_shared_project_id']} | "
            f"{row['n_a']} | {row['n_b']} |"
        )

    overlap = report.get("atlas_hub_id_overlap") or {}
    lines.extend(
        [
            "",
            "## Atlas study_ID vs Hub project_id",
            "",
            f"- Atlas study_ID count: {overlap.get('n_atlas_study_id')}",
            f"- Hub project_id count (union): {overlap.get('n_hub_project_id')}",
            f"- Hub GSE-like: {overlap.get('n_hub_gse_like')}",
            f"- Exact string equals: {overlap.get('n_exact_study_id_equals_project_id')}",
            f"- Note: {overlap.get('note')}",
            "",
            "## Related",
            "",
            "- Durable contracts: `docs/EWAS_METADATA.md`",
            "- Export path: `mbs.registry.sample_info`",
            "",
        ]
    )
    return "\n".join(lines)


def write_ewas_metadata_report(report: dict[str, Any], report_dir: Path) -> Path:
    report_dir = report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    (report_dir / "summary.md").write_text(
        render_ewas_metadata_markdown(report),
        encoding="utf-8",
    )
    return report_dir
