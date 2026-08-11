#!/usr/bin/env python3
"""Write per-platform probe→region annotation coverage report."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from mbs.annotation.coverage import ROLE_ORDER, compute_annotation_coverage

DEFAULT_GRAPH_ID = "graph-grch38-gencode38-five-role-v1"


def _pct(n: int, d: int) -> str:
    if d == 0:
        return "n/a"
    return f"{100.0 * n / d:.2f}%"


def write_report(
    *,
    project_root: Path,
    data_root: Path,
    graph_id: str = DEFAULT_GRAPH_ID,
) -> Path:
    ann = data_root / "canonical" / "annotations"
    graph = data_root / "canonical" / "graphs" / graph_id
    probes = pd.read_parquet(ann / "probes.parquet")
    edges = pd.read_parquet(ann / "probe_locus_edges.parquet")
    loci = pd.read_parquet(ann / "loci.parquet")
    lr = pd.read_parquet(graph / "locus_region_edges.parquet")
    regions = pd.read_parquet(graph / "regions.parquet")
    payload = compute_annotation_coverage(
        probes=probes,
        probe_locus_edges=edges,
        loci=loci,
        locus_region_edges=lr,
        regions=regions,
        graph_id=graph_id,
    )
    graph_v1 = project_root / "reports" / "inspection" / "annotation_graph_v1" / "summary.json"
    if graph_v1.is_file():
        ref = json.loads(graph_v1.read_text(encoding="utf-8"))
        ref_plat = {r["platform_id"]: r for r in ref.get("platform_coverage", [])}
        mismatches: list[str] = []
        for row in payload["probe_level"]["platforms"]:
            r = ref_plat.get(row["platform_id"])
            if r is None:
                continue
            if int(r["n_probes"]) != row["n_probes"] or int(r["n_mapped"]) != row["n_mapped"]:
                mismatches.append(
                    f"{row['platform_id']}: coverage={row['n_probes']}/{row['n_mapped']} "
                    f"vs graph_v1={r['n_probes']}/{r['n_mapped']}"
                )
        payload["cross_check_annotation_graph_v1"] = {
            "path": str(graph_v1),
            "ok": not mismatches,
            "mismatches": mismatches,
            "n_unassigned_loci_ref": ref.get("n_unassigned_loci"),
            "n_unassigned_loci_here": payload["locus_level"]["n_unassigned_loci"],
        }

    out_dir = project_root / "reports" / "inspection" / "annotation_coverage_v1"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "summary.md").write_text(_markdown(payload), encoding="utf-8")
    return out_dir


def _markdown(payload: dict) -> str:
    loc = payload["locus_level"]
    pr = payload["probe_level"]
    assigned_pct = _pct(loc["n_assigned_loci"], loc["n_loci"])
    unassigned_pct = _pct(loc["n_unassigned_loci"], loc["n_loci"])
    unmapped_pct = _pct(pr["n_unmapped_probes"], pr["n_probes"])
    lines = [
        "# Annotation coverage (probe → locus → five-role region)",
        "",
        f"- graph_id: `{payload['graph_id']}`",
        f"- loci: {loc['n_loci']}",
        f"- assigned loci: {loc['n_assigned_loci']} ({assigned_pct})",
        f"- unassigned loci: {loc['n_unassigned_loci']} ({unassigned_pct})",
        f"- multi-gene loci: {loc['n_loci_multi_gene']}",
        f"- locus-region edges: {loc['n_locus_region_edges']}",
        f"- probes (all platforms): {pr['n_probes']}",
        f"- unmapped probes: {pr['n_unmapped_probes']} ({unmapped_pct})",
        "",
        "## Figures",
        "",
        "![Locus assigned vs unassigned](figures/locus_assigned_pie.png)",
        "",
        "![Per-array mapped probe assignment](figures/platform_assigned_vs_unassigned.png)",
        "",
        "![Loci by regulatory role](figures/loci_by_role.png)",
        "",
        "![CpG-island context](figures/cpg_island_context.png)",
        "",
        "## Definitions",
        "",
        f"- **Unmapped probe:** {payload['definitions']['unmapped_probe']}",
        f"- **Unassigned locus:** {payload['definitions']['unassigned_locus']}",
        f"- **Atlas:** {payload['definitions']['atlas_probe_annotations']}",
        "",
        "## Locus-level roles (unique loci with ≥1 edge of that type)",
        "",
    ]
    lines.extend(
        [
            (
                f"- `{role}`: {loc['loci_by_role'][role]} loci "
                f"({loc['locus_region_edges_by_role'][role]} edges)"
            )
            for role in ROLE_ORDER
        ]
    )
    lines.extend(
        [
            "",
            "## Per-platform probe coverage",
            "",
            "| Platform | Probes | Mapped | Unmapped | Mapped→region | Mapped unassigned |",
            "|----------|-------:|-------:|---------:|--------------:|------------------:|",
        ]
    )
    lines.extend(
        [
            (
                "| {platform_id} | {n_probes} | {n_mapped} | {n_unmapped} "
                "({pct_unmapped}%) | {n_mapped_assigned} ({pct_mapped_assigned}%) | "
                "{n_mapped_unassigned} ({pct_mapped_unassigned}%) |"
            ).format(**row)
            for row in pr["platforms"]
        ]
    )
    lines.extend(["", "## Per-platform probes by role (unique probes)", ""])
    for row in pr["platforms"]:
        lines.append(f"### `{row['platform_id']}`")
        lines.append("")
        lines.extend(f"- `{role}`: {row['probes_by_role'][role]}" for role in ROLE_ORDER)
        lines.append("")
    cc = payload.get("cross_check_annotation_graph_v1")
    if cc:
        lines.extend(
            [
                "## Cross-check vs `annotation_graph_v1`",
                "",
                f"- ok: `{cc['ok']}`",
                (
                    f"- unassigned loci here / ref: {cc['n_unassigned_loci_here']} / "
                    f"{cc['n_unassigned_loci_ref']}"
                ),
            ]
        )
        if cc["mismatches"]:
            lines.append("- mismatches:")
            lines.extend(f"  - {m}" for m in cc["mismatches"])
        lines.append("")
    lines.extend(
        [
            "## Regenerate",
            "",
            "```bash",
            "uv run python scripts/write_annotation_coverage_report.py",
            "uv sync --extra analysis  # once, for matplotlib",
            "uv run python scripts/write_pipeline_doc_figures.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    data_root = Path(os.environ.get("MBS_DATA_ROOT", str(project_root / "data")))
    out = write_report(project_root=project_root, data_root=data_root)
    print(json.dumps({"report_dir": str(out)}, indent=2))  # noqa: T201


if __name__ == "__main__":
    main()
