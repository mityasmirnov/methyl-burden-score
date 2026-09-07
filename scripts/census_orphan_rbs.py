#!/usr/bin/env python3
"""Census orphan RBS regions by CpG cardinality (CPU; no qualification change)."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.manifest import write_json
from mbs.matrix.store import matrix_store_paths, read_locus_index
from mbs.paths import DataPaths
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.locus_gene import load_graph_tables

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "reports/inspection/stage0_7g_prime_matched_probe/orphan_rbs_census"


def _region_cpg_counts(assignment: Any) -> Counter[int]:
    """Count edges per orphan region_id; return Counter of cardinality → n_regions."""
    orphan_idx = set(int(i) for i in assignment.orphan_region_indices.tolist())
    if not orphan_idx:
        return Counter()
    per_region: Counter[int] = Counter()
    for reg in assignment.edge_region_index.tolist():
        r = int(reg)
        if r in orphan_idx:
            per_region[r] += 1
    card = Counter(int(n) for n in per_region.values())
    return card


def census_one(
    *,
    data_root: Path,
    matrix_id: str,
    graph_id: str,
    max_loci: int,
) -> dict[str, Any]:
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    graph_root = data_root / "canonical" / "graphs" / graph_id
    lr_edges, regions = load_graph_tables(graph_root)
    genes_path = graph_root / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
    graph_hash = None
    graph_manifest = graph_root / "graph_manifest.json"
    if graph_manifest.is_file():
        graph_hash = json.loads(graph_manifest.read_text(encoding="utf-8")).get("content_hash")
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation="explicit_only",
    )
    card = _region_cpg_counts(assignment)
    n_orphan = int(assignment.n_orphan_rbs)
    n_singleton = int(card.get(1, 0))
    n_multi = int(sum(n for k, n in card.items() if k >= 2))
    return {
        "matrix_id": matrix_id,
        "graph_id": graph_id,
        "graph_content_hash": graph_hash,
        "max_loci": max_loci,
        "gene_allocation": "explicit_only",
        "n_study_loci": int(assignment.n_study_loci),
        "n_genes": int(assignment.n_genes),
        "n_regions": int(assignment.n_regions),
        "n_direct": int(assignment.n_direct),
        "n_orphan_rbs": n_orphan,
        "n_orphan_singleton_cpg": n_singleton,
        "n_orphan_multi_cpg": n_multi,
        "orphan_cpg_cardinality": {str(k): int(v) for k, v in sorted(card.items())},
        "note": (
            "Product rule: qualify multi-CpG orphans only; singletons → direct. "
            "This census does not change build_cascade_assignment."
        ),
    }


def write_report(out_stem: Path, payload: dict[str, Any]) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    write_json(Path(str(out_stem) + ".json"), payload)
    cohorts = payload.get("cohorts") or []
    lines = [
        "# Orphan RBS qualification census",
        "",
        "CPU census only — does **not** change assignment. Product Stage B / OOF "
        "should keep one column per **qualified multi-CpG** orphan `region_id`; "
        "singletons belong on the direct path.",
        "",
        f"Graph: `{payload.get('graph_id')}`",
        f"max_loci: `{payload.get('max_loci')}`",
        f"gene_allocation: `{payload.get('gene_allocation')}`",
        "",
        "| Matrix | n_orphan | singleton (1 CpG) | multi (≥2 CpG) | n_direct | n_genes |",
        "|--------|---------:|------------------:|---------------:|---------:|--------:|",
    ]
    for c in cohorts:
        lines.append(
            f"| `{c['matrix_id']}` | {c['n_orphan_rbs']} | {c['n_orphan_singleton_cpg']} | "
            f"{c['n_orphan_multi_cpg']} | {c['n_direct']} | {c['n_genes']} |"
        )
    lines.extend(
        [
            "",
            "## Cardinality detail",
            "",
        ]
    )
    for c in cohorts:
        lines.append(f"### `{c['matrix_id']}`")
        lines.append("")
        lines.append("| CpGs/region | n_orphan_regions |")
        lines.append("|------------:|-----------------:|")
        for k, v in (c.get("orphan_cpg_cardinality") or {}).items():
            lines.append(f"| {k} | {v} |")
        lines.append("")
    Path(str(out_stem) + ".md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-loci", type=int, default=65536)
    parser.add_argument(
        "--graph-id",
        default="graph-grch38-gencode38-cgi-tile-v2",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Output stem (writes .json and .md)",
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    matrix_ids = [
        "matrix-hub-age-tissue-sex-full-v1",
        "matrix-hub-nine-pack-virtual-v1",
    ]
    cohorts: list[dict[str, Any]] = []
    for mid in matrix_ids:
        root = paths.data_root / "canonical" / "matrices" / mid
        if not (root / "locus_index.parquet").is_file():
            print(f"[orphan-census] skip missing matrix {mid}", flush=True)
            continue
        print(f"[orphan-census] {mid}…", flush=True)
        cohorts.append(
            census_one(
                data_root=paths.data_root,
                matrix_id=mid,
                graph_id=args.graph_id,
                max_loci=args.max_loci,
            )
        )
    if not cohorts:
        raise SystemExit("no matrices found for orphan census")
    # Compare locus prefixes when both present.
    locus_prefix_equal = None
    if len(cohorts) == 2:
        a = read_locus_index(
            paths.data_root
            / "canonical"
            / "matrices"
            / cohorts[0]["matrix_id"]
            / "locus_index.parquet"
        )
        b = read_locus_index(
            paths.data_root
            / "canonical"
            / "matrices"
            / cohorts[1]["matrix_id"]
            / "locus_index.parquet"
        )
        n = int(args.max_loci)
        ids_a = a["locus_id"].astype(str).tolist()[:n]
        ids_b = b["locus_id"].astype(str).tolist()[:n]
        locus_prefix_equal = ids_a == ids_b
    payload = {
        "graph_id": args.graph_id,
        "max_loci": args.max_loci,
        "gene_allocation": "explicit_only",
        "locus_prefix_equal": locus_prefix_equal,
        "cohorts": cohorts,
    }
    out = args.out if args.out.is_absolute() else paths.project_root / args.out
    write_report(out, payload)
    print(f"[orphan-census] wrote {out}.md", flush=True)


if __name__ == "__main__":
    main()
