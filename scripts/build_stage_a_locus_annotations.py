#!/usr/bin/env python3
"""Build versioned locus annotation artifacts for Stage A N-light benchmark.

Outputs (under reports/inspection/stage0_7g_gene_only_probe/locus_annotations/):
  locus_annotations.parquet        per-locus annotation + presence flags
  locus_gene_annotations.parquet   per-(locus, gene) edge annotations
  annotation_manifest.json         provenance + channel vocabulary + hashes
  annotation_qc.json               QC counts, distributions, coverage stats
  annotation_qc.md                 human-readable QC report

Usage:
  uv run python scripts/build_stage_a_locus_annotations.py [--out-dir DIR]

All paths resolve under $MBS_DATA_ROOT / $MBS_ROOT per filesystem policy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mbs.paths import DataPaths
from mbs.matrix.store import matrix_store_paths, read_locus_index
from mbs.training.cascade_assign import build_cascade_assignment, gene_linked_col_index
from mbs.training.flat_region_features import (
    GENE_ROLES,
    CPG_CONTEXTS,
    REGULATORY_CHANNELS,
    PRESENCE_FLAGS,
    build_flat_region_gene_index,
    count_other_gene_edges,
)
from mbs.training.locus_gene import load_graph_tables

# ── defaults ─────────────────────────────────────────────────────────────────

GRAPH_ID = "graph-grch38-gencode38-cgi-tile-v2"
MATRIX_ID = "matrix-hub-age-tissue-sex-full-v1"
SPLIT_ID = "hub-ats-7e-3fold-v1"
MAX_LOCI = 65536


# ── helpers ───────────────────────────────────────────────────────────────────


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def _utc_now() -> str:
    import datetime

    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _percentile(arr: np.ndarray, p: float) -> float:
    return float(np.percentile(arr, p)) if arr.size else 0.0


def _write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


# ── assertions ────────────────────────────────────────────────────────────────


def _assert_vocabulary() -> None:
    """Channel vocabulary must match flat_region_features constants exactly."""
    expected_roles = list(GENE_ROLES)
    expected_ctx = list(CPG_CONTEXTS)
    expected_reg = list(REGULATORY_CHANNELS)
    # These are the canonical lists — if flat_region_features changes, this script
    # will still match because it imports the same tuples.
    assert len(expected_roles) == 6, f"GENE_ROLES len changed: {expected_roles}"
    assert len(expected_ctx) == 7, f"CPG_CONTEXTS len changed: {expected_ctx}"
    assert len(expected_reg) == 6, f"REGULATORY_CHANNELS len changed: {expected_reg}"


def _assert_no_duplicate_edges(edges: pd.DataFrame) -> None:
    dups = edges.duplicated(subset=["locus_id", "gene_id"]).sum()
    if dups:
        raise AssertionError(f"{dups} duplicate (locus_id, gene_id) pairs found")


def _assert_required_edge_fields(edges: pd.DataFrame) -> None:
    for col in ("locus_id", "gene_id", "evidence_type"):
        null_count = edges[col].isna().sum()
        if null_count:
            raise AssertionError(f"edge column '{col}' has {null_count} null values")


def _assert_genome_build(loci: pd.DataFrame) -> None:
    builds = loci["genome_build"].dropna().unique()
    bad = [b for b in builds if str(b) != "GRCh38"]
    if bad:
        raise AssertionError(f"mixed genome builds found: {bad}")


def _assert_one_hot_exclusive(edges: pd.DataFrame) -> None:
    """Each edge must have exactly one active gene role and one cpg context."""
    role_cols = [f"role_{r}" for r in GENE_ROLES]
    ctx_cols = [f"ctx_{c}" for c in CPG_CONTEXTS]
    for cols, name in [(role_cols, "gene_role"), (ctx_cols, "cpg_context")]:
        present = [c for c in cols if c in edges.columns]
        if present:
            row_sums = edges[present].sum(axis=1)
            bad = (row_sums != 1).sum()
            if bad:
                raise AssertionError(
                    f"{bad} edges have != 1 active value in {name} one-hot block"
                )


def _assert_no_other_gene(edges: pd.DataFrame) -> None:
    col = "role_other_gene"
    if col in edges.columns and edges[col].any():
        n = int(edges[col].sum())
        raise AssertionError(
            f"{n} edges have other_gene role on five-role graph (expected 0)"
        )


def _assert_regulatory_all_zero(edges: pd.DataFrame) -> None:
    reg_cols = [f"reg_{r}" for r in REGULATORY_CHANNELS]
    for c in reg_cols:
        if c in edges.columns and edges[c].any():
            raise AssertionError(
                f"Regulatory channel '{c}' is non-zero — source not yet on disk. "
                "Remove this assertion once the source is loaded."
            )


def _assert_panel_loci_match(
    panel_locus_ids: set[int], edge_locus_ids: set[int]
) -> None:
    diff = panel_locus_ids.symmetric_difference(edge_locus_ids)
    if diff:
        raise AssertionError(
            f"{len(diff)} locus_ids differ between gene_cols panel and edge table"
        )


# ── QC stats ──────────────────────────────────────────────────────────────────


def _compute_qc(
    locus_ann: pd.DataFrame,
    edge_ann: pd.DataFrame,
    gene_cols: np.ndarray,
    folds_path: Path | None,
) -> dict:
    n_loci = len(locus_ann)
    n_edges = len(edge_ann)
    n_genes = int(edge_ann["gene_id"].nunique())
    n_multi = int(
        (edge_ann.groupby("locus_id")["gene_id"].nunique() > 1).sum()
    )

    # edges per gene
    epg = edge_ann.groupby("gene_id").size().values
    cpg_per_gene = (
        edge_ann.groupby("gene_id")["locus_id"].nunique().values
    )

    def _dist(arr: np.ndarray) -> dict:
        if not arr.size:
            return {}
        return {
            "min": int(arr.min()),
            "median": float(np.median(arr)),
            "p90": _percentile(arr, 90),
            "p99": _percentile(arr, 99),
            "max": int(arr.max()),
            "mean": float(arr.mean()),
        }

    # role counts
    role_counts = {
        r: int((edge_ann["gene_role"] == r).sum()) for r in GENE_ROLES
    }
    role_pct = {r: round(100 * v / n_edges, 2) if n_edges else 0 for r, v in role_counts.items()}

    # context counts (per locus, not per edge, since context is locus-level)
    ctx_counts = {c: int((locus_ann["cpg_context"] == c).sum()) for c in CPG_CONTEXTS}
    ctx_pct = {c: round(100 * v / n_loci, 2) if n_loci else 0 for c, v in ctx_counts.items()}

    # presence rates
    ctx_present_rate = round(
        100 * float(locus_ann["cpg_context_present"].mean()), 2
    ) if n_loci else 0.0
    unknown_rate = round(100 * float((locus_ann["cpg_context"] == "unknown").mean()), 2) if n_loci else 0.0
    open_sea_rate = round(100 * float((locus_ann["cpg_context"] == "open_sea").mean()), 2) if n_loci else 0.0

    # regulatory (all zero expected)
    reg_counts = {r: 0 for r in REGULATORY_CHANNELS}
    reg_nonzero = {r: False for r in REGULATORY_CHANNELS}

    # CpG context × gene role cross-tabulation
    ctx_role_cross: dict = {}
    for ctx in CPG_CONTEXTS:
        ctx_edges = edge_ann[edge_ann["cpg_context"] == ctx]
        ctx_role_cross[ctx] = {
            r: int((ctx_edges["gene_role"] == r).sum()) for r in GENE_ROLES
        }

    # genes with zero observed CpGs: not computable without a beta matrix here;
    # report the min CpGs-per-gene from the panel instead
    genes_one_cpg = int((cpg_per_gene == 1).sum())

    # Coverage by CpGs-per-gene quantile
    q_bins = [0, 5, 10, 25, 50, 100, 200, 999999]
    q_labels = ["1-4", "5-9", "10-24", "25-49", "50-99", "100-199", "200+"]
    cpg_per_gene_coverage: dict = {}
    for i, label in enumerate(q_labels):
        lo, hi = q_bins[i], q_bins[i + 1]
        mask = (cpg_per_gene >= lo) & (cpg_per_gene < hi)
        cpg_per_gene_coverage[label] = {
            "n_genes": int(mask.sum()),
            "n_cpgs": int(cpg_per_gene[mask].sum()),
        }

    # Fold coverage
    fold_coverage: dict = {}
    if folds_path and folds_path.is_file():
        try:
            folds_df = pd.read_parquet(folds_path)
            panel_loci_set = set(edge_ann["locus_id"].unique())
            for fold_id in sorted(folds_df["fold_id"].unique()):
                fold_coverage[str(fold_id)] = {}
                for split in ("train", "val", "test"):
                    split_samples = folds_df.loc[
                        (folds_df["fold_id"] == fold_id) & (folds_df["split"] == split),
                        "sample_id",
                    ].tolist()
                    fold_coverage[str(fold_id)][split] = {
                        "n_samples": len(split_samples),
                        "panel_loci": len(panel_loci_set),
                    }
        except Exception as exc:
            fold_coverage["error"] = str(exc)

    return {
        "n_unique_cpgs": n_loci,
        "n_locus_gene_edges": n_edges,
        "n_unique_genes": n_genes,
        "n_multi_gene_cpgs": n_multi,
        "edges_per_gene": _dist(epg),
        "cpgs_per_gene": _dist(cpg_per_gene),
        "genes_with_one_cpg": genes_one_cpg,
        "gene_role_counts": role_counts,
        "gene_role_pct": role_pct,
        "cpg_context_counts": ctx_counts,
        "cpg_context_pct": ctx_pct,
        "cpg_context_present_rate_pct": ctx_present_rate,
        "cpg_context_unknown_rate_pct": unknown_rate,
        "cpg_context_open_sea_rate_pct": open_sea_rate,
        "regulatory_channel_nonzero_counts": reg_counts,
        "regulatory_channel_nonzero_flags": reg_nonzero,
        "regulatory_note": "All regulatory channels are zero; cCRE/DHS/ChromHMM sources not on disk.",
        "regulatory_co_occurrence_matrix": "all_zero",
        "cpg_context_x_gene_role_cross_tabulation": ctx_role_cross,
        "annotation_coverage_by_cpgs_per_gene": cpg_per_gene_coverage,
        "fold_coverage": fold_coverage,
    }


def _write_qc_md(qc: dict, path: Path) -> None:
    lines = [
        "# Stage A N-light Annotation QC Report\n",
        f"Panel: {qc['n_unique_cpgs']:,} unique CpGs | {qc['n_locus_gene_edges']:,} locus-gene edges "
        f"| {qc['n_unique_genes']:,} genes | {qc['n_multi_gene_cpgs']:,} multi-gene CpGs\n",
        "",
        "## Gene Role Distribution (edges)\n",
        "| Role | Count | % |",
        "|------|------:|--:|",
    ]
    for r in GENE_ROLES:
        lines.append(
            f"| {r} | {qc['gene_role_counts'].get(r, 0):,} | {qc['gene_role_pct'].get(r, 0)} |"
        )
    lines += [
        "",
        "## CpG Context Distribution (loci)\n",
        "| Context | Count | % |",
        "|---------|------:|--:|",
    ]
    for c in CPG_CONTEXTS:
        lines.append(
            f"| {c} | {qc['cpg_context_counts'].get(c, 0):,} | {qc['cpg_context_pct'].get(c, 0)} |"
        )
    lines += [
        "",
        f"- Context present (UCSC CGI covered): {qc['cpg_context_present_rate_pct']}%",
        f"- Unknown/open_sea rate: {qc['cpg_context_open_sea_rate_pct']}%",
        "",
        "## Regulatory Channels\n",
        "> All zero — cCRE/DHS/ChromHMM source files not yet on disk (Stage A non-goal).\n",
        "",
        "## Edges per Gene\n",
        f"min={qc['edges_per_gene'].get('min')} "
        f"median={qc['edges_per_gene'].get('median')} "
        f"p90={qc['edges_per_gene'].get('p90')} "
        f"p99={qc['edges_per_gene'].get('p99')} "
        f"max={qc['edges_per_gene'].get('max')}\n",
        "",
        "## CpGs per Gene\n",
        f"min={qc['cpgs_per_gene'].get('min')} "
        f"median={qc['cpgs_per_gene'].get('median')} "
        f"p90={qc['cpgs_per_gene'].get('p90')} "
        f"p99={qc['cpgs_per_gene'].get('p99')} "
        f"max={qc['cpgs_per_gene'].get('max')}\n",
        "",
        "## CpG Context × Gene Role Cross-Tabulation\n",
        "| context \\ role | " + " | ".join(GENE_ROLES) + " |",
        "|" + "---|" * (len(GENE_ROLES) + 1),
    ]
    for ctx in CPG_CONTEXTS:
        row_vals = [str(qc["cpg_context_x_gene_role_cross_tabulation"].get(ctx, {}).get(r, 0)) for r in GENE_ROLES]
        lines.append(f"| {ctx} | " + " | ".join(row_vals) + " |")

    lines += [
        "",
        "## Coverage by CpGs-per-Gene Quantile\n",
        "| Bin | Genes | CpGs |",
        "|-----|------:|-----:|",
    ]
    for bin_label, stats in qc["annotation_coverage_by_cpgs_per_gene"].items():
        lines.append(f"| {bin_label} | {stats['n_genes']:,} | {stats['n_cpgs']:,} |")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ── main ─────────────────────────────────────────────────────────────────────


def build_locus_annotations(
    paths: DataPaths,
    *,
    graph_id: str,
    matrix_id: str,
    out_dir: Path,
    max_loci: int = MAX_LOCI,
) -> None:
    print(f"[annotation] graph={graph_id} matrix={matrix_id}", flush=True)

    # 1. Load source tables
    graph_root = paths.data_root / "canonical" / "graphs" / graph_id
    loci_ann_path = paths.data_root / "canonical" / "annotations" / "loci.parquet"
    graph_manifest_path = graph_root / "graph_manifest.json"
    matrix_paths = matrix_store_paths(
        paths.data_root / "canonical" / "matrices" / matrix_id
    )

    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(graph_root)
    genes_path = graph_root / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()

    loci_ann = pd.read_parquet(loci_ann_path)
    graph_manifest = json.loads(graph_manifest_path.read_text(encoding="utf-8"))

    # 2. Vocabulary assertion
    _assert_vocabulary()

    # 3. Build assignment + gene index
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation="explicit_only",
    )
    gene_cols = gene_linked_col_index(assignment)

    # cpg_context dict keyed by locus_id (str)
    merged_ctx = locus_index.merge(
        loci_ann[["locus_id", "cpg_context"]], on="locus_id", how="left"
    )
    cpg_context_by_locus = dict(
        zip(
            merged_ctx.loc[merged_ctx["cpg_context"].notna(), "locus_id"].astype(str),
            merged_ctx.loc[merged_ctx["cpg_context"].notna(), "cpg_context"].astype(str),
        )
    )

    region_index = build_flat_region_gene_index(
        assignment,
        locus_index=locus_index,
        cpg_context_by_locus=cpg_context_by_locus,
        allow_other_gene=False,
    )

    # 4. Build locus_annotations table (panel loci only)
    gene_edge = assignment.region_to_gene[assignment.edge_region_index] >= 0
    type_ids = assignment.region_type_id[assignment.edge_region_index[gene_edge]]
    n_other = count_other_gene_edges(type_ids, assignment.region_types)

    # panel locus_ids: loci that appear in gene_cols
    panel_col_set = set(gene_cols.tolist())
    panel_locus_mask = locus_index["col_index"].isin(panel_col_set)
    panel_locus_ids_ser = locus_index.loc[panel_locus_mask, "locus_id"]
    panel_locus_id_set = set(panel_locus_ids_ser.tolist())

    locus_ann_panel = loci_ann[loci_ann["locus_id"].isin(panel_locus_id_set)].copy()

    # gene_role_present: any edge for this locus in the gene-linked panel
    edge_locus_ids = set(region_index.edge_col_index.tolist())
    # map col_index → locus_id for edges
    col_to_lid = dict(zip(locus_index["col_index"].tolist(), locus_index["locus_id"].tolist()))
    edge_locus_id_set = {col_to_lid[c] for c in region_index.edge_col_index if c in col_to_lid}

    locus_ann_panel["gene_role_present"] = locus_ann_panel["locus_id"].isin(edge_locus_id_set)
    locus_ann_panel["cpg_context_present"] = locus_ann_panel["cpg_context"] != "unknown"
    locus_ann_panel["ccre_coverage_present"] = False
    locus_ann_panel["dhs_coverage_present"] = False
    locus_ann_panel["chromhmm_coverage_present"] = False

    # Keep required columns only
    locus_out = locus_ann_panel[[
        "locus_id", "chromosome", "position", "genome_build", "canonical_key",
        "cpg_context", "cpg_context_present", "gene_role_present",
        "ccre_coverage_present", "dhs_coverage_present", "chromhmm_coverage_present",
    ]].reset_index(drop=True)

    # 5. Build locus_gene_annotations table
    lid_arr = region_index.edge_col_index  # col indices
    gene_idx_arr = region_index.edge_gene_index
    role_id_arr = region_index.edge_role_id
    ctx_id_arr = region_index.edge_context_id

    # Map col_index → locus_id
    lid_series = pd.Series(lid_arr).map(col_to_lid)
    gene_ids_list = [region_index.gene_ids[int(g)] for g in gene_idx_arr]
    gene_roles_list = [GENE_ROLES[int(r)] for r in role_id_arr]
    ctx_list = [CPG_CONTEXTS[int(c)] for c in ctx_id_arr]

    edge_ann = pd.DataFrame({
        "locus_id": lid_series.values,
        "gene_id": gene_ids_list,
        "gene_role": gene_roles_list,
        "gene_role_present": region_index.edge_role_present,
        "cpg_context": ctx_list,
        "cpg_context_present": region_index.edge_context_present,
        "regulatory_present": region_index.edge_regulatory_present,
        "evidence_type": "interval_overlap",
    })

    # Regulatory multi-hot columns (all zero)
    reg_arr = region_index.edge_regulatory_multi_hot
    for i, ch in enumerate(REGULATORY_CHANNELS):
        edge_ann[f"reg_{ch}"] = reg_arr[:, i].astype(np.float32)

    # One-hot role/context columns for assertion checks
    for j, r in enumerate(GENE_ROLES):
        edge_ann[f"role_{r}"] = (role_id_arr == j).astype(np.int8)
    for j, c in enumerate(CPG_CONTEXTS):
        edge_ann[f"ctx_{c}"] = (ctx_id_arr == j).astype(np.int8)

    # 6. Run all assertions
    print("[annotation] running assertions...", flush=True)
    _assert_genome_build(loci_ann)
    _assert_no_duplicate_edges(edge_ann)
    _assert_required_edge_fields(edge_ann)
    _assert_one_hot_exclusive(edge_ann)
    _assert_no_other_gene(edge_ann)
    _assert_regulatory_all_zero(edge_ann)
    _assert_panel_loci_match(panel_locus_id_set, set(lid_series.dropna().astype(int).tolist()))
    if n_other != 0:
        raise AssertionError(f"other_gene edges: {n_other} (expected 0 for five-role graph)")
    print(f"[annotation] all assertions passed. n_loci={len(locus_out)} n_edges={len(edge_ann)}", flush=True)

    # Drop one-hot helper columns from edge output
    drop_cols = [f"role_{r}" for r in GENE_ROLES] + [f"ctx_{c}" for c in CPG_CONTEXTS]
    edge_out = edge_ann.drop(columns=drop_cols).reset_index(drop=True)

    # 7. Write parquet artifacts
    out_dir.mkdir(parents=True, exist_ok=True)
    locus_parquet = out_dir / "locus_annotations.parquet"
    edge_parquet = out_dir / "locus_gene_annotations.parquet"
    locus_out.to_parquet(locus_parquet, index=False)
    edge_out.to_parquet(edge_parquet, index=False)
    print(f"[annotation] wrote {locus_parquet}", flush=True)
    print(f"[annotation] wrote {edge_parquet}", flush=True)

    locus_hash = _sha256(locus_parquet)
    edge_hash = _sha256(edge_parquet)
    loci_ann_hash = _sha256(loci_ann_path)

    # graph CGI source hash from graph_manifest source_files
    graph_cgi_sha = None
    for sf in graph_manifest.get("source_files", []):
        if "cpg" in str(sf.get("name", "")).lower() or "cgi" in str(sf.get("uri", "")).lower():
            graph_cgi_sha = sf.get("sha256")
            break

    # 8. Write annotation_manifest.json
    manifest: dict = {
        "schema_version": "1",
        "preprocessing_version": "stage-a-annotation-v1",
        "genome_assembly": "GRCh38",
        "coordinate_convention": "1-based_cytosine_from_infinium_0based_CpG_beg",
        "interval_conventions": "1-based_inclusive",
        "sources": {
            "loci": {
                "path": str(loci_ann_path),
                "sha256": loci_ann_hash,
                "version": "Milestone-3-locus-registry",
                "n_rows": len(loci_ann),
            },
            "graph": {
                "graph_id": graph_id,
                "path": str(graph_root),
                "graph_manifest_sha256": _sha256(graph_manifest_path),
                "content_hash": graph_manifest.get("content_hash"),
            },
            "cpg_island": {
                "embedded_in_loci_annotation": True,
                "sha256_via_graph_manifest": graph_cgi_sha,
                "version": "UCSC_hg38",
                "shore_bp": 2000,
                "shelf_bp": 4000,
            },
            "ccre": {
                "status": "unavailable",
                "note": "ENCODE cCRE source not on disk; reserved slots remain zero",
            },
            "dhs": {
                "status": "unavailable",
                "note": "DHS source not on disk; reserved slots remain zero",
            },
            "chromhmm": {
                "status": "unavailable",
                "note": "ChromHMM source not on disk; reserved slots remain zero",
            },
        },
        "channel_vocabulary": {
            "gene_roles": list(GENE_ROLES),
            "cpg_contexts": list(CPG_CONTEXTS),
            "regulatory_channels": list(REGULATORY_CHANNELS),
            "presence_flags": list(PRESENCE_FLAGS),
            "feature_dim": 1 + len(GENE_ROLES) + len(CPG_CONTEXTS) + len(REGULATORY_CHANNELS) + len(PRESENCE_FLAGS) + 1,
        },
        "null_unknown_conventions": {
            "cpg_context_unknown": (
                "UCSC CGI table has no island within shelf_bp of this locus; "
                "context is recorded as 'open_sea'. 'unknown' means CpG context "
                "could not be determined from the source."
            ),
        },
        "missingness_semantics": {
            "gene_role_present": (
                "True when the locus-gene edge maps to a named GENCODE region type. "
                "False means the edge exists but uses the other_gene fallback (absent on five-role graph)."
            ),
            "cpg_context_present": (
                "True when context is assigned from UCSC CGI coverage (not open_sea or unknown fallback). "
                "False means locus is in open sea or context source was unavailable."
            ),
            "ccre_coverage_present": "Always False; cCRE source not yet on disk.",
            "dhs_coverage_present": "Always False; DHS source not yet on disk.",
            "chromhmm_coverage_present": "Always False; ChromHMM source not yet on disk.",
            "regulatory_annotation_present": "Always False; all regulatory sources unavailable.",
        },
        "row_counts": {
            "n_loci": len(locus_out),
            "n_locus_gene_edges": len(edge_out),
            "n_genes": int(edge_out["gene_id"].nunique()),
        },
        "content_hashes": {
            "locus_annotations_sha256": locus_hash,
            "locus_gene_annotations_sha256": edge_hash,
        },
        "artifact_paths": {
            "locus_annotations": str(locus_parquet),
            "locus_gene_annotations": str(edge_parquet),
        },
        "creation_command": f"uv run python {Path(__file__).name}",
        "repository_commit": _git_commit(),
        "created_at": _utc_now(),
    }
    manifest_path = out_dir / "annotation_manifest.json"
    _write_json(manifest_path, manifest)
    print(f"[annotation] wrote {manifest_path}", flush=True)

    # 9. QC stats
    folds_path = (
        paths.data_root / "canonical" / "splits" / SPLIT_ID / "folds.parquet"
    )
    qc = _compute_qc(locus_out, edge_ann, gene_cols, folds_path if folds_path.is_file() else None)
    qc_path = out_dir / "annotation_qc.json"
    _write_json(qc_path, qc)
    print(f"[annotation] wrote {qc_path}", flush=True)

    qc_md_path = out_dir / "annotation_qc.md"
    _write_qc_md(qc, qc_md_path)
    print(f"[annotation] wrote {qc_md_path}", flush=True)

    print("[annotation] done.", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: reports/inspection/stage0_7g_gene_only_probe/locus_annotations/)",
    )
    parser.add_argument("--graph-id", default=GRAPH_ID)
    parser.add_argument("--matrix-id", default=MATRIX_ID)
    parser.add_argument("--max-loci", type=int, default=MAX_LOCI)
    args = parser.parse_args()

    paths = DataPaths.from_env()
    out_dir = args.out_dir or (
        paths.project_root
        / "reports"
        / "inspection"
        / "stage0_7g_gene_only_probe"
        / "locus_annotations"
    )
    build_locus_annotations(
        paths,
        graph_id=args.graph_id,
        matrix_id=args.matrix_id,
        out_dir=out_dir,
        max_loci=args.max_loci,
    )


if __name__ == "__main__":
    main()
