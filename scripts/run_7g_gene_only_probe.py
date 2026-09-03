#!/usr/bin/env python3
"""Run 7G′ Stage A gene-only probe (P*-G + classical -G + orphan ablation)."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.manifest import write_json
from mbs.matrix.store import matrix_store_paths, read_locus_index
from mbs.paths import DataPaths
from mbs.training.cascade_assign import (
    GeneAllocationPolicy,
    build_cascade_assignment,
    gene_linked_col_index,
)
from mbs.training.classical_mvalue import run_classical_mvalue
from mbs.training.dev_cv import load_frozen_folds, samples_from_phenotype_table
from mbs.training.loop import load_experiment_config
from mbs.training.locus_gene import load_graph_tables

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_gene_only_probe.yaml"
DEFAULT_ONTOLOGY = "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml"


def _metric_from_fold(blob: dict[str, Any], metric_path: str) -> float | None:
    """Resolve dotted path like ``mbs_e2e.metrics.tissue.macro_f1`` on fold metrics."""
    parts = metric_path.split(".")
    eval_keys = ("mbs_e2e", "mbs_linear_probe", "mbs_enet", "fusion_full", "fusion_mbs_direct")
    if parts[0] in eval_keys:
        evaluations = blob.get("evaluations") or {}
        cur: Any = evaluations.get(parts[0])
        for key in parts[1:]:
            if not isinstance(cur, dict):
                return None
            cur = cur.get(key)
        if cur is None:
            return None
        return float(cur)
    cur: Any = blob
    for key in parts:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    if cur is None:
        return None
    return float(cur)


def _mean_metric(folds: list[dict[str, Any]], metric_path: str) -> float | None:
    vals = [_metric_from_fold(f, metric_path) for f in folds]
    nums = [v for v in vals if v is not None]
    if not nums:
        return None
    return float(np.mean(nums))


def load_cascade_arm_folds(paths: DataPaths, run_id: str) -> list[dict[str, Any]]:
    run_root = paths.artifact_root / "runs" / run_id
    folds: list[dict[str, Any]] = []
    for fold_dir in sorted(run_root.glob("fold_*")):
        metrics_path = fold_dir / "metrics.json"
        if metrics_path.is_file():
            folds.append(json.loads(metrics_path.read_text(encoding="utf-8")))
    return folds


def train_cascade_arm(
    *,
    paths: DataPaths,
    config_path: Path,
    run_id: str,
    device: str,
    staging_report_dir: Path,
) -> None:
    cmd = [
        "uv",
        "run",
        "mbs",
        "train",
        "cascade",
        "--config",
        str(config_path),
        "--run-id",
        run_id,
        "--device",
        device,
        "--report-dir",
        str(staging_report_dir),
        "--skip-if-done",
    ]
    print(f"[gene-probe train] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=paths.project_root, check=True)


def build_gene_cols(
    paths: DataPaths,
    *,
    matrix_id: str,
    graph_id: str,
    max_loci: int,
    gene_allocation: GeneAllocationPolicy = "explicit_only",
    max_nearest_gene_bp: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    matrix_paths = matrix_store_paths(paths.data_root / "canonical" / "matrices" / matrix_id)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(paths.data_root / "canonical" / "graphs" / graph_id)
    graph_root = paths.data_root / "canonical" / "graphs" / graph_id
    genes_path = graph_root / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation=gene_allocation,
        max_nearest_gene_bp=max_nearest_gene_bp,
    )
    gene_cols = gene_linked_col_index(assignment)
    if gene_cols.size == 0:
        raise RuntimeError("gene-linked CpG panel is empty")
    graph_manifest_path = graph_root / "graph_manifest.json"
    graph_hash = None
    if graph_manifest_path.is_file():
        graph_hash = json.loads(graph_manifest_path.read_text(encoding="utf-8")).get("content_hash")
    manifest = {
        "gene_allocation": gene_allocation,
        "max_nearest_gene_bp": max_nearest_gene_bp,
        "n_gene_cols": int(gene_cols.size),
        "graph_id": graph_id,
        "graph_content_hash": graph_hash,
        "matrix_id": matrix_id,
        "max_loci": int(max_loci),
        "gene_col_indices": gene_cols.astype(int).tolist(),
    }
    print(
        f"[gene-probe] gene_cols={gene_cols.size} allocation={gene_allocation} "
        f"max_loci={max_loci}",
        flush=True,
    )
    return gene_cols, manifest


def write_per_arm(report_dir: Path, arm_id: str, payload: dict[str, Any]) -> None:
    per_arm = report_dir / "per_arm"
    per_arm.mkdir(parents=True, exist_ok=True)
    write_json(per_arm / f"{arm_id}.json", payload)


def _slim_cascade_fold(blob: dict[str, Any]) -> dict[str, Any]:
    """Drop val_history and ROC curves from fold metrics for compact per_arm JSON."""
    out = dict(blob)
    ckpt = out.get("checkpoint_selection")
    if isinstance(ckpt, dict) and "val_history" in ckpt:
        out["checkpoint_selection"] = {k: v for k, v in ckpt.items() if k != "val_history"}
    evaluations = out.get("evaluations")
    if isinstance(evaluations, dict):
        slim_eval: dict[str, Any] = {}
        for key, ev in evaluations.items():
            if not isinstance(ev, dict):
                slim_eval[key] = ev
                continue
            slim_ev = dict(ev)
            metrics = slim_ev.get("metrics")
            if isinstance(metrics, dict):
                slim_m = dict(metrics)
                slim_m.pop("tissue_roc", None)
                slim_m.pop("tissue_by_study", None)
                sex = slim_m.get("sex")
                if isinstance(sex, dict):
                    slim_m["sex"] = {k: v for k, v in sex.items() if k not in {"fpr", "tpr"}}
                slim_ev["metrics"] = slim_m
            slim_eval[key] = slim_ev
        out["evaluations"] = slim_eval
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--skip-train", action="store_true")
    parser.add_argument("--arm", action="append", default=[])
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_rel = Path(str(cfg.get("report_dir", "reports/inspection/stage0_7g_gene_only_probe")))
    report_dir = report_rel if report_rel.is_absolute() else paths.project_root / report_rel
    report_dir.mkdir(parents=True, exist_ok=True)

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    max_loci = int(cfg.get("cv_budget", {}).get("max_loci", 65536))
    panel_cfg = cfg.get("gene_panel") or {}
    gene_allocation = str(panel_cfg.get("allocation", "explicit_only"))
    max_nearest_gene_bp = panel_cfg.get("max_nearest_gene_bp")
    if max_nearest_gene_bp is not None:
        max_nearest_gene_bp = int(max_nearest_gene_bp)
    arms = cfg.get("arms") or []
    pilot = cfg.get("pilot", {})
    matrix_id = str(pilot.get("matrix_id", "matrix-hub-age-tissue-sex-full-v1"))
    graph_id = str(pilot.get("graph_id", "graph-grch38-gencode38-cgi-tile-v2"))

    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)
    pheno_rel = Path(
        str(
            cfg.get("sample_phenotype_table")
            or "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    ont_path = paths.data_root / DEFAULT_ONTOLOGY
    _samples, phenotypes = samples_from_phenotype_table(pheno_path, ontology_path=ont_path)

    requested = set(args.arm)
    known = {str(a["id"]) for a in arms}
    unknown = sorted(requested - known)
    if unknown:
        raise ValueError(f"unknown --arm values: {unknown}")

    gene_cols, gene_panel_manifest = build_gene_cols(
        paths,
        matrix_id=matrix_id,
        graph_id=graph_id,
        max_loci=max_loci,
        gene_allocation=gene_allocation,  # type: ignore[arg-type]
        max_nearest_gene_bp=max_nearest_gene_bp,
    )
    write_json(report_dir / "gene_panel_manifest.json", gene_panel_manifest)

    completed: dict[str, list[dict[str, Any]]] = {}

    for arm in arms:
        arm_id = str(arm["id"])
        if requested and arm_id not in requested:
            continue
        kind = str(arm["kind"])
        if arm.get("optional"):
            gate_arm = str(arm.get("gate_arm", ""))
            ref_arm = str(arm.get("gate_reference_arm", ""))
            metric = str(arm.get("gate_metric", "mbs_e2e.metrics.tissue.macro_f1"))
            max_delta = float(arm.get("gate_max_delta", 0.03))
            gate_folds = completed.get(gate_arm) or load_cascade_arm_folds(
                paths, str((next(a for a in arms if str(a["id"]) == gate_arm))["run_id"])
            )
            ref_folds = completed.get(ref_arm) or load_cascade_arm_folds(
                paths, str((next(a for a in arms if str(a["id"]) == ref_arm))["run_id"])
            )
            gate_f1 = _mean_metric(gate_folds, metric)
            ref_f1 = _mean_metric(ref_folds, metric)
            if gate_f1 is None or ref_f1 is None:
                print(f"[gene-probe] skip optional {arm_id}: gate metrics missing", flush=True)
                continue
            if gate_f1 + max_delta < ref_f1:
                print(
                    f"[gene-probe] skip optional {arm_id}: {gate_arm} F1={gate_f1:.3f} "
                    f"not within {max_delta} of {ref_arm} F1={ref_f1:.3f}",
                    flush=True,
                )
                continue

        print(f"[gene-probe] arm {arm_id} kind={kind}", flush=True)
        if kind == "cascade_train":
            if not args.skip_train:
                rel_cfg = Path(str(arm["config"]))
                arm_cfg = rel_cfg if rel_cfg.is_absolute() else paths.project_root / rel_cfg
                staging = report_dir / f"_staging_{arm_id.replace('-', '_').lower()}_train"
                train_cascade_arm(
                    paths=paths,
                    config_path=arm_cfg,
                    run_id=str(arm["run_id"]),
                    device=args.device,
                    staging_report_dir=staging,
                )
            folds = load_cascade_arm_folds(paths, str(arm["run_id"]))
            completed[arm_id] = folds
            write_per_arm(
                report_dir,
                arm_id,
                {
                    "arm_id": arm_id,
                    "kind": kind,
                    "run_id": arm["run_id"],
                    "folds": [_slim_cascade_fold(f) for f in folds],
                    "mean_mbs_e2e_tissue_f1": _mean_metric(folds, "mbs_e2e.metrics.tissue.macro_f1"),
                    "mean_mbs_linear_probe_tissue_f1": _mean_metric(
                        folds, "mbs_linear_probe.metrics.tissue.macro_f1"
                    ),
                    "mean_mbs_enet_tissue_f1": _mean_metric(
                        folds, "mbs_enet.metrics.tissue.macro_f1"
                    ),
                    "mean_fusion_full_tissue_f1": _mean_metric(
                        folds, "fusion_full.metrics.tissue.macro_f1"
                    ),
                    "mean_fusion_mbs_direct_tissue_f1": _mean_metric(
                        folds, "fusion_mbs_direct.metrics.tissue.macro_f1"
                    ),
                },
            )
        elif kind == "classical_gene":
            payload = run_classical_mvalue(
                data_root=paths.data_root,
                fold_pack=fold_pack,
                phenotypes=phenotypes,
                max_loci=max_loci,
                matrix_id=matrix_id,
                gene_cols=gene_cols,
            )
            write_per_arm(report_dir, arm_id, {"arm_id": arm_id, "kind": kind, **payload})
        else:
            raise ValueError(f"unsupported arm kind: {kind!r}")

    summary = {
        "milestone": "7G-prime-stage-A",
        "report_dir": str(report_dir),
        "arms_run": sorted(completed.keys()) if completed else sorted(requested) if requested else [],
    }
    write_json(report_dir / "summary.json", summary)
    print(f"[gene-probe] wrote {report_dir / 'summary.json'}", flush=True)

    report_script = paths.project_root / "scripts" / "write_7g_gene_only_probe_report.py"
    if report_script.is_file():
        subprocess.run(
            ["uv", "run", "python", str(report_script), "--report-dir", str(report_dir)],
            cwd=paths.project_root,
            check=True,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
