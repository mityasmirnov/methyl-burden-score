#!/usr/bin/env python3
"""Run Milestone 7G cascade tissue probe arms P0–P5 (idempotent)."""

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
from mbs.paths import DataPaths
from mbs.training.cascade_scores import fusion_feature_matrix, load_cascade_score_blocks
from mbs.training.dev_cv import load_frozen_folds, samples_from_phenotype_table
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.loop import load_experiment_config
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.transparent_hub import run_hub_transparent_arm

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_cascade_tissue_probe.yaml"
DEFAULT_ONTOLOGY = "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml"


def _phenotype_arrays_for_ids(
    phenotypes: list[Any],
    sample_ids: list[str],
) -> dict[str, np.ndarray]:
    ph_by_id = {p.sample_id: p for p in phenotypes}
    ages = np.asarray([float(ph_by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64)
    tissue = np.asarray(
        [int(ph_by_id[s].class_index) if ph_by_id[s].tissue_mask else 0 for s in sample_ids],
        dtype=np.int64,
    )
    sex = np.asarray(
        [int(ph_by_id[s].sex_class_index or 0) if ph_by_id[s].sex_mask else 0 for s in sample_ids],
        dtype=np.int64,
    )
    age_mask = np.asarray([bool(ph_by_id[s].age_mask) for s in sample_ids], dtype=bool)
    tissue_mask = np.asarray([bool(ph_by_id[s].tissue_mask) for s in sample_ids], dtype=bool)
    sex_mask = np.asarray([bool(ph_by_id[s].sex_mask) for s in sample_ids], dtype=bool)
    studies = np.asarray([str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object)
    return {
        "ages": ages,
        "tissue": tissue,
        "sex": sex,
        "age_mask": age_mask,
        "tissue_mask": tissue_mask,
        "sex_mask": sex_mask,
        "studies": studies,
    }


def _fold_row_indices(
    sample_ids: list[str],
    fold: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Train and external-test row indices into ``sample_ids`` (cascade score order)."""
    sid_to_row = {sid: i for i, sid in enumerate(sample_ids)}
    train_ids = [s for s in fold["train_sample_ids"] if s in sid_to_row]
    test_ids = [s for s in (fold.get("external_test_sample_ids") or []) if s in sid_to_row]
    if not test_ids:
        test_ids = [s for s in fold.get("validation_sample_ids", []) if s in sid_to_row]
    if not train_ids or not test_ids:
        raise ValueError("fold missing train or test samples in score index")
    train_idx = np.asarray([sid_to_row[s] for s in train_ids], dtype=np.int64)
    test_idx = np.asarray([sid_to_row[s] for s in test_ids], dtype=np.int64)
    return train_idx, test_idx


def verify_baseline_run(paths: DataPaths, run_id: str, n_folds: int = 3) -> None:
    run_root = paths.artifact_root / "runs" / run_id
    for fold_i in range(n_folds):
        manifest = run_root / f"fold_{fold_i}" / "scores" / "score_manifest.json"
        metrics = run_root / f"fold_{fold_i}" / "metrics.json"
        if not manifest.is_file() or not metrics.is_file():
            raise FileNotFoundError(
                f"P0 baseline incomplete: missing {manifest} or {metrics}. "
                "Run 7G cascade first (scripts/run_7g_methylation_eval_driver.sh)."
            )


def load_cascade_arm_folds(
    paths: DataPaths,
    run_id: str,
    *,
    expected_folds: int | None = None,
) -> list[dict[str, Any]]:
    run_root = paths.artifact_root / "runs" / run_id
    folds: list[dict[str, Any]] = []
    for fold_dir in sorted(run_root.glob("fold_*")):
        metrics_path = fold_dir / "metrics.json"
        if metrics_path.is_file():
            folds.append(json.loads(metrics_path.read_text(encoding="utf-8")))
    if expected_folds is not None and len(folds) != expected_folds:
        raise RuntimeError(
            f"run {run_id!r} has {len(folds)} completed folds; expected {expected_folds}"
        )
    return folds


def run_fusion_only_arm(
    *,
    paths: DataPaths,
    source_run_id: str,
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
    class_names: list[str],
    fusion: dict[str, Any],
    max_loci: int,
) -> list[dict[str, Any]]:
    """P1: refusion on saved score Zarrs without retraining."""
    run_root = paths.artifact_root / "runs" / source_run_id
    matrix_id = "matrix-hub-age-tissue-sex-full-v1"
    graph_id = "graph-grch38-gencode38-cgi-tile-v2"
    del matrix_id, graph_id, max_loci
    out_folds: list[dict[str, Any]] = []
    for fold_i, fold in enumerate(fold_pack["folds"]):
        score_dir = run_root / f"fold_{fold_i}" / "scores"
        sample_index = pd.read_parquet(score_dir / "sample_index.parquet")
        sample_ids = sample_index["sample_id"].astype(str).tolist()
        blocks = load_cascade_score_blocks(score_dir)
        x = fusion_feature_matrix(blocks)
        train_idx, test_idx = _fold_row_indices(sample_ids, fold)
        ph = _phenotype_arrays_for_ids(phenotypes, sample_ids)
        fused = evaluate_late_fusion(
            scores_train=x[train_idx],
            scores_test=x[test_idx],
            age_train=ph["ages"][train_idx],
            age_mask_train=ph["age_mask"][train_idx],
            tissue_train=ph["tissue"][train_idx],
            tissue_mask_train=ph["tissue_mask"][train_idx],
            sex_train=ph["sex"][train_idx],
            sex_mask_train=ph["sex_mask"][train_idx],
            age_test=ph["ages"][test_idx],
            age_mask_test=ph["age_mask"][test_idx],
            tissue_test=ph["tissue"][test_idx],
            tissue_mask_test=ph["tissue_mask"][test_idx],
            sex_test=ph["sex"][test_idx],
            sex_mask_test=ph["sex_mask"][test_idx],
            study_ids_test=ph["studies"][test_idx],
            tissue_class_names=list(class_names) if class_names else None,
            fusion=fusion,
        )
        fused["fold_id"] = fold.get("fold_id", fold_i)
        fused["score_dir"] = str(score_dir)
        out_folds.append(fused)
    return out_folds


def run_transparent_region_arm(
    *,
    paths: DataPaths,
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
    max_loci: int,
    matrix_id: str,
    graph_id: str,
) -> list[dict[str, Any]]:
    """P3: tissue head on region-mean features (eval only)."""
    out_folds: list[dict[str, Any]] = []
    for fold in fold_pack["folds"]:
        blob = run_hub_transparent_arm(
            data_root=paths.data_root,
            fold=fold,
            phenotypes=phenotypes,
            arm="T-mean-region",
            max_loci=max_loci,
            matrix_id=matrix_id,
            graph_id=graph_id,
        )
        out_folds.append(blob)
    return out_folds


def train_cascade_arm(
    *,
    paths: DataPaths,
    config_path: Path,
    run_id: str,
    device: str,
    staging_report_dir: Path,
) -> None:
    """Train one Phase-2 cascade arm with its resolved config."""
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
    print(f"[probe train] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=paths.project_root, check=True)


def write_per_arm(
    report_dir: Path,
    arm_id: str,
    folds: list[dict[str, Any]],
    *,
    extra: dict[str, Any] | None = None,
) -> None:
    per_arm = report_dir / "per_arm"
    per_arm.mkdir(parents=True, exist_ok=True)

    def _slim_fold(blob: dict[str, Any]) -> dict[str, Any]:
        """Drop bulky ROC / per-study arrays; keep scalar metrics for audit."""
        out = dict(blob)
        metrics = out.get("metrics")
        if isinstance(metrics, dict):
            slim_m = dict(metrics)
            slim_m.pop("tissue_roc", None)
            slim_m.pop("tissue_by_study", None)
            sex = slim_m.get("sex")
            if isinstance(sex, dict):
                slim_sex = {k: v for k, v in sex.items() if k not in {"fpr", "tpr"}}
                slim_m["sex"] = slim_sex
            out["metrics"] = slim_m
        ckpt = out.get("checkpoint_selection")
        if isinstance(ckpt, dict) and "val_history" in ckpt:
            slim_ckpt = {k: v for k, v in ckpt.items() if k != "val_history"}
            out["checkpoint_selection"] = slim_ckpt
        return out

    payload: dict[str, Any] = {
        "arm_id": arm_id,
        "folds": [_slim_fold(f) for f in folds],
    }
    if extra:
        payload.update(extra)
    write_json(per_arm / f"{arm_id}.json", payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--device", default="cuda")
    parser.add_argument(
        "--skip-train",
        "--skip-p2",
        dest="skip_train",
        action="store_true",
        help="Skip all cascade_train arms (report/refusion smoke or no GPU)",
    )
    parser.add_argument(
        "--arm",
        action="append",
        default=[],
        help="Run only the named arm; repeat for multiple arms",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_rel = Path(
        str(cfg.get("report_dir", "reports/inspection/stage0_7g_cascade_tissue_probe"))
    )
    report_dir = report_rel if report_rel.is_absolute() else paths.project_root / report_rel
    report_dir.mkdir(parents=True, exist_ok=True)

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    max_loci = int(cfg.get("cv_budget", {}).get("max_loci", 65536))
    baseline_run_id = str(cfg.get("baseline_run_id", "stage0-7g-cascade-v1"))
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
    _, class_names = load_multitask_phenotypes(pheno_path)

    requested_arms = set(args.arm)
    known_arms = {str(arm["id"]) for arm in arms}
    unknown_arms = sorted(requested_arms - known_arms)
    if unknown_arms:
        raise ValueError(f"unknown --arm values: {unknown_arms}")

    for arm in arms:
        arm_id = str(arm["id"])
        if requested_arms and arm_id not in requested_arms:
            continue
        kind = str(arm["kind"])
        print(f"[probe] arm {arm_id} kind={kind}", flush=True)
        if kind == "cascade_replay":
            run_id = str(arm.get("run_id", baseline_run_id))
            verify_baseline_run(paths, run_id, n_folds=len(fold_pack["folds"]))
            folds = load_cascade_arm_folds(
                paths, run_id, expected_folds=len(fold_pack["folds"])
            )
            write_per_arm(report_dir, arm_id, folds, extra={"run_id": run_id, "kind": kind})
        elif kind == "fusion_only":
            source_run_id = str(arm.get("source_run_id", baseline_run_id))
            verify_baseline_run(paths, source_run_id, n_folds=len(fold_pack["folds"]))
            fusion = dict(arm.get("fusion") or {})
            folds = run_fusion_only_arm(
                paths=paths,
                source_run_id=source_run_id,
                fold_pack=fold_pack,
                phenotypes=phenotypes,
                class_names=list(class_names) if class_names else [],
                fusion=fusion,
                max_loci=max_loci,
            )
            write_per_arm(
                report_dir,
                arm_id,
                folds,
                extra={"source_run_id": source_run_id, "fusion": fusion, "kind": kind},
            )
        elif kind == "cascade_train":
            if args.skip_train:
                print(f"[probe] skipping {arm_id} (--skip-train)", flush=True)
                continue
            run_id = str(arm.get("run_id", "stage0-7g-tissue-probe-P2"))
            arm_config_rel = Path(
                str(arm.get("config", "configs/experiment/stage0_7g_cascade_tissue_probe_p2.yaml"))
            )
            arm_config_path = (
                arm_config_rel
                if arm_config_rel.is_absolute()
                else paths.project_root / arm_config_rel
            )
            staging = report_dir / f"_staging_{arm_id}_train"
            staging.mkdir(parents=True, exist_ok=True)
            train_cascade_arm(
                paths=paths,
                config_path=arm_config_path,
                run_id=run_id,
                device=args.device,
                staging_report_dir=staging,
            )
            folds = load_cascade_arm_folds(
                paths, run_id, expected_folds=len(fold_pack["folds"])
            )
            write_per_arm(report_dir, arm_id, folds, extra={"run_id": run_id, "kind": kind})
        elif kind == "transparent_region":
            folds = run_transparent_region_arm(
                paths=paths,
                fold_pack=fold_pack,
                phenotypes=phenotypes,
                max_loci=max_loci,
                matrix_id=matrix_id,
                graph_id=graph_id,
            )
            write_per_arm(
                report_dir,
                arm_id,
                folds,
                extra={"transparent_arm": arm.get("arm", "T-mean-region"), "kind": kind},
            )
        else:
            raise ValueError(f"unknown probe arm kind: {kind}")

    report_script = ROOT / "scripts" / "write_7g_cascade_tissue_probe_report.py"
    cmd = [
        sys.executable,
        str(report_script),
        "--config",
        str(config_path),
        "--report-dir",
        str(report_dir),
    ]
    print(f"[probe] report {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=paths.project_root, check=True)
    print(f"[probe] done report={report_dir}", flush=True)


if __name__ == "__main__":
    main()
