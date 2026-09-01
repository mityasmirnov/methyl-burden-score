#!/usr/bin/env python3
"""Run 7G′ Stage B: fold-safe enetS, N-cascade-S, N-light-type, fusion ablations."""

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
from mbs.matrix.store import matrix_store_paths, open_betas_zarr, read_locus_index, read_sample_index
from mbs.paths import DataPaths
from mbs.training.cascade_assign import assignment_col_subset, build_cascade_assignment
from mbs.training.cascade_loop import train_cascade_on_arrays
from mbs.training.classical_mvalue import run_classical_mvalue_enetS
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.flat_region_loop import train_flat_region_on_arrays
from mbs.training.fold_safe_panel import select_fold_panel
from mbs.training.loop import load_experiment_config
from mbs.training.locus_gene import load_graph_tables
from mbs.training.phenotypes import load_multitask_phenotypes

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_prime_stage_b.yaml"


def _phenotype_arrays(phenotypes: list[Any], sample_ids: list[str]) -> dict[str, np.ndarray]:
    by_id = {p.sample_id: p for p in phenotypes}
    tissue = np.asarray([int(by_id[s].class_index) for s in sample_ids], dtype=np.int64)
    tissue_mask = np.asarray([bool(by_id[s].tissue_mask) for s in sample_ids], dtype=bool)
    age = np.asarray([float(by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64)
    age_mask = np.asarray([bool(by_id[s].age_mask) for s in sample_ids], dtype=bool)
    sex = np.asarray([int(by_id[s].sex_class_index or 0) for s in sample_ids], dtype=np.int64)
    sex_mask = np.asarray([bool(by_id[s].sex_mask) for s in sample_ids], dtype=bool)
    return {
        "tissue": tissue,
        "tissue_mask": tissue_mask,
        "age": age,
        "age_mask": age_mask,
        "sex": sex,
        "sex_mask": sex_mask,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_rel = Path(str(cfg.get("report_dir", "reports/inspection/stage0_7g_prime_matched_probe")))
    report_dir = report_rel if report_rel.is_absolute() else paths.project_root / report_rel
    report_dir.mkdir(parents=True, exist_ok=True)

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    max_loci = int(cfg.get("cv_budget", {}).get("max_loci", 65536))
    max_epochs = int(cfg.get("cv_budget", {}).get("max_epochs", 15))
    pilot = cfg.get("pilot", {})
    matrix_id = str(pilot.get("matrix_id", "matrix-hub-age-tissue-sex-full-v1"))
    graph_id = str(pilot.get("graph_id", "graph-grch38-gencode38-cgi-tile-v2"))
    panel_cfg = cfg.get("panel") or {}
    max_seeds = int(panel_cfg.get("max_seeds", 10_000))

    pheno_rel = Path(str(cfg.get("sample_phenotype_table")))
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}

    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)

    matrix_paths = matrix_store_paths(paths.data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(paths.data_root / "canonical" / "graphs" / graph_id)
    genes_path = paths.data_root / "canonical" / "graphs" / graph_id / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
    )
    n_cols = assignment.n_study_loci
    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }
    betas_all = np.asarray(open_betas_zarr(matrix_paths.betas_path)[:, :n_cols], dtype=np.float32)
    locus_ids = locus_index["locus_id"].astype(str).tolist()[:n_cols]

    lock_path = paths.project_root / str(
        cfg.get("lock_from_stage_a", "reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json")
    )
    lock: dict[str, Any] = {}
    if lock_path.is_file():
        lock = json.loads(lock_path.read_text(encoding="utf-8"))

    results: dict[str, Any] = {"milestone": "7G-prime-stage-B", "lock_from_stage_a": lock, "folds": []}

    enetS_payload = run_classical_mvalue_enetS(
        data_root=paths.data_root,
        fold_pack=fold_pack,
        phenotypes=phenotypes,
        assignment=assignment,
        max_loci=max_loci,
        matrix_id=matrix_id,
        max_seeds=max_seeds,
    )
    write_json(report_dir / "per_arm" / "C-mvalue-enetS.json", enetS_payload)

    for fold_idx, fold in enumerate(fold_pack["folds"]):
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id and s in ph_by_id]
        external = fold.get("external_test_sample_ids") or []
        test_ids = [s for s in external if s in row_by_id and s in ph_by_id]
        if not test_ids:
            test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id and s in ph_by_id]
        sample_ids = train_ids + test_ids
        rows = np.asarray([row_by_id[s] for s in sample_ids], dtype=np.int64)
        betas = betas_all[rows]
        train_idx = np.arange(0, len(train_ids), dtype=np.int64)
        test_idx = np.arange(len(train_ids), len(sample_ids), dtype=np.int64)
        ph = _phenotype_arrays(phenotypes, sample_ids)
        x_tr = betas[train_idx]
        panel_info = select_fold_panel(
            x_train=x_tr[ph["tissue_mask"][train_idx]],
            y_train=ph["tissue"][train_idx][ph["tissue_mask"][train_idx]],
            assignment=assignment,
            max_seeds=max_seeds,
        )
        panel = np.asarray(panel_info["panel_cols"], dtype=np.int64)
        studies = np.asarray([str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object)

        fold_out: dict[str, Any] = {"fold": fold_idx, "panel": panel_info, "arms": {}}

        cascade_metrics = train_cascade_on_arrays(
            assignment=assignment_col_subset(assignment, panel),
            betas=betas,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ph["age"],
            tissue=ph["tissue"],
            sex=ph["sex"],
            study_ids=studies,
            sample_ids=sample_ids,
            class_names=class_names or ["A", "B"],
            out_dir=report_dir / f"_staging_N_cascade_S_fold_{fold_idx}",
            max_epochs=max_epochs,
            seed=42 + fold_idx,
            device_str=args.device,
            gene_linked_only=False,
            primary_evaluation="mbs_e2e",
            locus_ids=locus_ids,
        )
        fold_out["arms"]["N-cascade-S"] = cascade_metrics

        flat_metrics = train_flat_region_on_arrays(
            assignment=assignment,
            betas=betas,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ph["age"],
            tissue=ph["tissue"],
            sex=ph["sex"],
            study_ids=studies,
            sample_ids=sample_ids,
            class_names=class_names or ["A", "B"],
            out_dir=report_dir / f"_staging_N_light_type_fold_{fold_idx}",
            max_epochs=max_epochs,
            seed=42 + fold_idx,
            device_str=args.device,
            age_mask=ph["age_mask"],
            tissue_mask=ph["tissue_mask"],
            sex_mask=ph["sex_mask"],
            panel_cols=panel,
        )
        fold_out["arms"]["N-light-type"] = flat_metrics

        for arm_id, extra_modes in (
            ("N-full", ()),
            ("N-mbs-direct-only", ("mbs_direct",)),
        ):
            full_metrics = train_cascade_on_arrays(
                assignment=assignment_col_subset(assignment, panel),
                betas=betas,
                train_idx=train_idx,
                test_idx=test_idx,
                ages=ph["age"],
                tissue=ph["tissue"],
                sex=ph["sex"],
                study_ids=studies,
                sample_ids=sample_ids,
                class_names=class_names or ["A", "B"],
                out_dir=report_dir / f"_staging_{arm_id.replace('-', '_')}_fold_{fold_idx}",
                max_epochs=max_epochs,
                seed=42 + fold_idx,
                device_str=args.device,
                gene_linked_only=False,
                primary_evaluation="late_fusion",
                extra_fusion_modes=extra_modes,
                locus_ids=locus_ids,
            )
            fold_out["arms"][arm_id] = full_metrics

        results["folds"].append(fold_out)

    write_json(report_dir / "summary.json", results)
    report_script = paths.project_root / "scripts" / "write_7g_prime_stage_b_report.py"
    if report_script.is_file():
        subprocess.run(
            ["uv", "run", "python", str(report_script), "--report-dir", str(report_dir)],
            cwd=paths.project_root,
            check=True,
        )
    print(f"[stage-b] wrote {report_dir / 'summary.json'}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
