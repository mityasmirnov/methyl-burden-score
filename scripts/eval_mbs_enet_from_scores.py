#!/usr/bin/env python3
"""Post-hoc elastic-net readout on saved cascade MBS (no encoder retrain)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mbs.annotation.manifest import write_json
from mbs.paths import DataPaths
from mbs.training.cascade_loop import _evaluate_mbs_enet
from mbs.training.cascade_scores import load_cascade_score_blocks
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.loop import load_experiment_config
from mbs.training.phenotypes import load_multitask_phenotypes

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_gene_only_probe.yaml"
DEFAULT_RUN = "stage0-7g-gene-probe-P2-G-explicit"


def _fold_index(fold: dict, sample_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    train_set = {str(s) for s in fold["train_sample_ids"]}
    test_ids = fold.get("external_test_sample_ids") or fold.get("validation_sample_ids") or []
    test_set = {str(s) for s in test_ids}
    train_idx = np.asarray(
        [i for i, sid in enumerate(sample_ids) if sid in train_set], dtype=np.int64
    )
    test_idx = np.asarray(
        [i for i, sid in enumerate(sample_ids) if sid in test_set], dtype=np.int64
    )
    if train_idx.size == 0 or test_idx.size == 0:
        raise ValueError("empty train or test index for mbs_enet")
    return train_idx, test_idx


def patch_fold(
    *,
    fold_dir: Path,
    fold: dict,
    ph_by_id: dict,
    class_names: list[str],
) -> dict:
    metrics_path = fold_dir / "metrics.json"
    score_dir = fold_dir / "scores"
    sample_index = pd.read_parquet(score_dir / "sample_index.parquet")
    sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()
    train_idx, test_idx = _fold_index(fold, sample_ids)
    ages = np.asarray([float(ph_by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64)
    tissue = np.asarray(
        [int(ph_by_id[s].class_index) if ph_by_id[s].tissue_mask else 0 for s in sample_ids],
        dtype=np.int64,
    )
    sex = np.asarray(
        [
            int(ph_by_id[s].sex_class_index or 0) if ph_by_id[s].sex_mask else 0
            for s in sample_ids
        ],
        dtype=np.int64,
    )
    age_mask = np.asarray([bool(ph_by_id[s].age_mask) for s in sample_ids], dtype=bool)
    tissue_mask = np.asarray([bool(ph_by_id[s].tissue_mask) for s in sample_ids], dtype=bool)
    sex_mask = np.asarray([bool(ph_by_id[s].sex_mask) for s in sample_ids], dtype=bool)
    studies = np.asarray(
        [str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object
    )
    blocks = load_cascade_score_blocks(score_dir)
    blob = _evaluate_mbs_enet(
        blocks,
        train_idx=train_idx,
        test_idx=test_idx,
        ages=ages,
        age_mask=age_mask,
        tissue=tissue,
        tissue_mask=tissue_mask,
        sex=sex,
        sex_mask=sex_mask,
        study_ids=studies,
        class_names=class_names,
    )
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    evaluations = dict(metrics.get("evaluations") or {})
    evaluations["mbs_enet"] = blob
    metrics["evaluations"] = evaluations
    write_json(metrics_path, metrics)
    tissue_f1 = (blob.get("metrics") or {}).get("tissue", {}).get("macro_f1")
    age_mae = (blob.get("metrics") or {}).get("age", {}).get("mae")
    sex_auroc = (blob.get("metrics") or {}).get("sex", {}).get("auroc")
    print(
        f"[mbs_enet] {fold_dir.name} tissue_f1={tissue_f1} age_mae={age_mae} "
        f"sex_auroc={sex_auroc} n_features={blob.get('n_score_features')}",
        flush=True,
    )
    return blob


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    cfg = load_experiment_config(args.config)
    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    pheno_rel = Path(
        str(
            cfg.get("sample_phenotype_table")
            or "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}
    run_root = paths.artifact_root / "runs" / args.run_id
    if not run_root.is_dir():
        raise FileNotFoundError(f"run not found: {run_root}")
    for i, fold in enumerate(fold_pack["folds"]):
        fold_dir = run_root / f"fold_{i}"
        metrics_path = fold_dir / "metrics.json"
        if not metrics_path.is_file():
            print(f"[mbs_enet] skip {fold_dir.name}: metrics missing", flush=True)
            continue
        prior = json.loads(metrics_path.read_text(encoding="utf-8"))
        if not args.force and isinstance((prior.get("evaluations") or {}).get("mbs_enet"), dict):
            print(f"[mbs_enet] skip {fold_dir.name}: already present", flush=True)
            continue
        patch_fold(
            fold_dir=fold_dir,
            fold=fold,
            ph_by_id=ph_by_id,
            class_names=class_names,
        )


if __name__ == "__main__":
    main()
