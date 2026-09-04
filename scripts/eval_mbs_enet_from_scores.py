#!/usr/bin/env python3
"""Post-hoc elastic-net readout on saved MBS / gene-linked RBS (no encoder retrain).

Supports:
- Cascade layout: ``artifacts/runs/<run-id>/fold_{i}/`` (zarr score blocks)
  - ``--which mbs`` → ``mbs_enet``
  - ``--which rbs`` → ``rbs_enet`` on ``all_gene_rbs.zarr``
  - ``--which both`` → both
- Flat light layout: ``artifacts/runs/<run-prefix>-f{i}/`` (``scores/mbs.npy`` +
  ``stage_a_probe_inputs.npz`` when present; MBS only)
"""

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
from mbs.training.flat_stage_a_eval import evaluate_flat_mbs_enet
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


def _phenotype_arrays(sample_ids: list[str], ph_by_id: dict) -> dict[str, np.ndarray]:
    return {
        "age": np.asarray([float(ph_by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64),
        "age_mask": np.asarray([bool(ph_by_id[s].age_mask) for s in sample_ids], dtype=bool),
        "tissue": np.asarray(
            [int(ph_by_id[s].class_index) if ph_by_id[s].tissue_mask else 0 for s in sample_ids],
            dtype=np.int64,
        ),
        "tissue_mask": np.asarray(
            [bool(ph_by_id[s].tissue_mask) for s in sample_ids], dtype=bool
        ),
        "sex": np.asarray(
            [
                int(ph_by_id[s].sex_class_index or 0) if ph_by_id[s].sex_mask else 0
                for s in sample_ids
            ],
            dtype=np.int64,
        ),
        "sex_mask": np.asarray([bool(ph_by_id[s].sex_mask) for s in sample_ids], dtype=bool),
        "study_ids": np.asarray(
            [str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object
        ),
    }


def patch_cascade_fold(
    *,
    fold_dir: Path,
    fold: dict,
    ph_by_id: dict,
    class_names: list[str],
    which: str = "mbs",
    force: bool = False,
) -> dict:
    metrics_path = fold_dir / "metrics.json"
    score_dir = fold_dir / "scores"
    sample_index = pd.read_parquet(score_dir / "sample_index.parquet")
    sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()
    train_idx, test_idx = _fold_index(fold, sample_ids)
    arrays = _phenotype_arrays(sample_ids, ph_by_id)
    blocks = load_cascade_score_blocks(score_dir)
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    evaluations = dict(metrics.get("evaluations") or {})
    written: dict[str, dict] = {}

    modes: list[tuple[str, dict]] = []
    if which in ("mbs", "both"):
        modes.append(("mbs_enet", blocks))
    if which in ("rbs", "both"):
        if "all_gene_rbs" not in blocks or int(blocks["all_gene_rbs"].shape[1]) == 0:
            raise FileNotFoundError(f"missing all_gene_rbs under {score_dir}")
        modes.append(("rbs_enet", {"mbs": blocks["all_gene_rbs"]}))

    for eval_key, score_blocks in modes:
        if not force and isinstance(evaluations.get(eval_key), dict):
            print(f"[{eval_key}] skip {fold_dir.name}: already present", flush=True)
            continue
        blob = _evaluate_mbs_enet(
            score_blocks,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=arrays["age"],
            age_mask=arrays["age_mask"],
            tissue=arrays["tissue"],
            tissue_mask=arrays["tissue_mask"],
            sex=arrays["sex"],
            sex_mask=arrays["sex_mask"],
            study_ids=arrays["study_ids"],
            class_names=class_names,
        )
        blob["evaluation"] = eval_key
        evaluations[eval_key] = blob
        written[eval_key] = blob
        tissue_f1 = (blob.get("metrics") or {}).get("tissue", {}).get("macro_f1")
        age_mae = (blob.get("metrics") or {}).get("age", {}).get("mae")
        sex_auroc = (blob.get("metrics") or {}).get("sex", {}).get("auroc")
        print(
            f"[{eval_key}] {fold_dir.name} tissue_f1={tissue_f1} age_mae={age_mae} "
            f"sex_auroc={sex_auroc} n_features={blob.get('n_score_features')}",
            flush=True,
        )

    metrics["evaluations"] = evaluations
    write_json(metrics_path, metrics)
    return written


def patch_flat_run(
    *,
    run_dir: Path,
    fold: dict,
    ph_by_id: dict,
    class_names: list[str],
) -> dict:
    """Post-hoc enet for FlatDeepSetRegion runs (``…-f{i}``)."""
    metrics_path = run_dir / "metrics.json"
    score_dir = run_dir / "scores"
    mbs_path = score_dir / "mbs.npy"
    if not mbs_path.is_file():
        raise FileNotFoundError(f"missing flat MBS scores: {mbs_path}")
    mbs = np.load(mbs_path)
    inputs_path = score_dir / "stage_a_probe_inputs.npz"
    if inputs_path.is_file():
        blob_np = np.load(inputs_path, allow_pickle=True)
        train_idx = np.asarray(blob_np["train_idx"], dtype=np.int64)
        test_idx = np.asarray(blob_np["test_idx"], dtype=np.int64)
        arrays = {
            "age": blob_np["age"],
            "age_mask": blob_np["age_mask"],
            "tissue": blob_np["tissue"],
            "tissue_mask": blob_np["tissue_mask"],
            "sex": blob_np["sex"],
            "sex_mask": blob_np["sex_mask"],
            "study_ids": blob_np["study_ids"],
        }
        meta_path = score_dir / "stage_a_probe_meta.json"
        if meta_path.is_file():
            class_names = list(
                json.loads(meta_path.read_text(encoding="utf-8")).get("class_names")
                or class_names
            )
    else:
        sample_ids = json.loads((score_dir / "sample_ids.json").read_text(encoding="utf-8"))
        train_ids = [str(s) for s in fold["train_sample_ids"]]
        val_ids = [str(s) for s in (fold.get("validation_sample_ids") or [])]
        test_ids = [str(s) for s in (fold.get("external_test_sample_ids") or [])]
        n_tr, n_va, n_te = len(train_ids), len(val_ids), len(test_ids)
        if len(sample_ids) != n_tr + n_va + n_te:
            train_idx, test_idx = _fold_index(fold, sample_ids)
        else:
            train_idx = np.arange(0, n_tr, dtype=np.int64)
            test_idx = np.arange(n_tr + n_va, n_tr + n_va + n_te, dtype=np.int64)
        arrays = _phenotype_arrays([str(s) for s in sample_ids], ph_by_id)
    out = evaluate_flat_mbs_enet(
        mbs_all=mbs,
        train_idx=train_idx,
        test_idx=test_idx,
        arrays=arrays,
        class_names=list(class_names),
    )
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    evaluations = dict(metrics.get("evaluations") or {})
    evaluations["mbs_enet"] = out
    metrics["evaluations"] = evaluations
    write_json(metrics_path, metrics)
    tissue_f1 = (out.get("metrics") or {}).get("tissue", {}).get("macro_f1")
    age_mae = (out.get("metrics") or {}).get("age", {}).get("mae")
    sex_auroc = (out.get("metrics") or {}).get("sex", {}).get("auroc")
    print(
        f"[mbs_enet] {run_dir.name} tissue_f1={tissue_f1} age_mae={age_mae} "
        f"sex_auroc={sex_auroc} n_features={out.get('n_score_features')}",
        flush=True,
    )
    return out


def _resolve_targets(
    *,
    paths: DataPaths,
    run_id: str | None,
    run_prefix: str | None,
    n_folds: int,
) -> list[tuple[str, Path]]:
    """Return ``(kind, path)`` where kind is ``cascade_fold`` or ``flat_run``."""
    if run_prefix:
        out: list[tuple[str, Path]] = []
        for i in range(n_folds):
            p = paths.artifact_root / "runs" / f"{run_prefix}-f{i}"
            out.append(("flat_run", p))
        return out
    if not run_id:
        raise ValueError("provide --run-id (cascade) or --run-prefix (flat)")
    run_root = paths.artifact_root / "runs" / run_id
    if (run_root / "fold_0").is_dir():
        return [("cascade_fold", run_root / f"fold_{i}") for i in range(n_folds)]
    # Single flat run id ending in -f0 style passed as run-id
    if run_id.endswith(("-f0", "-f1", "-f2")) or (run_root / "scores" / "mbs.npy").is_file():
        return [("flat_run", run_root)]
    raise FileNotFoundError(f"run not found or unsupported layout: {run_root}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=None, help="Cascade run dir or single flat run id")
    parser.add_argument(
        "--run-prefix",
        default=None,
        help="Flat arm prefix, e.g. stage0-7g-gene-probe-light-mean → …-f0/f1/f2",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--which",
        choices=("mbs", "rbs", "both"),
        default="mbs",
        help="Cascade only: evaluate mbs_enet, rbs_enet (all_gene_rbs), or both",
    )
    args = parser.parse_args()
    if not args.run_id and not args.run_prefix:
        args.run_id = DEFAULT_RUN
    paths = DataPaths.from_environment()
    cfg = load_experiment_config(args.config)
    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    folds = fold_pack["folds"]
    pheno_rel = Path(
        str(
            cfg.get("sample_phenotype_table")
            or (cfg.get("pilot") or {}).get("sample_phenotype_table")
            or "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}
    targets = _resolve_targets(
        paths=paths,
        run_id=args.run_id,
        run_prefix=args.run_prefix,
        n_folds=len(folds),
    )
    for i, (kind, target) in enumerate(targets):
        metrics_path = target / "metrics.json"
        if not metrics_path.is_file():
            print(f"[enet] skip {target.name}: metrics missing", flush=True)
            continue
        fold = folds[min(i, len(folds) - 1)]
        if kind == "cascade_fold":
            patch_cascade_fold(
                fold_dir=target,
                fold=fold,
                ph_by_id=ph_by_id,
                class_names=class_names,
                which=args.which,
                force=args.force,
            )
        else:
            if args.which == "rbs":
                print(f"[enet] skip {target.name}: --which rbs is cascade-only", flush=True)
                continue
            prior = json.loads(metrics_path.read_text(encoding="utf-8"))
            if not args.force and isinstance((prior.get("evaluations") or {}).get("mbs_enet"), dict):
                print(f"[mbs_enet] skip {target.name}: already present", flush=True)
                continue
            patch_flat_run(
                run_dir=target,
                fold=fold,
                ph_by_id=ph_by_id,
                class_names=class_names,
            )


if __name__ == "__main__":
    main()
