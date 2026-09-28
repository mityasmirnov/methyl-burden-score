#!/usr/bin/env python3
"""Post-hoc nested elastic-net for G1 gene-holdout (train-gene vs heldout-gene MBS).

Supports:
- Flat N-light: ``artifacts/runs/<run-id>/scores/{mbs.npy,mbs_heldout.npy}``
- Cascade fold: ``artifacts/runs/<run-id>/fold_{i}/scores/{mbs.zarr,mbs_heldout.npy}``

Fits the same nested enet readout on each column set and writes
``scores/gene_holdout_eval.json`` (also patches ``metrics.json``).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from mbs.annotation.manifest import write_json
from mbs.paths import DataPaths
from mbs.training.cascade_scores import load_cascade_score_blocks
from mbs.training.transparent_baselines import run_nested_elasticnet_multitask

ROOT = Path(__file__).resolve().parents[1]


def _summary(out: dict[str, Any]) -> dict[str, Any]:
    metrics = out.get("metrics") or {}
    return {
        "tissue_macro_f1": (metrics.get("tissue") or {}).get("macro_f1"),
        "age_mae": (metrics.get("age") or {}).get("mae"),
        "sex_auroc": (metrics.get("sex") or {}).get("auroc"),
        "n_score_features": out.get("n_score_features") or out.get("n_features"),
    }


def _nested_on_mbs(
    mbs: np.ndarray,
    *,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    arrays: dict[str, np.ndarray],
    class_names: list[str],
    evaluation: str,
) -> dict[str, Any]:
    out = run_nested_elasticnet_multitask(
        x_train=mbs[train_idx],
        x_test=mbs[test_idx],
        age_train=arrays["age"][train_idx],
        age_mask_train=arrays["age_mask"][train_idx],
        tissue_train=arrays["tissue"][train_idx],
        tissue_mask_train=arrays["tissue_mask"][train_idx],
        sex_train=arrays["sex"][train_idx],
        sex_mask_train=arrays["sex_mask"][train_idx],
        age_test=arrays["age"][test_idx],
        age_mask_test=arrays["age_mask"][test_idx],
        tissue_test=arrays["tissue"][test_idx],
        tissue_mask_test=arrays["tissue_mask"][test_idx],
        sex_test=arrays["sex"][test_idx],
        sex_mask_test=arrays["sex_mask"][test_idx],
        study_ids_train=arrays["study_ids"][train_idx],
        study_ids_test=arrays["study_ids"][test_idx],
        tissue_class_names=list(class_names) if class_names else None,
    )
    out["evaluation"] = evaluation
    out["eval_split"] = "test"
    out["n_eval_samples"] = int(test_idx.size)
    out["n_score_features"] = int(out.get("n_features", mbs.shape[1]))
    return out


def _resolve_score_dir(run_dir: Path, fold: int | None) -> Path:
    flat = run_dir / "scores"
    if (flat / "mbs_heldout.npy").is_file():
        return flat
    fold_i = 0 if fold is None else int(fold)
    cascade = run_dir / f"fold_{fold_i}" / "scores"
    if (cascade / "mbs_heldout.npy").is_file():
        return cascade
    raise FileNotFoundError(
        f"no gene-holdout scores under {run_dir} (tried scores/ and fold_{fold_i}/scores/)"
    )


def _load_train_mbs(score_dir: Path) -> np.ndarray:
    npy = score_dir / "mbs.npy"
    if npy.is_file():
        return np.load(npy)
    zarr_path = score_dir / "mbs.zarr"
    if zarr_path.exists():
        return load_cascade_score_blocks(score_dir)["mbs"]
    raise FileNotFoundError(f"missing train MBS under {score_dir}")


def _load_split_arrays(
    score_dir: Path,
    run_dir: Path,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], list[str]]:
    """Return train_idx, test_idx, phenotype arrays, class_names."""
    inputs_path = score_dir / "stage_a_probe_inputs.npz"
    meta_path = score_dir / "stage_a_probe_meta.json"
    if inputs_path.is_file():
        blob = np.load(inputs_path, allow_pickle=True)
        arrays = {
            "age": blob["age"],
            "age_mask": blob["age_mask"],
            "tissue": blob["tissue"],
            "tissue_mask": blob["tissue_mask"],
            "sex": blob["sex"],
            "sex_mask": blob["sex_mask"],
            "study_ids": blob["study_ids"],
        }
        class_names = ["tissue"]
        if meta_path.is_file():
            class_names = list(
                json.loads(meta_path.read_text(encoding="utf-8")).get("class_names")
                or class_names
            )
        return (
            np.asarray(blob["train_idx"], dtype=np.int64),
            np.asarray(blob["test_idx"], dtype=np.int64),
            arrays,
            class_names,
        )

    pheno_path = score_dir / "gene_holdout_pheno.npz"
    if not pheno_path.is_file():
        raise FileNotFoundError(
            f"need stage_a_probe_inputs.npz or gene_holdout_pheno.npz under {score_dir}"
        )
    blob = np.load(pheno_path, allow_pickle=True)
    arrays = {
        "age": blob["age"],
        "age_mask": blob["age_mask"],
        "tissue": blob["tissue"],
        "tissue_mask": blob["tissue_mask"],
        "sex": blob["sex"],
        "sex_mask": blob["sex_mask"],
        "study_ids": blob["study_ids"],
    }
    class_names = (
        [str(x) for x in blob["class_names"].tolist()]
        if "class_names" in blob.files
        else ["tissue"]
    )
    return (
        np.asarray(blob["train_idx"], dtype=np.int64),
        np.asarray(blob["test_idx"], dtype=np.int64),
        arrays,
        class_names,
    )


def evaluate_gene_holdout_run(
    run_dir: Path,
    *,
    force: bool = False,
    fold: int | None = None,
) -> dict[str, Any]:
    score_dir = _resolve_score_dir(run_dir, fold)
    metrics_path = score_dir.parent / "metrics.json"
    holdout_meta_path = score_dir / "gene_holdout.json"
    mbs_held_path = score_dir / "mbs_heldout.npy"

    for p in (mbs_held_path, holdout_meta_path, metrics_path):
        if not p.is_file() and not (p.exists() if p.suffix == ".zarr" else False):
            raise FileNotFoundError(f"missing gene-holdout artifact: {p}")

    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    existing = ((metrics.get("evaluations") or {}).get("gene_holdout_nested")) or {}
    if existing and not force:
        print(  # noqa: T201
            f"[gene_holdout_nested] skip {run_dir.name}: already present", flush=True
        )
        return existing

    train_idx, test_idx, arrays, class_names = _load_split_arrays(score_dir, run_dir)
    mbs_train = _load_train_mbs(score_dir)
    mbs_held = np.load(mbs_held_path)
    if mbs_train.shape[0] != mbs_held.shape[0]:
        raise ValueError(
            f"row mismatch train MBS {mbs_train.shape} vs heldout {mbs_held.shape}"
        )

    train_out = _nested_on_mbs(
        mbs_train,
        train_idx=train_idx,
        test_idx=test_idx,
        arrays=arrays,
        class_names=class_names,
        evaluation="mbs_enet_nested_train_genes",
    )
    held_out = _nested_on_mbs(
        mbs_held,
        train_idx=train_idx,
        test_idx=test_idx,
        arrays=arrays,
        class_names=class_names,
        evaluation="mbs_enet_nested_heldout_genes",
    )
    holdout_meta = json.loads(holdout_meta_path.read_text(encoding="utf-8"))
    report = {
        "partition": {
            "method": holdout_meta.get("method"),
            "seed": holdout_meta.get("seed"),
            "heldout_fraction": holdout_meta.get("heldout_fraction"),
            "actual_heldout_fraction": holdout_meta.get("actual_heldout_fraction"),
            "heldout_chromosomes": holdout_meta.get("heldout_chromosomes"),
            "n_train_genes": holdout_meta.get("n_train_genes_scored")
            or len(holdout_meta.get("train_gene_ids") or []),
            "n_heldout_genes": holdout_meta.get("n_heldout_genes_scored")
            or len(holdout_meta.get("heldout_gene_ids") or []),
            "caveat": holdout_meta.get("caveat"),
        },
        "train_genes": _summary(train_out),
        "heldout_genes": _summary(held_out),
        "evaluations": {
            "mbs_enet_nested_train_genes": train_out,
            "mbs_enet_nested_heldout_genes": held_out,
        },
    }
    out_path = score_dir / "gene_holdout_eval.json"
    write_json(out_path, report)

    evaluations = dict(metrics.get("evaluations") or {})
    evaluations["gene_holdout_nested"] = {
        "train_genes": report["train_genes"],
        "heldout_genes": report["heldout_genes"],
        "partition": report["partition"],
        "artifact": str(out_path),
    }
    evaluations["mbs_enet_nested_train_genes"] = train_out
    evaluations["mbs_enet_nested_heldout_genes"] = held_out
    metrics["evaluations"] = evaluations
    write_json(metrics_path, metrics)

    tg = report["train_genes"]
    hg = report["heldout_genes"]
    print(  # noqa: T201
        f"[gene_holdout_nested] {run_dir.name} ({score_dir.parent.name})\n"
        f"  train_genes   tissue={tg['tissue_macro_f1']} age={tg['age_mae']} "
        f"sex={tg['sex_auroc']} n={tg['n_score_features']}\n"
        f"  heldout_genes tissue={hg['tissue_macro_f1']} age={hg['age_mae']} "
        f"sex={hg['sex_auroc']} n={hg['n_score_features']}",
        flush=True,
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-id",
        required=True,
        help="Flat run id or cascade run id under artifacts/runs/",
    )
    parser.add_argument("--fold", type=int, default=None, help="Cascade fold index")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    run_dir = paths.artifact_root / "runs" / args.run_id
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    evaluate_gene_holdout_run(run_dir, force=bool(args.force), fold=args.fold)


if __name__ == "__main__":
    main()
