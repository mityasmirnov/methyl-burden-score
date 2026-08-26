"""Milestone 7E development CV: shared folds, arm matrix, report."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from mbs.annotation.manifest import write_json
from mbs.evaluation.splits import assert_no_study_leakage, build_outer_study_grouped_folds
from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)
from mbs.training.branch import train_branch_arm
from mbs.training.controls import fit_metadata_only
from mbs.training.hier_loop import train_hierarchical_baseline
from mbs.training.late_fusion import concatenate_score_blocks, evaluate_late_fusion
from mbs.training.locus_gene import build_locus_gene_index, load_graph_tables
from mbs.training.loop import train_flat_baseline
from mbs.training.phenotype_table import load_tissue_ontology
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.run_artifacts import run_dir
from mbs.training.transparent_baselines import (
    fit_elasticnet_phenotype,
    presence_aware_means,
    run_elasticnet_baseline,
    run_mean_baseline,
)

DEFAULT_SPLIT_ID = "hub-ats-7e-3fold-v1"
DEFAULT_RESTART_SEEDS = (42, 43)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def freeze_outer_folds(
    samples: list[dict[str, Any]],
    *,
    out_dir: Path,
    n_folds: int = 3,
    seed: int = 42,
    val_fraction: float = 0.15,
    split_id: str = DEFAULT_SPLIT_ID,
) -> dict[str, Any]:
    """Build and persist the shared 7E outer-fold manifest (do not touch 5d freeze)."""
    pack = build_outer_study_grouped_folds(
        samples,
        n_folds=n_folds,
        seed=seed,
        val_fraction=val_fraction,
        split_id=split_id,
    )
    for fold in pack["folds"]:
        assert_no_study_leakage(fold)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "folds.json"
    text = json.dumps(pack, indent=2, sort_keys=True, default=str)
    path.write_text(text + "\n", encoding="utf-8")
    digest = _sha256_bytes(text.encode("utf-8"))
    write_json(out_dir / "folds.sha256.json", {"path": str(path), "sha256": digest})
    pack["sha256"] = digest
    pack["path"] = str(path)
    return pack


def load_frozen_folds(path: Path) -> dict[str, Any]:
    """Load a previously frozen outer-fold pack."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if "folds" not in payload:
        raise ValueError(f"not an outer-fold pack: {path}")
    for fold in payload["folds"]:
        assert_no_study_leakage(fold)
    return payload


def samples_from_phenotype_table(
    table_path: Path,
    *,
    ontology_path: Path,
) -> tuple[list[dict[str, Any]], list[Any]]:
    """Build split sample rows + multitask phenotypes from Hub ATS table."""
    ontology = load_tissue_ontology(ontology_path)
    table = pd.read_parquet(table_path)
    sample_ids = table["sample_id"].astype(str).tolist()
    phenotypes, _class_names = load_multitask_phenotypes(
        table_path,
        sample_ids=sample_ids,
        class_names=ontology.class_names,
    )
    rows = [
        {
            "sample_id": p.sample_id,
            "study_id": str(p.study_id or p.sample_id),
            "platform": p.platform,
            "group_id": str(p.study_id or p.sample_id),
        }
        for p in phenotypes
    ]
    return rows, phenotypes


def inject_fold_into_config(
    config: dict[str, Any],
    fold: dict[str, Any],
    *,
    seed: int,
) -> dict[str, Any]:
    """Copy config and pin study lists to one outer fold (disable auto_split)."""
    cfg = deepcopy(config)
    cfg.setdefault("experiment", {})["seed"] = int(seed)
    pilot = cfg.setdefault("pilot", {})
    pilot["auto_split"] = False
    pilot["reuse_flat_split"] = False
    pilot["train_studies"] = list(fold["train_studies"])
    pilot["validation_studies"] = list(fold["validation_studies"])
    pilot["external_test_studies"] = list(fold["external_test_studies"])
    pilot["split_id"] = str(fold.get("split_id") or fold.get("outer_fold"))
    return cfg


def arm_run_id(arm: str, fold_idx: int, restart: int, *, tag: str = "") -> str:
    suffix = f"-{tag}" if tag else ""
    return f"stage0-7e-{arm}-f{fold_idx}-r{restart}{suffix}"


def _metrics_done(run_root: Path) -> bool:
    return (run_root / "metrics.json").is_file()


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _phenotype_arrays(phenotypes: list[Any], sample_ids: list[str]) -> dict[str, Any]:
    by_id = {p.sample_id: p for p in phenotypes}
    ordered = [by_id[s] for s in sample_ids if s in by_id]
    return {
        "sample_ids": [p.sample_id for p in ordered],
        "age": np.array([float(p.age or 0.0) for p in ordered], dtype=np.float64),
        "age_mask": np.array([bool(p.age_mask) for p in ordered], dtype=bool),
        "tissue": np.array([int(p.class_index) for p in ordered], dtype=np.int64),
        "tissue_mask": np.array([bool(p.tissue_mask) for p in ordered], dtype=bool),
        "sex": np.array([int(p.sex_class_index or 0) for p in ordered], dtype=np.int64),
        "sex_mask": np.array([bool(p.sex_mask) for p in ordered], dtype=bool),
        "study_ids": np.array([str(p.study_id or p.sample_id) for p in ordered], dtype=object),
        "platforms": np.array([str(p.platform or "unknown") for p in ordered], dtype=object),
        "tissues": np.array([str(p.cell_type or "unknown") for p in ordered], dtype=object),
    }


def run_transparent_on_synthetic(
    *,
    fold: dict[str, Any],
    seed: int,
) -> dict[str, Any]:
    """Fixture-scale transparent baselines (unit / smoke; not Hub Done-when)."""
    rng = np.random.default_rng(seed)
    n_train = max(8, len(fold["train_sample_ids"]))
    n_test = max(4, len(fold.get("external_test_sample_ids") or fold["validation_sample_ids"]))
    n_genes = 6
    n_loci = 12
    gene_index = np.array([i % n_genes for i in range(n_loci)], dtype=np.int64)
    betas_tr = rng.uniform(0.1, 0.9, size=(n_train, n_loci))
    obs_tr = rng.random((n_train, n_loci)) > 0.1
    betas_te = rng.uniform(0.1, 0.9, size=(n_test, n_loci))
    obs_te = rng.random((n_test, n_loci)) > 0.1
    x_tr = presence_aware_means(betas_tr, obs_tr, gene_index, n_groups=n_genes)
    x_te = presence_aware_means(betas_te, obs_te, gene_index, n_groups=n_genes)
    age_tr = 20.0 + 5.0 * x_tr.mean(axis=1) + rng.normal(0, 0.1, size=n_train)
    age_te = 20.0 + 5.0 * x_te.mean(axis=1) + rng.normal(0, 0.1, size=n_test)
    tissue_tr = (x_tr[:, 0] > 0.5).astype(np.int64)
    tissue_te = (x_te[:, 0] > 0.5).astype(np.int64)
    mean_out = run_mean_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=age_tr,
        age_mask_train=np.ones(n_train, dtype=bool),
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_train, dtype=bool),
        sex_train=None,
        sex_mask_train=None,
        age_test=age_te,
        age_mask_test=np.ones(n_test, dtype=bool),
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_test, dtype=bool),
        sex_test=None,
        sex_mask_test=None,
        study_ids_test=np.array(["S0"] * n_test, dtype=object),
        kind="gene",
    )
    enet = run_elasticnet_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=age_tr,
        age_mask_train=np.ones(n_train, dtype=bool),
        tissue_train=tissue_tr,
        tissue_mask_train=np.ones(n_train, dtype=bool),
        age_test=age_te,
        age_mask_test=np.ones(n_test, dtype=bool),
        tissue_test=tissue_te,
        tissue_mask_test=np.ones(n_test, dtype=bool),
        study_ids_test=np.array(["S0"] * n_test, dtype=object),
    )
    return {"mean_gene": mean_out, "elasticnet": enet}


def run_metadata_control(
    phenotypes: list[Any],
    fold: dict[str, Any],
) -> dict[str, Any]:
    """Metadata-only ceiling on the held-out fold (same folds as neural arms)."""
    test_ids = set(fold.get("external_test_sample_ids") or [])
    held = [p for p in phenotypes if p.sample_id in test_ids]
    if not held:
        held = [p for p in phenotypes if p.sample_id in set(fold["validation_sample_ids"])]
    out: dict[str, Any] = {}
    age_m = [p for p in held if p.age_mask and p.age is not None]
    if age_m:
        out["age"] = fit_metadata_only(
            study_ids=[str(p.study_id or p.sample_id) for p in age_m],
            platforms=[p.platform for p in age_m],
            tissues=[p.cell_type for p in age_m],
            y=np.array([float(p.age) for p in age_m]),
            task="regression",
        )
    tissue_m = [p for p in held if p.tissue_mask]
    if tissue_m:
        out["tissue"] = fit_metadata_only(
            study_ids=[str(p.study_id or p.sample_id) for p in tissue_m],
            platforms=[p.platform for p in tissue_m],
            tissues=[p.cell_type for p in tissue_m],
            y=np.array([p.class_index for p in tissue_m]),
            task="multiclass",
        )
    return out


# Arms that may win Milestone 7 (exclude intermediate RBS/TBS branch trains).
_WINNER_ARM_PREFIXES = (
    "N-flat-gene",
    "N-hier-gene",
    "N-gene-direct",
    "N-multipath",
)


def _is_winner_candidate(arm: str) -> bool:
    return any(arm.startswith(p) for p in _WINNER_ARM_PREFIXES)


def _select_winner(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Pick architecture by mean tissue macro-F1 then mean age RMSE across folds."""
    by_arm: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if str(row.get("family")) != "neural":
            continue
        arm = str(row.get("arm") or "")
        if not _is_winner_candidate(arm):
            continue
        by_arm.setdefault(arm, []).append(row)
    if not by_arm:
        return {"arm": None, "reason": "no completed neural architecture arms"}

    def aggregate(arm: str) -> dict[str, Any]:
        items = by_arm[arm]
        f1s = [float(x["tissue_macro_f1"]) for x in items if x.get("tissue_macro_f1") is not None]
        maes = [float(x["age_mae"]) for x in items if x.get("age_mae") is not None]
        rmses = [float(x["age_rmse"]) for x in items if x.get("age_rmse") is not None]
        return {
            "arm": arm,
            "tissue_macro_f1": float(np.mean(f1s)) if f1s else None,
            "age_mae": float(np.mean(maes)) if maes else None,
            "age_rmse": float(np.mean(rmses)) if rmses else None,
            "n_cells": len(items),
        }

    scored = [aggregate(a) for a in by_arm]
    scored.sort(
        key=lambda r: (
            -(r["tissue_macro_f1"] if r["tissue_macro_f1"] is not None else -1.0),
            # Prefer destandardized MAE (years) when present; RMSE units differ
            # across neural (z) vs late-fusion (years).
            r["age_mae"]
            if r["age_mae"] is not None
            else (r["age_rmse"] if r["age_rmse"] is not None else 1e9),
        )
    )
    best = scored[0]
    best["reason"] = (
        "max mean tissue macro-F1, then min mean age MAE (years) among architecture arms"
    )
    return best


def run_hub_transparent_arm(
    *,
    data_root: Path,
    fold: dict[str, Any],
    phenotypes: list[Any],
    kind: str,
    max_loci: int | None,
    seed: int,
) -> dict[str, Any]:
    """Gene-mean (or elastic-net on means) baseline on Hub ATS fold."""
    matrix_id = "matrix-hub-age-tissue-sex-full-v1"
    graph_id = "graph-grch38-gencode38-five-role-v1"
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
    locus_gene = build_locus_gene_index(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        max_loci=max_loci,
        region_systems=("gene",),
    )
    betas = open_betas_zarr(matrix_paths.betas_path)
    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }
    train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id]
    test_ids = [s for s in (fold.get("external_test_sample_ids") or []) if s in row_by_id]
    if not test_ids:
        test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id]
    n_cols = (
        locus_gene.n_study_loci if max_loci is None else min(int(max_loci), locus_gene.n_study_loci)
    )
    cols = locus_gene.edge_col_index
    genes = locus_gene.edge_gene_index
    # Compact locus→gene map for first n_cols study columns used in edges.
    keep = cols < n_cols
    cols_k = cols[keep]
    genes_k = genes[keep]

    def _matrix_for(ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
        rows = []
        obs = []
        for sid in ids:
            beta_row = np.asarray(betas[row_by_id[sid], :n_cols], dtype=np.float64)
            vals = beta_row[cols_k]
            finite = np.isfinite(vals)
            rows.append(np.where(finite, vals, 0.5))
            obs.append(finite)
        return np.stack(rows, axis=0), np.stack(obs, axis=0)

    betas_tr, obs_tr = _matrix_for(train_ids)
    betas_te, obs_te = _matrix_for(test_ids)
    x_tr = presence_aware_means(
        betas_tr, obs_tr, genes_k, n_groups=locus_gene.n_genes, empty_fill=0.5
    )
    x_te = presence_aware_means(
        betas_te, obs_te, genes_k, n_groups=locus_gene.n_genes, empty_fill=0.5
    )
    ph_tr = _phenotype_arrays(phenotypes, train_ids)
    ph_te = _phenotype_arrays(phenotypes, test_ids)
    if kind == "elasticnet":
        out = run_elasticnet_baseline(
            x_train=x_tr,
            x_test=x_te,
            age_train=ph_tr["age"],
            age_mask_train=ph_tr["age_mask"],
            tissue_train=ph_tr["tissue"],
            tissue_mask_train=ph_tr["tissue_mask"],
            age_test=ph_te["age"],
            age_mask_test=ph_te["age_mask"],
            tissue_test=ph_te["tissue"],
            tissue_mask_test=ph_te["tissue_mask"],
            study_ids_test=ph_te["study_ids"],
            platforms_test=ph_te["platforms"],
        )
    else:
        out = run_mean_baseline(
            x_train=x_tr,
            x_test=x_te,
            age_train=ph_tr["age"],
            age_mask_train=ph_tr["age_mask"],
            tissue_train=ph_tr["tissue"],
            tissue_mask_train=ph_tr["tissue_mask"],
            sex_train=ph_tr["sex"],
            sex_mask_train=ph_tr["sex_mask"],
            age_test=ph_te["age"],
            age_mask_test=ph_te["age_mask"],
            tissue_test=ph_te["tissue"],
            tissue_mask_test=ph_te["tissue_mask"],
            sex_test=ph_te["sex"],
            sex_mask_test=ph_te["sex_mask"],
            study_ids_test=ph_te["study_ids"],
            platforms_test=ph_te["platforms"],
            kind="gene",
        )
    out["seed"] = seed
    out["n_train"] = len(train_ids)
    out["n_test"] = len(test_ids)
    out["max_loci"] = n_cols
    return out


def run_hub_late_fusion_arm(
    *,
    data_root: Path,
    fold: dict[str, Any],
    phenotypes: list[Any],
    max_loci: int | None,
    seed: int,
    systems: tuple[str, ...] = ("gene", "rbs", "tbs"),
    include_direct: bool = False,
) -> dict[str, Any]:
    """Late-fuse presence-aware means (+ optional direct enet preds) with linear heads.

    Independent neural branch trains remain separate arms; this supplies the
    multipath / gene+direct fusion cell on the same frozen folds (graph-v2).
    """
    matrix_id = "matrix-hub-age-tissue-sex-full-v1"
    graph_id = "graph-grch38-gencode38-cgi-tile-v2"
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
    betas = open_betas_zarr(matrix_paths.betas_path)
    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }
    train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id]
    test_ids = [s for s in (fold.get("external_test_sample_ids") or []) if s in row_by_id]
    if not test_ids:
        test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id]
    n_cols = int(max_loci) if max_loci is not None else int(locus_index.shape[0])
    n_cols = min(n_cols, int(locus_index.shape[0]))

    def _betas_for(ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
        rows = []
        obs = []
        for sid in ids:
            beta_row = np.asarray(betas[row_by_id[sid], :n_cols], dtype=np.float64)
            finite = np.isfinite(beta_row)
            rows.append(np.where(finite, beta_row, 0.5))
            obs.append(finite)
        return np.stack(rows, axis=0), np.stack(obs, axis=0)

    betas_tr, obs_tr = _betas_for(train_ids)
    betas_te, obs_te = _betas_for(test_ids)
    blocks_tr: list[np.ndarray] = []
    blocks_te: list[np.ndarray] = []
    for system in systems:
        locus_gene = build_locus_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=n_cols,
            region_systems=(system,),
        )
        cols = locus_gene.edge_col_index
        groups = locus_gene.edge_gene_index
        keep = cols < n_cols
        cols_k = cols[keep]
        groups_k = groups[keep]
        if cols_k.size == 0 or locus_gene.n_genes == 0:
            continue
        blocks_tr.append(
            presence_aware_means(
                betas_tr[:, cols_k],
                obs_tr[:, cols_k],
                groups_k,
                n_groups=locus_gene.n_genes,
                empty_fill=0.5,
            )
        )
        blocks_te.append(
            presence_aware_means(
                betas_te[:, cols_k],
                obs_te[:, cols_k],
                groups_k,
                n_groups=locus_gene.n_genes,
                empty_fill=0.5,
            )
        )
    if not blocks_tr:
        raise ValueError(f"late fusion found no loci for systems={systems}")
    x_tr = concatenate_score_blocks(blocks_tr)
    x_te = concatenate_score_blocks(blocks_te)
    ph_tr = _phenotype_arrays(phenotypes, train_ids)
    ph_te = _phenotype_arrays(phenotypes, test_ids)
    if include_direct:
        # Direct path: append train-fold elastic-net predictions as score columns.
        age_tr_col = np.zeros((x_tr.shape[0], 1), dtype=np.float32)
        age_te_col = np.zeros((x_te.shape[0], 1), dtype=np.float32)
        if ph_tr["age"] is not None and bool(np.asarray(ph_tr["age_mask"]).any()):
            m = np.asarray(ph_tr["age_mask"], dtype=bool)
            age_model = fit_elasticnet_phenotype(
                x_tr[m], np.asarray(ph_tr["age"])[m], task="regression"
            )
            age_tr_col[:, 0] = age_model.predict(x_tr).astype(np.float32)
            age_te_col[:, 0] = age_model.predict(x_te).astype(np.float32)
        tissue_tr_col = np.zeros((x_tr.shape[0], 1), dtype=np.float32)
        tissue_te_col = np.zeros((x_te.shape[0], 1), dtype=np.float32)
        if ph_tr["tissue"] is not None and bool(np.asarray(ph_tr["tissue_mask"]).any()):
            m = np.asarray(ph_tr["tissue_mask"], dtype=bool)
            tissue_model = fit_elasticnet_phenotype(
                x_tr[m], np.asarray(ph_tr["tissue"])[m], task="multiclass"
            )
            tissue_tr_col[:, 0] = tissue_model.predict(x_tr).astype(np.float32)
            tissue_te_col[:, 0] = tissue_model.predict(x_te).astype(np.float32)
        x_tr = concatenate_score_blocks([x_tr, age_tr_col, tissue_tr_col])
        x_te = concatenate_score_blocks([x_te, age_te_col, tissue_te_col])
    out = evaluate_late_fusion(
        scores_train=x_tr,
        scores_test=x_te,
        age_train=ph_tr["age"],
        age_mask_train=ph_tr["age_mask"],
        tissue_train=ph_tr["tissue"],
        tissue_mask_train=ph_tr["tissue_mask"],
        sex_train=ph_tr["sex"],
        sex_mask_train=ph_tr["sex_mask"],
        age_test=ph_te["age"],
        age_mask_test=ph_te["age_mask"],
        tissue_test=ph_te["tissue"],
        tissue_mask_test=ph_te["tissue_mask"],
        sex_test=ph_te["sex"],
        sex_mask_test=ph_te["sex_mask"],
        study_ids_test=ph_te["study_ids"],
        platforms_test=ph_te["platforms"],
    )
    out["late_fusion"] = True
    out["systems"] = list(systems)
    out["include_direct"] = include_direct
    out["seed"] = seed
    out["n_train"] = len(train_ids)
    out["n_test"] = len(test_ids)
    out["max_loci"] = n_cols
    return out


def _extract_holdout_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    # Transparent / late-fusion nested under metrics.metrics
    nested = metrics.get("metrics") if isinstance(metrics.get("metrics"), dict) else None
    if nested and ("age" in nested or "tissue" in nested):
        age = nested.get("age") or {}
        tissue = nested.get("tissue") or {}
        return {
            "tissue_macro_f1": tissue.get("macro_f1"),
            "tissue_balanced_accuracy": tissue.get("balanced_accuracy"),
            "age_rmse": age.get("rmse"),
            "age_mae": age.get("mae"),
            "age_r2": age.get("r2"),
            "age_pearson_r": age.get("pearson_r"),
        }
    ext = metrics.get("external_test") or {}
    if isinstance(ext, dict) and "age" in ext and isinstance(ext.get("age"), dict):
        age = ext.get("age") or {}
        tissue = ext.get("tissue") or {}
        return {
            "tissue_macro_f1": tissue.get("macro_f1"),
            "tissue_balanced_accuracy": tissue.get("balanced_accuracy"),
            "age_rmse": age.get("rmse"),
            "age_mae": age.get("mae"),
            "age_r2": age.get("r2"),
            "age_pearson_r": age.get("pearson_r"),
        }
    return {
        "tissue_macro_f1": ext.get("macro_f1"),
        "tissue_balanced_accuracy": ext.get("balanced_accuracy"),
        "age_rmse": ext.get("rmse"),
        "age_mae": ext.get("mae"),
        "age_r2": ext.get("r2"),
        "age_pearson_r": ext.get("pearson_r"),
        "auroc": ext.get("auroc"),
        "auprc": ext.get("auprc"),
    }


def write_dev_cv_report(
    *,
    report_dir: Path,
    fold_pack: dict[str, Any],
    results: list[dict[str, Any]],
    metadata_controls: list[dict[str, Any]],
) -> dict[str, Any]:
    """Aggregate arm metrics and name a Milestone 7 winner."""
    report_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for item in results:
        metrics = item.get("metrics") or {}
        hold = _extract_holdout_metrics(metrics if isinstance(metrics, dict) else {})
        rows.append(
            {
                "arm": item.get("arm"),
                "family": item.get("family", "neural"),
                "fold": item.get("fold"),
                "restart": item.get("restart"),
                "run_id": item.get("run_id"),
                "level1": item.get("level1"),
                "cpgpt": item.get("cpgpt"),
                **hold,
            }
        )
    winner = _select_winner(rows)
    summary = {
        "milestone": "7E",
        "split_id": fold_pack.get("split_id"),
        "n_folds": fold_pack.get("n_folds"),
        "fold_sha256": fold_pack.get("sha256"),
        "n_results": len(rows),
        "winner": winner,
        "results": rows,
        "metadata_controls": metadata_controls,
        "selection_rule": (
            "Among neural arms: highest mean tissue macro-F1 on held-out studies, "
            "ties broken by lowest age RMSE. Transparent and metadata-only are ceilings."
        ),
    }
    write_json(report_dir / "summary.json", summary)
    lines = [
        "# Milestone 7E development CV",
        "",
        f"Split: `{summary['split_id']}` ({summary['n_folds']} outer folds).",
        "",
        f"**Winner for Milestone 7:** `{winner.get('arm')}` "
        f"(tissue macro-F1={winner.get('tissue_macro_f1')}, "
        f"age RMSE={winner.get('age_rmse')}).",
        "",
        f"Selection: {summary['selection_rule']}",
        "",
        "## Arms",
        "",
    ]
    lines.extend(
        [
            f"- `{row.get('arm')}` fold={row.get('fold')} restart={row.get('restart')} "
            f"L1={row.get('level1')} CpGPT={row.get('cpgpt')}: "
            f"tissue_f1={row.get('tissue_macro_f1')} age_rmse={row.get('age_rmse')}"
            for row in rows
        ]
    )
    lines.extend(["", "## Metadata-only controls", ""])
    lines.extend(
        f"- fold={mc.get('fold')}: `{json.dumps(mc.get('metrics'), default=str)}`"
        for mc in metadata_controls
    )
    (report_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def _apply_late_fusion_overlay(
    *,
    metrics: dict[str, Any],
    arm_name: str,
    fusion_from: Any,
    data_root: Path,
    fold: dict[str, Any],
    phenotypes: list[Any],
    max_loci: int | None,
    seed: int,
    run_root: Path,
) -> dict[str, Any]:
    """Replace neural metrics with late-fusion phenotype metrics; keep neural under key."""
    if arm_name.startswith("N-multipath") or fusion_from:
        systems: tuple[str, ...] = ("gene", "rbs", "tbs")
        include_direct = True
    else:
        systems = ("gene",)
        include_direct = True
    fused = run_hub_late_fusion_arm(
        data_root=data_root,
        fold=fold,
        phenotypes=phenotypes,
        max_loci=max_loci,
        seed=seed,
        systems=systems,
        include_direct=include_direct,
    )
    fused["neural_branch_metrics"] = metrics
    if fusion_from:
        fused["fusion_from"] = list(fusion_from)
    run_root.mkdir(parents=True, exist_ok=True)
    write_json(run_root / "metrics.json", fused)
    write_json(run_root / "neural_branch_metrics.json", metrics)
    return fused


def run_dev_cv(
    *,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    bakeoff_config: Path,
    report_dir: Path,
    device: str = "cuda",
    max_epochs: int | None = None,
    max_loci: int | None = None,
    fixture: bool = False,
) -> dict[str, Any]:
    """Orchestrate the 7E arm matrix with resume on completed run metrics."""
    spec = _load_yaml(bakeoff_config)
    split_id = str(spec.get("split_id", DEFAULT_SPLIT_ID))
    folds_dir = Path(str(spec.get("folds_dir") or (artifact_root / "splits" / split_id)))
    if not folds_dir.is_absolute():
        folds_dir = artifact_root / folds_dir
    restart_seeds = [int(x) for x in spec.get("restart_seeds", list(DEFAULT_RESTART_SEEDS))]
    n_folds = int(spec.get("n_folds", 3))
    raw_budget = spec.get("cv_budget")
    budget: dict[str, Any] = raw_budget if isinstance(raw_budget, dict) else {}
    if max_epochs is None and budget.get("max_epochs") is not None:
        max_epochs = int(budget["max_epochs"])
    if max_loci is None and budget.get("max_loci") is not None:
        max_loci = int(budget["max_loci"])

    table_rel = Path(
        str(
            spec.get(
                "sample_phenotype_table",
                "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet",
            )
        )
    )
    ont_rel = Path(
        str(
            spec.get(
                "tissue_ontology",
                "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml",
            )
        )
    )
    table_path = table_rel if table_rel.is_absolute() else data_root / table_rel
    ont_path = ont_rel if ont_rel.is_absolute() else data_root / ont_rel

    if fixture:
        # Synthetic studies for plumbing only.
        samples = [
            {"sample_id": f"s{i}", "study_id": f"GSE{i % 9}", "platform": "HM450"}
            for i in range(90)
        ]
        phenotypes: list[Any] = []
        fold_pack = freeze_outer_folds(
            samples, out_dir=folds_dir, n_folds=n_folds, seed=42, split_id=split_id
        )
    else:
        samples, phenotypes = samples_from_phenotype_table(table_path, ontology_path=ont_path)
        folds_path = folds_dir / "folds.json"
        if folds_path.is_file():
            fold_pack = load_frozen_folds(folds_path)
            fold_pack["sha256"] = json.loads(
                (folds_dir / "folds.sha256.json").read_text(encoding="utf-8")
            ).get("sha256")
        else:
            fold_pack = freeze_outer_folds(
                samples, out_dir=folds_dir, n_folds=n_folds, seed=42, split_id=split_id
            )

    results: list[dict[str, Any]] = []
    metadata_controls: list[dict[str, Any]] = []
    arms = list(spec.get("arms") or [])

    for fold_idx, fold in enumerate(fold_pack["folds"]):
        if phenotypes:
            metadata_controls.append(
                {"fold": fold_idx, "metrics": run_metadata_control(phenotypes, fold)}
            )
        for restart_i, seed in enumerate(restart_seeds):
            if fixture:
                results.append(
                    {
                        "arm": "T-mean-gene",
                        "family": "transparent",
                        "fold": fold_idx,
                        "restart": restart_i,
                        "run_id": arm_run_id("T-mean-gene", fold_idx, restart_i, tag="fixture"),
                        "level1": False,
                        "cpgpt": False,
                        "metrics": run_transparent_on_synthetic(fold=fold, seed=seed),
                    }
                )
                continue
            for arm_spec in arms:
                arm_name = str(arm_spec["name"])
                family = str(arm_spec.get("family", "neural"))
                rid = arm_run_id(
                    arm_name,
                    fold_idx,
                    restart_i,
                    tag=str(arm_spec.get("tag") or ""),
                )
                run_root = run_dir(artifact_root, rid)
                fusion_from = arm_spec.get("fusion_from")
                tag = str(arm_spec.get("tag") or "")
                wants_late_fusion = bool(fusion_from) or (
                    arm_name.startswith("N-gene-direct") and tag == "direct"
                )

                if _metrics_done(run_root):
                    metrics = json.loads((run_root / "metrics.json").read_text(encoding="utf-8"))
                    if wants_late_fusion and not metrics.get("late_fusion"):
                        # Neural branch finished earlier; apply fusion overlay.
                        metrics = _apply_late_fusion_overlay(
                            metrics=metrics,
                            arm_name=arm_name,
                            fusion_from=fusion_from,
                            data_root=data_root,
                            fold=fold,
                            phenotypes=phenotypes,
                            max_loci=max_loci,
                            seed=seed,
                            run_root=run_root,
                        )
                    results.append(
                        {
                            "arm": arm_name,
                            "family": family,
                            "fold": fold_idx,
                            "restart": restart_i,
                            "run_id": rid,
                            "level1": arm_spec.get("level1"),
                            "cpgpt": arm_spec.get("cpgpt"),
                            "metrics": metrics,
                            "resumed": True,
                        }
                    )
                    continue

                if family == "control":
                    continue

                if family == "transparent":
                    # Means/enet are deterministic given fold features; one restart.
                    if restart_i > 0:
                        continue
                    kind = "elasticnet" if "enet" in arm_name.lower() else "gene"
                    metrics = run_hub_transparent_arm(
                        data_root=data_root,
                        fold=fold,
                        phenotypes=phenotypes,
                        kind=kind,
                        max_loci=max_loci,
                        seed=seed,
                    )
                    run_root.mkdir(parents=True, exist_ok=True)
                    write_json(run_root / "metrics.json", metrics)
                    results.append(
                        {
                            "arm": arm_name,
                            "family": family,
                            "fold": fold_idx,
                            "restart": restart_i,
                            "run_id": rid,
                            "level1": arm_spec.get("level1"),
                            "cpgpt": False,
                            "metrics": metrics,
                        }
                    )
                    continue

                cfg_path = project_root / str(arm_spec["config"])
                cfg = inject_fold_into_config(_load_yaml(cfg_path), fold, seed=seed)
                if max_epochs is not None:
                    cfg.setdefault("training", {})["max_epochs"] = int(max_epochs)
                model_type = str(cfg.get("model", {}).get("type", "flat_deepset"))
                branch_arm = arm_spec.get("branch_arm")
                if branch_arm:
                    out = train_branch_arm(
                        arm=str(branch_arm),
                        project_root=project_root,
                        data_root=data_root,
                        artifact_root=artifact_root,
                        config=cfg,
                        run_id=rid,
                        device=device,
                        overfit_fixture=False,
                        max_epochs=max_epochs,
                        max_loci=max_loci,
                    )
                    metrics = out.get("metrics") or out
                elif model_type == "hierarchical_deepset":
                    result = train_hierarchical_baseline(
                        project_root=project_root,
                        data_root=data_root,
                        artifact_root=artifact_root,
                        config=cfg,
                        run_id=rid,
                        device_str=device,
                        max_epochs=max_epochs,
                        max_loci=max_loci,
                    )
                    metrics = result.metrics
                else:
                    result = train_flat_baseline(
                        project_root=project_root,
                        data_root=data_root,
                        artifact_root=artifact_root,
                        config=cfg,
                        run_id=rid,
                        device_str=device,
                        max_epochs=max_epochs,
                        max_loci=max_loci,
                    )
                    metrics = result.metrics

                if wants_late_fusion:
                    metrics = _apply_late_fusion_overlay(
                        metrics=metrics,
                        arm_name=arm_name,
                        fusion_from=fusion_from,
                        data_root=data_root,
                        fold=fold,
                        phenotypes=phenotypes,
                        max_loci=max_loci,
                        seed=seed,
                        run_root=run_root,
                    )

                results.append(
                    {
                        "arm": arm_name,
                        "family": family,
                        "fold": fold_idx,
                        "restart": restart_i,
                        "run_id": rid,
                        "level1": arm_spec.get("level1"),
                        "cpgpt": arm_spec.get("cpgpt"),
                        "metrics": metrics,
                    }
                )
                if fusion_from:
                    results[-1]["fusion_children"] = list(fusion_from)

    return write_dev_cv_report(
        report_dir=report_dir,
        fold_pack=fold_pack,
        results=results,
        metadata_controls=metadata_controls,
    )
