#!/usr/bin/env python3
"""Milestone 10c: freeze-and-reuse disease/cancer probes on the finalized P2-G
encoder (DeepRVAT precedent -- freeze the gene-invariant encoder, fit
lightweight per-trait probes on frozen scores rather than joint retraining).

Neither disease nor cancer is in the current champion config's loss at all.
Real case/control labels already exist (12,194 / 9,077 usable samples) via
``sample_phenotype_table_hub_nine_pack_v2.parquet``'s ``disease_label_status``
/ ``cancer_label_status`` (case/control/adjacent_normal, built from
``SAMPLE_TYPE_CASE_CONTROL`` -- "unknown" is never treated as "control").

For each of the 3 study-grouped nine-pack folds: load that fold's own
converged P2-G checkpoint (no cross-fold reuse -- avoids leaking that fold's
train split into this fold's evaluation), score ALL nine-pack samples once
(frozen forward pass, no gradient), then fit a logistic-regression probe on
MBS (gene-pooled) and on RBS (region-level, pre-pooling) using only that
fold's own train_idx, scored on that fold's own test_idx.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import yaml
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score

from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_for_matrix,
    read_locus_index,
    read_sample_index,
)
from mbs.models import CascadeDeepSet
from mbs.paths import DataPaths
from mbs.training.cascade_assign import assignment_gene_linked_only, build_cascade_assignment
from mbs.training.cascade_loop import score_samples
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.locus_gene import load_graph_tables

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "experiment" / "stage0_7h_nine_pack_p2_g.yaml"
CHECKPOINT_RUN_ID = "stage0-7h-nine-pack-P2-G"
REPORT_DIR = (
    ROOT / "reports" / "inspection" / "stage0_7h_nine_pack_smoke" / "disease_cancer_frozen_probes"
)

TRAITS = {
    "disease": ("disease_label_status", "disease_case_control_mask"),
    "cancer": ("cancer_label_status", "cancer_case_control_mask"),
}


def _load_assignment(
    paths: DataPaths, cfg: dict[str, Any], max_loci: int
) -> tuple[Any, Path]:
    pilot = cfg["pilot"]
    graph_id = str(pilot["graph_id"])
    graph_dir = paths.data_root / "canonical" / "graphs" / graph_id
    matrix_id = str(pilot["matrix_id"])
    matrix_root = paths.data_root / "canonical" / "matrices" / matrix_id
    locus_index = read_locus_index(matrix_store_paths(matrix_root).locus_index_path)
    lr_edges, regions = load_graph_tables(graph_dir)
    genes_path = graph_dir / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation="explicit_only",
    )
    assignment = assignment_gene_linked_only(assignment)
    return assignment, matrix_root


def _load_model_from_checkpoint(
    ckpt_path: Path, n_region_types: int, device: torch.device
) -> CascadeDeepSet:
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    model = CascadeDeepSet(
        1,
        n_region_types,
        cpg_hidden_dim=int(ckpt["cpg_hidden_dim"]),
        region_hidden_dim=int(ckpt["region_hidden_dim"]),
        cpg_pool=ckpt["cpg_pool"],
        region_pool=ckpt["region_pool"],
        gene_aggregation=ckpt["gene_aggregation"],
        # Must match the training-time dropout rate -- SharedMLP only inserts
        # a Dropout layer when dropout > 0, which shifts state_dict key
        # indices (network.3 vs network.4). Harmless either way since
        # model.eval() makes nn.Dropout a no-op regardless of rate.
        dropout=0.1,
        activation="gelu",
        layer_norm=True,
    )
    model.load_state_dict(ckpt["model"])
    model.to(device)
    model.eval()
    return model


def _fit_eval_probe(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> dict[str, Any]:
    train_sel = train_idx[mask[train_idx]]
    test_sel = test_idx[mask[test_idx]]
    if train_sel.size < 20 or test_sel.size < 5:
        return {"n_train": int(train_sel.size), "n_test": int(test_sel.size), "skipped": True}
    y_train = y[train_sel]
    y_test = y[test_sel]
    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return {"n_train": int(train_sel.size), "n_test": int(test_sel.size), "skipped": True}
    clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.1)
    clf.fit(x[train_sel], y_train)
    proba = clf.predict_proba(x[test_sel])[:, 1]
    return {
        "n_train": int(train_sel.size),
        "n_test": int(test_sel.size),
        "auroc": float(roc_auc_score(y_test, proba)),
        "auprc": float(average_precision_score(y_test, proba)),
        "skipped": False,
    }


def _fit_eval_regression_probe(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> dict[str, Any]:
    train_sel = train_idx[mask[train_idx]]
    test_sel = test_idx[mask[test_idx]]
    if train_sel.size < 20 or test_sel.size < 5:
        return {"n_train": int(train_sel.size), "n_test": int(test_sel.size), "skipped": True}
    reg = Ridge(alpha=10.0)
    reg.fit(x[train_sel], y[train_sel])
    pred = reg.predict(x[test_sel])
    return {
        "n_train": int(train_sel.size),
        "n_test": int(test_sel.size),
        "mae": float(mean_absolute_error(y[test_sel], pred)),
        "r2": float(r2_score(y[test_sel], pred)),
        "skipped": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    max_loci = int(cfg["cv_budget"]["max_loci"])
    split_id = str(cfg["split_id"])
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")

    print(f"[disease-cancer-probe] loading assignment max_loci={max_loci}", flush=True)
    assignment, matrix_root = _load_assignment(paths, cfg, max_loci)
    n_region_types = max(len(assignment.region_types), 1)

    sample_index = read_sample_index(matrix_store_paths(matrix_root).sample_index_path)
    sample_ids_all = sample_index["sample_id"].astype(str).tolist()
    row_by_id = {sid: i for i, sid in enumerate(sample_ids_all)}

    print("[disease-cancer-probe] loading betas into RAM…", flush=True)
    betas_all = np.asarray(open_betas_for_matrix(matrix_root)[:, :max_loci], dtype=np.float32)
    print(f"[disease-cancer-probe] betas shape={betas_all.shape}", flush=True)

    pheno = pd.read_parquet(
        paths.data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_hub_nine_pack_v2.parquet"
    )
    pheno = pheno.set_index("sample_id")

    n = len(sample_ids_all)
    trait_arrays: dict[str, dict[str, np.ndarray]] = {}
    for trait, (status_col, mask_col) in TRAITS.items():
        y = np.zeros(n, dtype=np.int64)
        mask = np.zeros(n, dtype=bool)
        for i, sid in enumerate(sample_ids_all):
            if sid not in pheno.index:
                continue
            row = pheno.loc[sid]
            if bool(row.get(mask_col, False)) and row.get(status_col) in ("case", "control"):
                mask[i] = True
                y[i] = 1 if row[status_col] == "case" else 0
        trait_arrays[trait] = {"y": y, "mask": mask}
        print(
            f"[disease-cancer-probe] {trait}: usable={int(mask.sum())} case={int(y[mask].sum())}",
            flush=True,
        )

    bmi = pd.read_parquet(
        paths.data_root / "canonical" / "phenotypes" / "bmi_sample_info.parquet"
    ).set_index("sample_id")
    bmi_y = np.zeros(n, dtype=np.float64)
    bmi_mask = np.zeros(n, dtype=bool)
    for i, sid in enumerate(sample_ids_all):
        if sid in bmi.index:
            bmi_mask[i] = True
            bmi_y[i] = float(bmi.loc[sid, "phenotype_value_numeric"])
    print(f"[disease-cancer-probe] bmi: usable={int(bmi_mask.sum())}", flush=True)

    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)
    folds = fold_pack["folds"]

    results: dict[str, list[dict[str, Any]]] = {t: [] for t in (*TRAITS, "bmi")}
    for fold_i, fold in enumerate(folds):
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id]
        test_ids = [s for s in (fold.get("external_test_sample_ids") or []) if s in row_by_id]
        if not test_ids:
            test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id]
        train_idx = np.asarray([row_by_id[s] for s in train_ids], dtype=np.int64)
        test_idx = np.asarray([row_by_id[s] for s in test_ids], dtype=np.int64)

        ckpt_path = ROOT / "artifacts" / "runs" / CHECKPOINT_RUN_ID / f"fold_{fold_i}" / "best.pt"
        print(f"[disease-cancer-probe] fold {fold_i}: loading {ckpt_path}", flush=True)
        model = _load_model_from_checkpoint(ckpt_path, n_region_types, device)

        mbs, present, _orphan_rbs, all_rbs, all_rbs_present = score_samples(
            model, assignment, betas_all, device=device, batch_size=1
        )
        mbs_centered = np.where(present, mbs, 0.5)
        rbs_centered = np.where(all_rbs_present, all_rbs, 0.5)

        for trait in TRAITS:
            y = trait_arrays[trait]["y"]
            mask = trait_arrays[trait]["mask"]
            mbs_res = _fit_eval_probe(mbs_centered, y, mask, train_idx, test_idx)
            rbs_res = _fit_eval_probe(rbs_centered, y, mask, train_idx, test_idx)
            print(
                f"[disease-cancer-probe] fold {fold_i} {trait} mbs={mbs_res} rbs={rbs_res}",
                flush=True,
            )
            results[trait].append({"fold": fold_i, "mbs": mbs_res, "rbs": rbs_res})

        bmi_mbs_res = _fit_eval_regression_probe(mbs_centered, bmi_y, bmi_mask, train_idx, test_idx)
        bmi_rbs_res = _fit_eval_regression_probe(rbs_centered, bmi_y, bmi_mask, train_idx, test_idx)
        print(
            f"[disease-cancer-probe] fold {fold_i} bmi mbs={bmi_mbs_res} rbs={bmi_rbs_res}",
            flush=True,
        )
        results["bmi"].append({"fold": fold_i, "mbs": bmi_mbs_res, "rbs": bmi_rbs_res})

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = {
        "generated": datetime.now(UTC).isoformat(),
        "checkpoint_run_id": CHECKPOINT_RUN_ID,
        "split_id": split_id,
        "results": results,
    }
    (REPORT_DIR / "results.json").write_text(json.dumps(out, indent=2))
    print(f"[disease-cancer-probe] wrote {REPORT_DIR / 'results.json'}", flush=True)


if __name__ == "__main__":
    main()
