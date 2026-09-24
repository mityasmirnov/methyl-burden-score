#!/usr/bin/env python3
"""Milestone 10c follow-ups on the frozen P2-G encoder:

1. Alzheimer's vs. control probe (disease_sample_info's ``phenotype_value``
   == "Alzheimer's disease", 945 samples -- the only individual diagnosis
   clearing the n>=600 threshold) -- matched against the broader disease
   pack's ``label_status=control`` pool, not a same-study-only comparison.
2. Cancer subtype multiclass probe (predict which subtype, among cases
   only, restricted to subtypes with >=200 samples) -- tests whether
   subtypes are distinguishable from each other, not just cancer-vs-normal.
3. Case/control homogeneity diagnostic for disease and cancer: within-group
   feature variance (do controls cluster tighter than cases? do cases
   cluster tighter than each other?) vs. between-group separation. Explains
   *why* a case/control probe works or doesn't, independent of any specific
   classifier.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, f1_score, roc_auc_score

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
MIN_SUBTYPE_N = 200


def _load_assignment(paths: DataPaths, cfg: dict[str, Any], max_loci: int) -> tuple[Any, Path]:
    pilot = cfg["pilot"]
    graph_dir = paths.data_root / "canonical" / "graphs" / str(pilot["graph_id"])
    matrix_root = paths.data_root / "canonical" / "matrices" / str(pilot["matrix_id"])
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
    return assignment_gene_linked_only(assignment), matrix_root


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
        dropout=0.1,  # must match training-time rate; harmless under eval()
        activation="gelu",
        layer_norm=True,
    )
    model.load_state_dict(ckpt["model"])
    model.to(device)
    model.eval()
    return model


def _homogeneity(x: np.ndarray, case_idx: np.ndarray, control_idx: np.ndarray) -> dict[str, float]:
    within_case = float(np.mean(np.var(x[case_idx], axis=0)))
    within_control = float(np.mean(np.var(x[control_idx], axis=0)))
    between = float(np.mean((x[case_idx].mean(axis=0) - x[control_idx].mean(axis=0)) ** 2))
    return {
        "within_case_var": within_case,
        "within_control_var": within_control,
        "between_group_var": between,
        "case_more_variable_than_control": within_case > within_control,
        "separation_over_control_spread": between / within_control if within_control > 0 else None,
        "separation_over_case_spread": between / within_case if within_case > 0 else None,
    }


def main() -> None:
    paths = DataPaths.from_environment()
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    max_loci = int(cfg["cv_budget"]["max_loci"])
    split_id = str(cfg["split_id"])
    device = torch.device("cpu")

    print(f"[subtype-homog] loading assignment max_loci={max_loci}", flush=True)
    assignment, matrix_root = _load_assignment(paths, cfg, max_loci)
    n_region_types = max(len(assignment.region_types), 1)

    sample_index = read_sample_index(matrix_store_paths(matrix_root).sample_index_path)
    sample_ids_all = sample_index["sample_id"].astype(str).tolist()
    row_by_id = {sid: i for i, sid in enumerate(sample_ids_all)}
    n = len(sample_ids_all)

    print("[subtype-homog] loading betas into RAM…", flush=True)
    betas_all = np.asarray(open_betas_for_matrix(matrix_root)[:, :max_loci], dtype=np.float32)

    nine_pack = pd.read_parquet(
        paths.data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_hub_nine_pack_v2.parquet"
    ).set_index("sample_id")
    # Both source tables have duplicate sample_id rows (2,283 / 740 -- same
    # issue found and fixed the same way in write_nine_pack_label_status.py).
    disease_info = (
        pd.read_parquet(paths.data_root / "canonical" / "phenotypes" / "disease_sample_info.parquet")
        .drop_duplicates(subset="sample_id")
        .set_index("sample_id")
    )
    cancer_info = (
        pd.read_parquet(paths.data_root / "canonical" / "phenotypes" / "cancer_sample_info.parquet")
        .drop_duplicates(subset="sample_id")
        .set_index("sample_id")
    )

    # Alzheimer's cases vs. the broader disease-pack control pool.
    ad_y = np.zeros(n, dtype=np.int64)
    ad_mask = np.zeros(n, dtype=bool)
    for i, sid in enumerate(sample_ids_all):
        is_ad = (
            sid in disease_info.index
            and disease_info.loc[sid, "phenotype_value"] == "Alzheimer's disease"
        )
        if is_ad:
            ad_mask[i] = True
            ad_y[i] = 1
        elif sid in nine_pack.index and nine_pack.loc[sid, "disease_label_status"] == "control":
            ad_mask[i] = True
            ad_y[i] = 0
    print(
        f"[subtype-homog] AD: usable={int(ad_mask.sum())} case={int(ad_y[ad_mask].sum())}",
        flush=True,
    )

    # Cancer subtype labels (multiclass, cases only, restricted to n>=MIN_SUBTYPE_N).
    subtype_counts = cancer_info["phenotype_value"].value_counts()
    kept_subtypes = subtype_counts[subtype_counts >= MIN_SUBTYPE_N].index.tolist()
    subtype_to_idx = {s: i for i, s in enumerate(kept_subtypes)}
    subtype_y = np.full(n, -1, dtype=np.int64)
    subtype_mask = np.zeros(n, dtype=bool)
    for i, sid in enumerate(sample_ids_all):
        if sid in cancer_info.index:
            val = cancer_info.loc[sid, "phenotype_value"]
            if val in subtype_to_idx:
                subtype_mask[i] = True
                subtype_y[i] = subtype_to_idx[val]
    print(
        f"[subtype-homog] cancer subtypes kept={kept_subtypes} usable={int(subtype_mask.sum())}",
        flush=True,
    )

    # Disease / cancer case & control indices for the homogeneity diagnostic.
    disease_case = np.zeros(n, dtype=bool)
    disease_control = np.zeros(n, dtype=bool)
    cancer_case = np.zeros(n, dtype=bool)
    cancer_control = np.zeros(n, dtype=bool)
    for i, sid in enumerate(sample_ids_all):
        if sid in nine_pack.index:
            row = nine_pack.loc[sid]
            if row.get("disease_case_control_mask"):
                if row["disease_label_status"] == "case":
                    disease_case[i] = True
                elif row["disease_label_status"] == "control":
                    disease_control[i] = True
            if row.get("cancer_case_control_mask"):
                if row["cancer_label_status"] == "case":
                    cancer_case[i] = True
                elif row["cancer_label_status"] == "control":
                    cancer_control[i] = True

    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)
    folds = fold_pack["folds"]

    ad_results: list[dict[str, Any]] = []
    subtype_results: list[dict[str, Any]] = []
    homogeneity_results: dict[str, list[dict[str, Any]]] = {"disease": [], "cancer": []}

    for fold_i, fold in enumerate(folds):
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id]
        test_ids = [s for s in (fold.get("external_test_sample_ids") or []) if s in row_by_id]
        if not test_ids:
            test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id]
        train_idx = np.asarray([row_by_id[s] for s in train_ids], dtype=np.int64)
        test_idx = np.asarray([row_by_id[s] for s in test_ids], dtype=np.int64)

        ckpt_path = ROOT / "artifacts" / "runs" / CHECKPOINT_RUN_ID / f"fold_{fold_i}" / "best.pt"
        print(f"[subtype-homog] fold {fold_i}: loading {ckpt_path}", flush=True)
        model = _load_model_from_checkpoint(ckpt_path, n_region_types, device)

        mbs, present, _orphan_rbs, _all_rbs, _all_rbs_present = score_samples(
            model, assignment, betas_all, device=device, batch_size=1
        )
        mbs_c = np.where(present, mbs, 0.5)

        # AD vs. control.
        train_sel = train_idx[ad_mask[train_idx]]
        test_sel = test_idx[ad_mask[test_idx]]
        if train_sel.size >= 20 and test_sel.size >= 5 and len(np.unique(ad_y[train_sel])) > 1:
            clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.1)
            clf.fit(mbs_c[train_sel], ad_y[train_sel])
            proba = clf.predict_proba(mbs_c[test_sel])[:, 1]
            ad_results.append(
                {
                    "fold": fold_i,
                    "n_train": int(train_sel.size),
                    "n_test": int(test_sel.size),
                    "auroc": float(roc_auc_score(ad_y[test_sel], proba)),
                    "auprc": float(average_precision_score(ad_y[test_sel], proba)),
                }
            )
            print(f"[subtype-homog] fold {fold_i} AD result: {ad_results[-1]}", flush=True)

        # Cancer subtype multiclass (cases only).
        train_sel = train_idx[subtype_mask[train_idx]]
        test_sel = test_idx[subtype_mask[test_idx]]
        if train_sel.size >= 50 and test_sel.size >= 20:
            # Newer sklearn dropped multi_class= (LogisticRegression now always
            # fits a genuine multinomial model for >2 classes by default).
            clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.1)
            clf.fit(mbs_c[train_sel], subtype_y[train_sel])
            pred = clf.predict(mbs_c[test_sel])
            subtype_results.append(
                {
                    "fold": fold_i,
                    "n_train": int(train_sel.size),
                    "n_test": int(test_sel.size),
                    "macro_f1": float(f1_score(subtype_y[test_sel], pred, average="macro")),
                    "n_classes": len(kept_subtypes),
                }
            )
            print(
                f"[subtype-homog] fold {fold_i} subtype result: {subtype_results[-1]}",
                flush=True,
            )

        # Homogeneity diagnostic (MBS features; whole-cohort, not fold-split --
        # this is a descriptive statistic, not a held-out predictive claim).
        d_case_idx = np.where(disease_case)[0]
        d_control_idx = np.where(disease_control)[0]
        c_case_idx = np.where(cancer_case)[0]
        c_control_idx = np.where(cancer_control)[0]
        homogeneity_results["disease"].append(
            {"fold": fold_i, **_homogeneity(mbs_c, d_case_idx, d_control_idx)}
        )
        homogeneity_results["cancer"].append(
            {"fold": fold_i, **_homogeneity(mbs_c, c_case_idx, c_control_idx)}
        )
        print(
            f"[subtype-homog] fold {fold_i} disease homogeneity: {homogeneity_results['disease'][-1]}",
            flush=True,
        )
        print(
            f"[subtype-homog] fold {fold_i} cancer homogeneity: {homogeneity_results['cancer'][-1]}",
            flush=True,
        )

    out = {
        "generated": datetime.now(UTC).isoformat(),
        "checkpoint_run_id": CHECKPOINT_RUN_ID,
        "split_id": split_id,
        "alzheimers_vs_control": ad_results,
        "cancer_subtypes_kept": kept_subtypes,
        "cancer_subtype_multiclass": subtype_results,
        "homogeneity": homogeneity_results,
    }
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = REPORT_DIR / "subtype_and_homogeneity.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"[subtype-homog] wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
