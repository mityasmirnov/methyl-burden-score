"""Hub transparent arms for Milestone 7G (T-mean-gene / T-mean-region / T-enet)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import numpy as np

from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)
from mbs.training.dev_cv import _phenotype_arrays
from mbs.training.locus_gene import build_locus_gene_index, load_graph_tables
from mbs.training.locus_region_gene import build_locus_region_gene_index
from mbs.training.transparent_baselines import (
    presence_aware_means,
    run_elasticnet_baseline,
    run_mean_baseline,
)

TransparentArm = Literal["T-mean-gene", "T-mean-region", "T-enet"]


def run_hub_transparent_arm(
    *,
    data_root: Path,
    fold: dict[str, Any],
    phenotypes: list[Any],
    arm: TransparentArm,
    max_loci: int,
    matrix_id: str = "matrix-hub-age-tissue-sex-full-v1",
    graph_id: str = "graph-grch38-gencode38-cgi-tile-v2",
    seed: int = 42,
) -> dict[str, Any]:
    """Gene-mean, region-mean, or elastic-net-on-gene-means on one frozen fold."""
    del seed  # reserved for future subsample
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
    n_cols = min(int(max_loci), int(locus_index.shape[0]))

    if arm == "T-mean-region":
        index = build_locus_region_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=n_cols,
            region_systems=("gene", "rbs"),
        )
        cols = index.edge_col_index
        groups = index.edge_region_index
        n_groups = index.n_regions
        mean_kind: Literal["gene", "region"] = "region"
    else:
        locus_gene = build_locus_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=n_cols,
            region_systems=("gene",),
        )
        cols = locus_gene.edge_col_index
        groups = locus_gene.edge_gene_index
        n_groups = locus_gene.n_genes
        mean_kind = "gene"

    keep = cols < n_cols
    cols_k = cols[keep]
    groups_k = groups[keep]

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
        betas_tr, obs_tr, groups_k, n_groups=max(n_groups, 1), empty_fill=0.5
    )
    x_te = presence_aware_means(
        betas_te, obs_te, groups_k, n_groups=max(n_groups, 1), empty_fill=0.5
    )
    ph_tr = _phenotype_arrays(phenotypes, train_ids)
    ph_te = _phenotype_arrays(phenotypes, test_ids)
    class_names = sorted({str(t) for t in ph_tr["tissues"].tolist() if t}) or None

    if arm == "T-enet":
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
        out["arm"] = "T-enet"
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
            kind=mean_kind,
            tissue_class_names=class_names,
        )
        out["arm"] = arm
    out["n_train"] = len(train_ids)
    out["n_test"] = len(test_ids)
    out["max_loci"] = n_cols
    out["n_groups"] = int(n_groups)
    return out


def run_all_transparent_arms(
    *,
    data_root: Path,
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
    max_loci: int = 65536,
) -> dict[str, Any]:
    """Run T-mean-gene, T-mean-region, T-enet on every fold."""
    arms: list[TransparentArm] = ["T-mean-gene", "T-mean-region", "T-enet"]
    payload: dict[str, Any] = {"arms": arms, "folds": [], "max_loci": max_loci}
    for fold_idx, fold in enumerate(fold_pack["folds"]):
        fold_out: dict[str, Any] = {"fold": fold_idx, "arms": {}}
        for arm in arms:
            print(f"[transparent] fold {fold_idx} {arm}…", flush=True)
            fold_out["arms"][arm] = run_hub_transparent_arm(
                data_root=data_root,
                fold=fold,
                phenotypes=phenotypes,
                arm=arm,
                max_loci=max_loci,
            )
        payload["folds"].append(fold_out)
    return payload
