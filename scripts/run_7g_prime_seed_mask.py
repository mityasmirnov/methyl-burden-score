#!/usr/bin/env python3
"""Run 7G′ age-primary seed-mask screen (G0/G1/G2/G3/C0/C2).

Stage B CpG-panel GPU is blocked until this screen and typed-RBS diagnostics
land. First seed source: ``internal_fold`` only (ADR 0011).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from mbs.annotation.manifest import write_json
from mbs.matrix.store import matrix_store_paths, open_betas_zarr, read_locus_index, read_sample_index
from mbs.paths import DataPaths
from mbs.training.cascade_assign import (
    CascadeAssignment,
    assignment_col_subset,
    assignment_gene_linked_only,
    build_cascade_assignment,
    gene_linked_col_index,
)
from mbs.training.cascade_loop import (
    AGE_PRIMARY_SELECTION,
    train_cascade_on_arrays,
)
from mbs.training.classical_mvalue import EPSILON, fit_eval_mvalue_fold
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.features import beta_to_m_value
from mbs.training.loop import load_experiment_config
from mbs.training.locus_gene import load_graph_tables
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.seed_panel import (
    build_internal_fold_seed_panel,
    gene_mask_tensor,
    matched_random_gene_panel,
    write_seed_panel,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_prime_seed_mask.yaml"


def _phenotype_arrays(phenotypes: list[Any], sample_ids: list[str]) -> dict[str, np.ndarray]:
    by_id = {p.sample_id: p for p in phenotypes}
    return {
        "tissue": np.asarray([int(by_id[s].class_index) for s in sample_ids], dtype=np.int64),
        "tissue_mask": np.asarray([bool(by_id[s].tissue_mask) for s in sample_ids], dtype=bool),
        "age": np.asarray([float(by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64),
        "age_mask": np.asarray([bool(by_id[s].age_mask) for s in sample_ids], dtype=bool),
        "sex": np.asarray([int(by_id[s].sex_class_index or 0) for s in sample_ids], dtype=np.int64),
        "sex_mask": np.asarray([bool(by_id[s].sex_mask) for s in sample_ids], dtype=bool),
    }


def _remap_assignment_cols(assignment: CascadeAssignment, col_map: dict[int, int]) -> CascadeAssignment:
    """Rewrite edge/direct column indices into a dense local column space."""
    edge_col = np.asarray(
        [col_map[int(c)] for c in assignment.edge_col_index.tolist() if int(c) in col_map],
        dtype=np.int64,
    )
    edge_reg = np.asarray(
        [
            int(r)
            for c, r in zip(
                assignment.edge_col_index.tolist(),
                assignment.edge_region_index.tolist(),
                strict=True,
            )
            if int(c) in col_map
        ],
        dtype=np.int64,
    )
    direct = np.asarray(
        [col_map[int(c)] for c in assignment.direct_col_index.tolist() if int(c) in col_map],
        dtype=np.int64,
    )
    return CascadeAssignment(
        gene_ids=list(assignment.gene_ids),
        region_ids=list(assignment.region_ids),
        region_type_id=np.asarray(assignment.region_type_id, dtype=np.int64),
        region_to_gene=np.asarray(assignment.region_to_gene, dtype=np.int64),
        orphan_region_mask=np.asarray(assignment.orphan_region_mask, dtype=bool),
        edge_col_index=edge_col,
        edge_region_index=edge_reg,
        direct_col_index=direct,
        region_types=assignment.region_types,
        n_study_loci=len(col_map),
        allocated_gene_id=list(assignment.allocated_gene_id),
    )


def _gene_indices_for_trait(artifacts: Any, trait: str, gene_ids: list[str]) -> list[int]:
    genes = artifacts.genes
    trait_genes = genes.loc[genes["trait"].astype(str) == trait, "gene_id"].astype(str).tolist()
    id_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    return [id_to_idx[g] for g in trait_genes if g in id_to_idx]


def _masks_from_indices(
    age_idx: list[int],
    tissue_idx: list[int],
    sex_idx: list[int],
    *,
    n_genes: int,
    n_tissue: int,
) -> dict[str, torch.Tensor]:
    return {
        "age_seed_mask": torch.as_tensor(
            gene_mask_tensor(age_idx, n_genes, n_outputs=1), dtype=torch.float32
        ),
        "tissue_seed_mask": torch.as_tensor(
            gene_mask_tensor(tissue_idx, n_genes, n_outputs=n_tissue), dtype=torch.float32
        ),
        "sex_seed_mask": torch.as_tensor(
            gene_mask_tensor(sex_idx, n_genes, n_outputs=2), dtype=torch.float32
        ),
    }


def _ones_masks(n_genes: int, n_tissue: int) -> dict[str, torch.Tensor]:
    return {
        "age_seed_mask": torch.ones(1, n_genes, dtype=torch.float32),
        "tissue_seed_mask": torch.ones(n_tissue, n_genes, dtype=torch.float32),
        "sex_seed_mask": torch.ones(2, n_genes, dtype=torch.float32),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--arm", action="append", default=[], help="Subset of arm ids")
    parser.add_argument("--dry-run-panels", action="store_true", help="Build panels only")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_rel = Path(str(cfg.get("report_dir", "reports/inspection/stage0_7g_prime_seed_mask")))
    report_dir = report_rel if report_rel.is_absolute() else paths.project_root / report_rel
    report_dir.mkdir(parents=True, exist_ok=True)
    panel_root = report_dir / "seed_panels"
    panel_root.mkdir(parents=True, exist_ok=True)

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    cv = cfg.get("cv_budget") or {}
    max_loci = int(cv.get("max_loci", 65536))
    fold_ids = [int(x) for x in (cv.get("folds") or [0])]
    seeds = [int(x) for x in (cv.get("seeds") or [42, 43])]
    n_seed_genes = int((cfg.get("seed_panel") or {}).get("n_genes", cv.get("n_seed_genes", 256)))
    min_genes = int((cfg.get("seed_panel") or {}).get("min_genes", 32))
    pilot = cfg.get("pilot") or {}
    matrix_id = str(pilot.get("matrix_id", "matrix-hub-age-tissue-sex-full-v1"))
    graph_id = str(pilot.get("graph_id", "graph-grch38-gencode38-cgi-tile-v2"))
    train_cfg = cfg.get("training") or {}
    arm_filter = {str(a) for a in args.arm} if args.arm else None
    arms_cfg = cfg.get("arms") or []
    ckpt_mode = str(train_cfg.get("checkpoint_selection", AGE_PRIMARY_SELECTION))

    pheno_rel = Path(str(cfg.get("sample_phenotype_table")))
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}
    n_tissue = max(len(class_names or []), 2)

    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    matrix_paths = matrix_store_paths(paths.data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    graph_root = paths.data_root / "canonical" / "graphs" / graph_id
    lr_edges, regions = load_graph_tables(graph_root)
    genes_path = graph_root / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()

    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation="explicit_only",
    )
    gene_linked = assignment_gene_linked_only(assignment)
    gene_cols = gene_linked_col_index(gene_linked)
    col_map = {int(c): i for i, c in enumerate(gene_cols.tolist())}
    local_assignment = _remap_assignment_cols(gene_linked, col_map)
    gene_ids = list(local_assignment.gene_ids)
    n_genes = len(gene_ids)
    n_cols = int(assignment.n_study_loci)

    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }
    print(f"[seed-mask] loading betas[:, :{n_cols}] …", flush=True)
    betas_all = np.asarray(open_betas_zarr(matrix_paths.betas_path)[:, :n_cols], dtype=np.float32)
    m_all = np.asarray(
        beta_to_m_value(np.clip(betas_all, 0, 1), epsilon=EPSILON),
        dtype=np.float32,
    )
    locus_chrom = (
        locus_index["chrom"].astype(str).to_numpy()
        if "chrom" in locus_index.columns
        else None
    )

    summary: dict[str, Any] = {
        "milestone": "7G-prime-seed-mask",
        "selection_source": "internal_fold",
        "n_seed_genes": n_seed_genes,
        "folds": fold_ids,
        "seeds": seeds,
        "checkpoint_selection": ckpt_mode,
        "arms": {},
    }

    for fold_idx in fold_ids:
        fold = fold_pack["folds"][fold_idx]
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id and s in ph_by_id]
        external = fold.get("external_test_sample_ids") or []
        test_ids = [s for s in external if s in row_by_id and s in ph_by_id]
        if not test_ids:
            test_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id and s in ph_by_id]
        val_ids = [s for s in fold["validation_sample_ids"] if s in row_by_id and s in ph_by_id]
        ordered: list[str] = []
        seen: set[str] = set()
        for s in train_ids + val_ids + test_ids:
            if s not in seen:
                seen.add(s)
                ordered.append(s)
        sample_ids = ordered
        id_to_pos = {s: i for i, s in enumerate(sample_ids)}
        rows = np.asarray([row_by_id[s] for s in sample_ids], dtype=np.int64)
        betas = betas_all[rows]
        m_vals = m_all[rows]
        train_idx = np.asarray([id_to_pos[s] for s in train_ids], dtype=np.int64)
        val_idx = np.asarray([id_to_pos[s] for s in val_ids], dtype=np.int64)
        test_idx = np.asarray([id_to_pos[s] for s in test_ids], dtype=np.int64)
        ph = _phenotype_arrays(phenotypes, sample_ids)
        studies = np.asarray([str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object)
        excluded = sorted({str(ph_by_id[s].study_id or "NA") for s in test_ids})

        print(f"[seed-mask] fold={fold_idx} building internal_fold panel K={n_seed_genes}", flush=True)
        artifacts = build_internal_fold_seed_panel(
            x_train=m_vals[train_idx][:, gene_cols],
            age=ph["age"][train_idx],
            age_mask=ph["age_mask"][train_idx],
            sex=ph["sex"][train_idx],
            sex_mask=ph["sex_mask"][train_idx],
            tissue=ph["tissue"][train_idx],
            tissue_mask=ph["tissue_mask"][train_idx],
            study_ids=studies[train_idx],
            assignment=local_assignment,
            locus_chrom=None if locus_chrom is None else locus_chrom[gene_cols],
            n_genes=n_seed_genes,
            fold_id=fold_idx,
            excluded_study_ids=excluded,
            graph_id=graph_id,
            matrix_id=matrix_id,
            seed=42 + fold_idx,
            min_genes=min_genes,
        )
        panel_dir = panel_root / f"fold_{fold_idx}"
        write_seed_panel(panel_dir, artifacts)

        age_idx = _gene_indices_for_trait(artifacts, "age", gene_ids)
        tissue_idx = _gene_indices_for_trait(artifacts, "tissue", gene_ids)
        sex_idx = _gene_indices_for_trait(artifacts, "sex", gene_ids)
        seed_masks = _masks_from_indices(
            age_idx, tissue_idx, sex_idx, n_genes=n_genes, n_tissue=n_tissue
        )
        # locus_col is local to gene_cols; map back to matrix columns.
        local_cols = sorted({int(c) for c in artifacts.loci["locus_col"].astype(int).tolist()})
        seed_cols = gene_cols[np.asarray(local_cols, dtype=np.int64)]

        seed_gene_set = sorted(set(artifacts.genes["gene_id"].astype(str).tolist()))
        gene_cpg_counts: dict[str, int] = {gid: 0 for gid in gene_ids}
        for col, reg in zip(
            local_assignment.edge_col_index.tolist(),
            local_assignment.edge_region_index.tolist(),
            strict=True,
        ):
            g = int(local_assignment.region_to_gene[int(reg)])
            if 0 <= g < n_genes:
                gene_cpg_counts[gene_ids[g]] += 1
        rng = np.random.default_rng(123 + fold_idx)
        random_genes = matched_random_gene_panel(
            seed_gene_set,
            candidate_gene_ids=gene_ids,
            gene_cpg_counts=gene_cpg_counts,
            gene_length_bp=None,
            gene_role_coverage=None,
            rng=rng,
        )
        random_idx = [gene_ids.index(g) for g in random_genes if g in gene_ids]
        random_masks = _masks_from_indices(
            random_idx, random_idx, random_idx, n_genes=n_genes, n_tissue=n_tissue
        )
        random_set = set(random_idx)
        random_local = sorted(
            {
                int(c)
                for c, reg in zip(
                    local_assignment.edge_col_index.tolist(),
                    local_assignment.edge_region_index.tolist(),
                    strict=True,
                )
                if int(local_assignment.region_to_gene[int(reg)]) in random_set
            }
        )
        random_cols = (
            gene_cols[np.asarray(random_local, dtype=np.int64)]
            if random_local
            else seed_cols
        )

        if args.dry_run_panels:
            print(f"[seed-mask] dry-run panels written under {panel_dir}", flush=True)
            continue

        for arm in arms_cfg:
            arm_id = str(arm.get("id"))
            if arm_filter is not None and arm_id not in arm_filter:
                continue
            kind = str(arm.get("kind"))
            summary["arms"].setdefault(arm_id, {"folds": []})
            for seed in seeds:
                run_id = f"stage0-7g-seed-{arm_id}-f{fold_idx}-s{seed}"
                out_dir = paths.artifact_root / "runs" / run_id
                print(f"[seed-mask] arm={arm_id} seed={seed} kind={kind}", flush=True)

                if kind.startswith("classical"):
                    cols = gene_cols if kind == "classical_enet_all" else seed_cols
                    metrics = fit_eval_mvalue_fold(
                        m_vals[train_idx][:, cols],
                        m_vals[test_idx][:, cols],
                        {k: ph[k][train_idx] for k in ph},
                        {k: ph[k][test_idx] for k in ph},
                        "enet",
                        impute=True,
                    )
                    out_dir.mkdir(parents=True, exist_ok=True)
                    write_json(out_dir / "metrics.json", {"arm_id": arm_id, "metrics": metrics})
                    summary["arms"][arm_id]["folds"].append(
                        {"fold": fold_idx, "seed": seed, "metrics": metrics}
                    )
                    continue

                if kind == "cascade_all_genes":
                    assign = gene_linked
                    masks = _ones_masks(n_genes, n_tissue)
                elif kind == "cascade_all_cpgs_seed_heads":
                    assign = gene_linked
                    masks = seed_masks
                elif kind == "cascade_seed_cpgs_seed_heads":
                    # Keep original matrix column indices; pass full betas.
                    assign = assignment_col_subset(assignment, seed_cols)
                    masks = seed_masks
                elif kind == "cascade_matched_random":
                    assign = assignment_col_subset(assignment, random_cols)
                    masks = random_masks
                else:
                    raise ValueError(f"unknown arm kind: {kind}")

                metrics = train_cascade_on_arrays(
                    assignment=assign,
                    betas=betas,
                    train_idx=train_idx,
                    test_idx=test_idx,
                    ages=ph["age"],
                    tissue=ph["tissue"],
                    sex=ph["sex"],
                    study_ids=studies,
                    sample_ids=sample_ids,
                    class_names=class_names or ["A", "B"],
                    out_dir=out_dir,
                    max_epochs=int(train_cfg.get("max_epochs", 15)),
                    seed=seed,
                    device_str=args.device,
                    age_mask=ph["age_mask"],
                    tissue_mask=ph["tissue_mask"],
                    sex_mask=ph["sex_mask"],
                    cpg_pool=str((cfg.get("model") or {}).get("cpg_pool", "max")),
                    region_pool=str((cfg.get("model") or {}).get("region_pool", "max")),
                    gene_allocation_policy="explicit_only",
                    val_idx=val_idx,
                    age_loss_weight=float(train_cfg.get("age_loss_weight", 1.0)),
                    tissue_loss_weight=float(train_cfg.get("tissue_loss_weight", 0.3)),
                    sex_loss_weight=float(train_cfg.get("sex_loss_weight", 0.1)),
                    primary_evaluation="mbs_e2e",
                    gene_linked_only=False,
                    include_mbs_enet=False,
                    age_seed_mask=masks["age_seed_mask"],
                    tissue_seed_mask=masks["tissue_seed_mask"],
                    sex_seed_mask=masks["sex_seed_mask"],
                    checkpoint_selection_mode=ckpt_mode,
                )
                summary["arms"][arm_id]["folds"].append(
                    {"fold": fold_idx, "seed": seed, "run_id": run_id, "metrics": metrics}
                )

    write_json(report_dir / "summary.json", summary)
    analysis = [
        "# 7G′ age-primary seed-mask screen",
        "",
        "Selection: validation age MAE primary; tissue F1 secondary; sex AUROC tertiary.",
        "Seed source: **internal_fold** (ADR 0011). P2-G topology is a reference, not a lock.",
        "",
        f"Folds: {fold_ids}; seeds: {seeds}; K={n_seed_genes}.",
        "",
        "Arms: G0 all-gene control; G1 head masks; G2 seed CpGs+masks; G3 matched random;",
        "C0 classical all-gene; C2 classical on G2 CpGs.",
        "",
        f"See `{report_dir / 'summary.json'}`.",
        "",
    ]
    (report_dir / "analysis.md").write_text("\n".join(analysis), encoding="utf-8")
    print(f"[seed-mask] wrote {report_dir / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
