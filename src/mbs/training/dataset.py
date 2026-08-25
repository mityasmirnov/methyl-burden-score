"""Flat sample records and batch tensors for FlatDeepSet training."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import torch
from torch import Tensor

from mbs.training.features import (
    SampleFeatureBundle,
    assemble_cpg_features,
    beta_to_m_value,
    cpg_input_dim,
    gather_sample_features,
)
from mbs.training.level1_norm import (
    Level1NormParams,
    apply_level1_robust_z,
    fit_level1_robust_z,
)
from mbs.training.locus_gene import LocusGeneIndex
from mbs.training.phenotypes import SamplePhenotype


@dataclass(frozen=True, slots=True)
class FlatSampleRecord:
    sample_id: str
    donor_id: str | None
    class_index: int
    features: SampleFeatureBundle


@dataclass(slots=True)
class FlatBatch:
    """Single-sample (or small) flat batch ready for ``FlatDeepSet``."""

    sample_ids: list[str]
    cpg_features: Tensor
    cpg_to_gene: Tensor
    n_genes: int
    tissue_target: Tensor
    tissue_mask: Tensor
    age_target: Tensor | None
    age_mask: Tensor
    sex_target: Tensor | None = None
    sex_mask: Tensor | None = None
    disease_target: Tensor | None = None
    disease_mask: Tensor | None = None
    cancer_target: Tensor | None = None
    cancer_mask: Tensor | None = None

    def to(self, device: torch.device | str, *, non_blocking: bool = False) -> FlatBatch:
        sex_mask = self.sex_mask
        if sex_mask is None:
            sex_mask = torch.zeros(len(self.sample_ids), dtype=torch.bool)
        sex_target = self.sex_target
        if sex_target is None:
            sex_target = torch.zeros(len(self.sample_ids), dtype=torch.long)
        return FlatBatch(
            sample_ids=list(self.sample_ids),
            cpg_features=self.cpg_features.to(device, non_blocking=non_blocking),
            cpg_to_gene=self.cpg_to_gene.to(device, non_blocking=non_blocking),
            n_genes=self.n_genes,
            tissue_target=self.tissue_target.to(device, non_blocking=non_blocking),
            tissue_mask=self.tissue_mask.to(device, non_blocking=non_blocking),
            age_target=(
                None
                if self.age_target is None
                else self.age_target.to(device, non_blocking=non_blocking)
            ),
            age_mask=self.age_mask.to(device, non_blocking=non_blocking),
            sex_target=sex_target.to(device, non_blocking=non_blocking),
            sex_mask=sex_mask.to(device, non_blocking=non_blocking),
            disease_target=(
                None
                if self.disease_target is None
                else self.disease_target.to(device, non_blocking=non_blocking)
            ),
            disease_mask=(
                None
                if self.disease_mask is None
                else self.disease_mask.to(device, non_blocking=non_blocking)
            ),
            cancer_target=(
                None
                if self.cancer_target is None
                else self.cancer_target.to(device, non_blocking=non_blocking)
            ),
            cancer_mask=(
                None
                if self.cancer_mask is None
                else self.cancer_mask.to(device, non_blocking=non_blocking)
            ),
        )


def build_flat_sample(
    *,
    phenotype: SamplePhenotype,
    beta_row: np.ndarray,
    static_by_col: np.ndarray,
    static_valid: np.ndarray,
    locus_gene: LocusGeneIndex,
    epsilon: float = 0.001,
    include_m_value: bool = True,
    include_robust_z: bool = False,
    level1_params: Any | None = None,
) -> FlatSampleRecord:
    features = gather_sample_features(
        beta_row=beta_row,
        static_by_col=static_by_col,
        static_valid=static_valid,
        locus_gene=locus_gene,
        epsilon=epsilon,
        include_m_value=include_m_value,
        include_robust_z=include_robust_z,
        level1_params=level1_params,
    )
    if features.n_observed_edges == 0:
        raise ValueError(f"sample {phenotype.sample_id!r} has zero observed gene-mapped CpGs")
    return FlatSampleRecord(
        sample_id=phenotype.sample_id,
        donor_id=phenotype.donor_id,
        class_index=phenotype.class_index,
        features=features,
    )


def record_to_batch(
    record: FlatSampleRecord,
    *,
    n_genes: int,
    age_value: float | None = None,
    age_enabled: bool = False,
    tissue_enabled: bool = True,
    sex_enabled: bool = False,
    sex_class_index: int = 0,
    disease_target: np.ndarray | None = None,
    disease_mask: np.ndarray | None = None,
    cancer_target: np.ndarray | None = None,
    cancer_mask: np.ndarray | None = None,
) -> FlatBatch:
    feats = record.features
    tissue_on = bool(tissue_enabled)
    age_on = bool(age_enabled and age_value is not None)
    sex_on = bool(sex_enabled)
    age_target: Tensor | None = None
    if age_on:
        if age_value is None:
            raise RuntimeError("age_enabled set but age_value is None")
        age_target = torch.tensor([float(age_value)], dtype=torch.float32)
    dis_t = (
        None
        if disease_target is None
        else torch.as_tensor(disease_target, dtype=torch.float32).unsqueeze(0)
    )
    dis_m = (
        None
        if disease_mask is None
        else torch.as_tensor(disease_mask, dtype=torch.bool).unsqueeze(0)
    )
    can_t = (
        None
        if cancer_target is None
        else torch.as_tensor(cancer_target, dtype=torch.float32).unsqueeze(0)
    )
    can_m = (
        None if cancer_mask is None else torch.as_tensor(cancer_mask, dtype=torch.bool).unsqueeze(0)
    )
    return FlatBatch(
        sample_ids=[record.sample_id],
        cpg_features=torch.from_numpy(feats.cpg_features),
        cpg_to_gene=torch.from_numpy(feats.cpg_to_gene),
        n_genes=n_genes,
        tissue_target=torch.tensor([record.class_index], dtype=torch.long),
        tissue_mask=torch.tensor([tissue_on]),
        age_target=age_target,
        age_mask=torch.tensor([age_on]),
        sex_target=torch.tensor([int(sex_class_index)], dtype=torch.long),
        sex_mask=torch.tensor([sex_on]),
        disease_target=dis_t,
        disease_mask=dis_m,
        cancer_target=can_t,
        cancer_mask=can_m,
    )


def pack_records_to_batch(
    records: list[FlatSampleRecord],
    *,
    n_genes: int,
    age_values: list[float | None],
    age_enabled: list[bool],
    tissue_enabled: list[bool],
    sex_enabled: list[bool] | None = None,
    sex_class_indices: list[int] | None = None,
    disease_targets: list[np.ndarray | None] | None = None,
    disease_masks: list[np.ndarray | None] | None = None,
    cancer_targets: list[np.ndarray | None] | None = None,
    cancer_masks: list[np.ndarray | None] | None = None,
) -> FlatBatch:
    """Pack ragged samples into one FlatDeepSet forward via gene-index offsets.

    Gene ids for sample ``i`` become ``gene + i * n_genes`` so a single
    ``segment_pool`` over ``B * n_genes`` segments yields ``[B, n_genes]`` MBS.
    """
    if not records:
        raise ValueError("pack_records_to_batch requires at least one record")
    if not (len(records) == len(age_values) == len(age_enabled) == len(tissue_enabled)):
        raise ValueError("records and label lists must have equal length")
    sex_flags = sex_enabled if sex_enabled is not None else [False] * len(records)
    sex_idxs = sex_class_indices if sex_class_indices is not None else [0] * len(records)
    if len(sex_flags) != len(records) or len(sex_idxs) != len(records):
        raise ValueError("sex label lists must match records length")

    feat_parts: list[Tensor] = []
    gene_parts: list[Tensor] = []
    tissue_targets: list[int] = []
    tissue_masks: list[bool] = []
    age_targets: list[float] = []
    age_masks: list[bool] = []
    sex_targets: list[int] = []
    sex_masks: list[bool] = []
    sample_ids: list[str] = []

    for i, record in enumerate(records):
        feats = record.features
        feat_parts.append(torch.from_numpy(feats.cpg_features))
        gene_parts.append(torch.from_numpy(feats.cpg_to_gene).to(torch.long) + i * int(n_genes))
        tissue_on = bool(tissue_enabled[i])
        age_val = age_values[i]
        age_on = bool(age_enabled[i] and age_val is not None)
        sex_on = bool(sex_flags[i])
        tissue_targets.append(int(record.class_index))
        tissue_masks.append(tissue_on)
        age_masks.append(age_on)
        age_targets.append(0.0 if age_val is None else float(age_val))
        sex_masks.append(sex_on)
        sex_targets.append(int(sex_idxs[i]))
        sample_ids.append(record.sample_id)

    def _stack_optional(
        rows: list[np.ndarray | None] | None,
        *,
        dtype: torch.dtype,
    ) -> Tensor | None:
        if rows is None or all(r is None for r in rows):
            return None
        if any(r is None for r in rows):
            raise ValueError("optional multilabel rows must be all-or-none for a batch")
        stacked = np.stack([np.asarray(r) for r in rows if r is not None], axis=0)
        return torch.as_tensor(stacked, dtype=dtype)

    return FlatBatch(
        sample_ids=sample_ids,
        cpg_features=torch.cat(feat_parts, dim=0),
        cpg_to_gene=torch.cat(gene_parts, dim=0),
        n_genes=n_genes,
        tissue_target=torch.tensor(tissue_targets, dtype=torch.long),
        tissue_mask=torch.tensor(tissue_masks, dtype=torch.bool),
        age_target=torch.tensor(age_targets, dtype=torch.float32),
        age_mask=torch.tensor(age_masks, dtype=torch.bool),
        sex_target=torch.tensor(sex_targets, dtype=torch.long),
        sex_mask=torch.tensor(sex_masks, dtype=torch.bool),
        disease_target=_stack_optional(disease_targets, dtype=torch.float32),
        disease_mask=_stack_optional(disease_masks, dtype=torch.bool),
        cancer_target=_stack_optional(cancer_targets, dtype=torch.float32),
        cancer_mask=_stack_optional(cancer_masks, dtype=torch.bool),
    )


def make_synthetic_overfit_bundle(
    *,
    n_samples: int = 12,
    n_cpgs: int = 15,
    n_genes: int = 3,
    n_classes: int = 3,
    static_dim: int = 4,
    seed: int = 0,
    include_m_value: bool = True,
    include_robust_z: bool = False,
) -> dict[str, Any]:
    """Tiny in-memory fixture that is linearly separable by class-tied patterns.

    Each class lights up exactly one gene (high beta on that gene's CpGs).
    When ``include_robust_z``, Level-1 params are fit on all fixture samples
    (overfit path); study-holdout trains must re-fit on the train fold only.
    """
    if n_genes < n_classes:
        raise ValueError("n_genes must be >= n_classes for the overfit fixture")
    if include_robust_z and not include_m_value:
        raise ValueError("include_robust_z requires include_m_value")
    rng = np.random.default_rng(seed)
    # Equal CpGs per gene
    cpg_to_gene = np.array([i % n_genes for i in range(n_cpgs)], dtype=np.int64)
    gene_ids = [f"GENE{i}" for i in range(n_genes)]
    class_names = [f"class_{i}" for i in range(n_classes)]
    edge_cols = np.arange(n_cpgs, dtype=np.int64)

    beta_rows: list[np.ndarray] = []
    static_rows: list[np.ndarray] = []
    records: list[FlatSampleRecord] = []
    ages: list[float] = []
    for s in range(n_samples):
        cls = s % n_classes
        betas = np.full(n_cpgs, 0.05, dtype=np.float32)
        active = cpg_to_gene == cls
        betas[active] = 0.95
        # Tiny noise that cannot flip the active gene
        betas = np.clip(betas + rng.normal(0, 0.005, size=n_cpgs), 0.01, 0.99).astype(np.float32)
        static = rng.normal(0, 0.01, size=(n_cpgs, static_dim)).astype(np.float32)
        static[active, 0] = 1.0 + float(cls)
        beta_rows.append(betas)
        static_rows.append(static)
        ages.append(20.0 + 10.0 * cls)

    level1 = None
    if include_robust_z:
        m_mat = np.stack([beta_to_m_value(b) for b in beta_rows], axis=0)
        level1 = fit_level1_robust_z(m_mat, sigma_min=1e-6)

    for s in range(n_samples):
        betas = beta_rows[s]
        static = static_rows[s]
        robust_z = None
        norm_present = None
        if include_robust_z and level1 is not None:
            robust_z, norm_present = apply_level1_robust_z(
                beta_to_m_value(betas), level1, col_index=edge_cols
            )
        features = assemble_cpg_features(
            betas=betas,
            static_rows=static,
            static_present=np.ones(n_cpgs, dtype=np.float32),
            include_m_value=include_m_value,
            include_robust_z=include_robust_z,
            robust_z=robust_z,
            norm_present=norm_present,
        )
        bundle = SampleFeatureBundle(
            cpg_features=features,
            cpg_to_gene=cpg_to_gene.copy(),
            n_observed_edges=n_cpgs,
            n_dropped_nan_beta=0,
            n_dropped_no_static=0,
            edge_col_index=edge_cols.copy(),
        )
        records.append(
            FlatSampleRecord(
                sample_id=f"SYN{s:03d}",
                donor_id=str((s % 4) + 1),
                class_index=s % n_classes,
                features=bundle,
            )
        )

    return {
        "records": records,
        "gene_ids": gene_ids,
        "class_names": class_names,
        "ages": ages,
        "n_genes": n_genes,
        "n_cpgs": n_cpgs,
        "static_dim": static_dim,
        "input_dim": cpg_input_dim(
            static_dim, include_m_value=include_m_value, include_robust_z=include_robust_z
        ),
        "n_classes": n_classes,
        "level1_params": level1,
        "include_m_value": include_m_value,
        "include_robust_z": include_robust_z,
    }


def make_synthetic_study_holdout_bundle(
    *,
    studies: tuple[str, ...] = ("GSE_A", "GSE_B", "GSE_C"),
    samples_per_study: int = 6,
    n_cpgs: int = 15,
    n_genes: int = 3,
    n_classes: int = 3,
    static_dim: int = 4,
    seed: int = 0,
    task: str = "tissue",
    include_m_value: bool = True,
    include_robust_z: bool = False,
) -> dict[str, Any]:
    """Multi-study synthetic fixture for study-grouped holdout benchmarks.

    ``task='tissue'`` → multiclass labels; ``task='age'`` → ages tied to class
    with continuous targets (class_index still set for shared head smoke).

    When ``include_robust_z``, the base bundle fits Level-1 on **all** studies;
    callers that need leakage-safe params must re-fit on the train fold via
    ``refit_level1_on_flat_records``.
    """
    base = make_synthetic_overfit_bundle(
        n_samples=len(studies) * samples_per_study,
        n_cpgs=n_cpgs,
        n_genes=n_genes,
        n_classes=n_classes,
        static_dim=static_dim,
        seed=seed,
        include_m_value=include_m_value,
        include_robust_z=include_robust_z,
    )
    records: list[FlatSampleRecord] = []
    sample_rows: list[dict[str, Any]] = []
    ages: list[float] = []
    for i, rec in enumerate(base["records"]):
        study = studies[i // samples_per_study]
        new_id = f"{study}_{rec.sample_id}"
        records.append(
            FlatSampleRecord(
                sample_id=new_id,
                donor_id=f"{study}_{rec.donor_id}",
                class_index=rec.class_index,
                features=rec.features,
            )
        )
        age = float(base["ages"][i]) + (0.1 if task == "age" else 0.0)
        ages.append(age)
        sample_rows.append(
            {
                "sample_id": new_id,
                "study_id": study,
                "platform": "HM450",
                "group_id": study,
                "class_index": rec.class_index,
                "age": age,
            }
        )
    return {
        **base,
        "records": records,
        "ages": ages,
        "sample_rows": sample_rows,
        "studies": list(studies),
        "task": task,
    }


def refit_level1_on_flat_records(
    train_records: list[FlatSampleRecord],
    all_records: list[FlatSampleRecord],
    *,
    include_m_value: bool = True,
    sigma_min: float = 1e-6,
    epsilon: float = 0.001,
    fold_id: str | None = None,
    run_id: str | None = None,
) -> tuple[Any, list[FlatSampleRecord]]:
    """Fit Level-1 on train records' M columns; rebuild all records with channel B.

    Expects channel-A features (beta, M, static..., static_present) and
    ``edge_col_index`` on each bundle. Returns params and new record list.
    """
    if not include_m_value:
        raise ValueError("Level-1 refit requires m_value")
    if not train_records:
        raise ValueError("cannot fit Level-1 without train records")

    n_loci = 0
    for rec in train_records:
        cols = rec.features.edge_col_index
        if cols is None:
            raise ValueError("edge_col_index required for Level-1 refit")
        n_loci = max(n_loci, int(cols.max()) + 1 if cols.size else 0)
    for rec in all_records:
        cols = rec.features.edge_col_index
        if cols is None:
            raise ValueError("edge_col_index required for Level-1 refit")
        n_loci = max(n_loci, int(cols.max()) + 1 if cols.size else 0)

    m_mat = np.full((len(train_records), n_loci), np.nan, dtype=np.float64)
    for i, rec in enumerate(train_records):
        feats = rec.features.cpg_features
        cols = rec.features.edge_col_index
        if cols is None:
            raise ValueError("edge_col_index required for Level-1 refit")
        # Channel A: beta=0, M=1
        m_vals = feats[:, 1]
        m_mat[i, cols] = m_vals

    params: Level1NormParams = fit_level1_robust_z(
        m_mat,
        sigma_min=sigma_min,
        epsilon=epsilon,
        fold_id=fold_id,
        run_id=run_id,
    )

    rebuilt: list[FlatSampleRecord] = []
    for rec in all_records:
        feats = rec.features.cpg_features
        cols = rec.features.edge_col_index
        if cols is None:
            raise ValueError("edge_col_index required for Level-1 refit")
        betas = feats[:, 0]
        m_vals = feats[:, 1]
        # static starts after beta+M; ends before static_present
        static = feats[:, 2:-1]
        static_present = feats[:, -1]
        z, norm_present = apply_level1_robust_z(m_vals, params, col_index=cols)
        new_feats = assemble_cpg_features(
            betas=betas,
            static_rows=static,
            static_present=static_present,
            include_m_value=True,
            include_robust_z=True,
            robust_z=z,
            norm_present=norm_present,
            epsilon=epsilon,
        )
        bundle = SampleFeatureBundle(
            cpg_features=new_feats,
            cpg_to_gene=rec.features.cpg_to_gene.copy(),
            n_observed_edges=rec.features.n_observed_edges,
            n_dropped_nan_beta=rec.features.n_dropped_nan_beta,
            n_dropped_no_static=rec.features.n_dropped_no_static,
            n_missing_static=rec.features.n_missing_static,
            edge_col_index=cols.copy(),
        )
        rebuilt.append(
            FlatSampleRecord(
                sample_id=rec.sample_id,
                donor_id=rec.donor_id,
                class_index=rec.class_index,
                features=bundle,
            )
        )
    return params, rebuilt
