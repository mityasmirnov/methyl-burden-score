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
    cpg_input_dim,
    gather_sample_features,
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
) -> FlatSampleRecord:
    features = gather_sample_features(
        beta_row=beta_row,
        static_by_col=static_by_col,
        static_valid=static_valid,
        locus_gene=locus_gene,
        epsilon=epsilon,
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
    )


def make_synthetic_overfit_bundle(
    *,
    n_samples: int = 12,
    n_cpgs: int = 15,
    n_genes: int = 3,
    n_classes: int = 3,
    static_dim: int = 4,
    seed: int = 0,
) -> dict[str, Any]:
    """Tiny in-memory fixture that is linearly separable by class-tied patterns.

    Each class lights up exactly one gene (high beta on that gene's CpGs).
    """
    if n_genes < n_classes:
        raise ValueError("n_genes must be >= n_classes for the overfit fixture")
    rng = np.random.default_rng(seed)
    # Equal CpGs per gene
    cpg_to_gene = np.array([i % n_genes for i in range(n_cpgs)], dtype=np.int64)
    gene_ids = [f"GENE{i}" for i in range(n_genes)]
    class_names = [f"class_{i}" for i in range(n_classes)]

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

        features = assemble_cpg_features(
            betas=betas,
            static_rows=static,
            static_present=np.ones(n_cpgs, dtype=np.float32),
        )
        bundle = SampleFeatureBundle(
            cpg_features=features,
            cpg_to_gene=cpg_to_gene.copy(),
            n_observed_edges=n_cpgs,
            n_dropped_nan_beta=0,
            n_dropped_no_static=0,
        )
        records.append(
            FlatSampleRecord(
                sample_id=f"SYN{s:03d}",
                donor_id=str((s % 4) + 1),
                class_index=cls,
                features=bundle,
            )
        )
        ages.append(20.0 + 10.0 * cls)

    return {
        "records": records,
        "gene_ids": gene_ids,
        "class_names": class_names,
        "ages": ages,
        "n_genes": n_genes,
        "input_dim": cpg_input_dim(static_dim),
        "n_classes": n_classes,
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
) -> dict[str, Any]:
    """Multi-study synthetic fixture for study-grouped holdout benchmarks.

    ``task='tissue'`` → multiclass labels; ``task='age'`` → ages tied to class
    with continuous targets (class_index still set for shared head smoke).
    """
    base = make_synthetic_overfit_bundle(
        n_samples=len(studies) * samples_per_study,
        n_cpgs=n_cpgs,
        n_genes=n_genes,
        n_classes=n_classes,
        static_dim=static_dim,
        seed=seed,
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
