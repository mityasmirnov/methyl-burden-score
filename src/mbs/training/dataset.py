"""Flat sample records and batch tensors for FlatDeepSet training."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import torch
from torch import Tensor

from mbs.training.features import SampleFeatureBundle, beta_to_m_value, gather_sample_features
from mbs.training.locus_gene import LocusGeneIndex
from mbs.training.phenotypes import SamplePhenotype


@dataclass(frozen=True, slots=True)
class FlatSampleRecord:
    sample_id: str
    donor_id: str
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

    def to(self, device: torch.device | str, *, non_blocking: bool = False) -> FlatBatch:
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
) -> FlatBatch:
    feats = record.features
    return FlatBatch(
        sample_ids=[record.sample_id],
        cpg_features=torch.from_numpy(feats.cpg_features),
        cpg_to_gene=torch.from_numpy(feats.cpg_to_gene),
        n_genes=n_genes,
        tissue_target=torch.tensor([record.class_index], dtype=torch.long),
        tissue_mask=torch.tensor([True]),
        age_target=(
            torch.tensor([float(age_value)], dtype=torch.float32)
            if age_enabled and age_value is not None
            else None
        ),
        age_mask=torch.tensor([bool(age_enabled and age_value is not None)]),
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

        m_vals = beta_to_m_value(betas)
        features = np.concatenate(
            [betas.reshape(-1, 1), m_vals.reshape(-1, 1), static],
            axis=1,
        ).astype(np.float32)
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
        "input_dim": 2 + static_dim,
        "n_classes": n_classes,
    }
