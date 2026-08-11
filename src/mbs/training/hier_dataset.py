"""Hierarchical sample records and packed batches for HierarchicalDeepSet."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import torch
from torch import Tensor

from mbs.training.features import beta_to_m_value
from mbs.training.locus_region_gene import (
    HIER_REGION_TYPES,
    REGION_TYPE_TO_ID,
    UNASSIGNED_GENE_ID,
    UNASSIGNED_REGION_TYPE,
    LocusRegionGeneIndex,
)
from mbs.training.phenotypes import SamplePhenotype


@dataclass(frozen=True, slots=True)
class HierSampleFeatureBundle:
    cpg_features: np.ndarray  # float32 [n_edges_obs, feat_dim]
    cpg_to_region: np.ndarray  # int64 [n_edges_obs]
    n_observed_edges: int
    n_dropped_nan_beta: int
    n_dropped_no_static: int


@dataclass(frozen=True, slots=True)
class HierSampleRecord:
    sample_id: str
    donor_id: str
    class_index: int
    features: HierSampleFeatureBundle


@dataclass(slots=True)
class HierBatch:
    """Packed multi-sample batch for ``HierarchicalDeepSet``."""

    sample_ids: list[str]
    cpg_features: Tensor
    cpg_to_region: Tensor
    region_type: Tensor
    region_to_gene: Tensor
    n_regions: int
    n_genes: int
    tissue_target: Tensor
    tissue_mask: Tensor
    age_target: Tensor | None
    age_mask: Tensor
    sex_target: Tensor | None = None
    sex_mask: Tensor | None = None

    def to(self, device: torch.device | str, *, non_blocking: bool = False) -> HierBatch:
        sex_mask = self.sex_mask
        if sex_mask is None:
            sex_mask = torch.zeros(len(self.sample_ids), dtype=torch.bool)
        sex_target = self.sex_target
        if sex_target is None:
            sex_target = torch.zeros(len(self.sample_ids), dtype=torch.long)
        return HierBatch(
            sample_ids=list(self.sample_ids),
            cpg_features=self.cpg_features.to(device, non_blocking=non_blocking),
            cpg_to_region=self.cpg_to_region.to(device, non_blocking=non_blocking),
            region_type=self.region_type.to(device, non_blocking=non_blocking),
            region_to_gene=self.region_to_gene.to(device, non_blocking=non_blocking),
            n_regions=self.n_regions,
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
        )


def gather_hier_sample_features(
    *,
    beta_row: np.ndarray,
    static_by_col: np.ndarray,
    static_valid: np.ndarray,
    locus_region: LocusRegionGeneIndex,
    epsilon: float = 0.001,
    include_m_value: bool = True,
) -> HierSampleFeatureBundle:
    """Build ragged hierarchical features for one sample."""
    beta_row = np.asarray(beta_row, dtype=np.float32)
    if beta_row.ndim != 1:
        raise ValueError("beta_row must be 1-D")
    if beta_row.shape[0] < locus_region.n_study_loci:
        raise ValueError(
            f"beta_row length {beta_row.shape[0]} < n_study_loci {locus_region.n_study_loci}"
        )

    cols = locus_region.edge_col_index
    regions = locus_region.edge_region_index
    betas = beta_row[cols]
    finite = np.isfinite(betas)
    static_ok = static_valid[cols]
    keep = finite & static_ok
    n_dropped_nan = int((~finite).sum())
    n_dropped_static = int((finite & ~static_ok).sum())

    cols_k = cols[keep]
    regions_k = regions[keep]
    betas_k = betas[keep]
    static_k = static_by_col[cols_k]

    parts: list[np.ndarray] = [betas_k.reshape(-1, 1)]
    if include_m_value:
        m_vals = beta_to_m_value(betas_k, epsilon=epsilon).reshape(-1, 1)
        parts.append(m_vals)
    parts.append(static_k.astype(np.float32, copy=False))
    features = np.concatenate(parts, axis=1).astype(np.float32, copy=False)

    return HierSampleFeatureBundle(
        cpg_features=features,
        cpg_to_region=regions_k.astype(np.int64, copy=False),
        n_observed_edges=int(features.shape[0]),
        n_dropped_nan_beta=n_dropped_nan,
        n_dropped_no_static=n_dropped_static,
    )


def build_hier_sample(
    *,
    phenotype: SamplePhenotype,
    beta_row: np.ndarray,
    static_by_col: np.ndarray,
    static_valid: np.ndarray,
    locus_region: LocusRegionGeneIndex,
    epsilon: float = 0.001,
) -> HierSampleRecord:
    features = gather_hier_sample_features(
        beta_row=beta_row,
        static_by_col=static_by_col,
        static_valid=static_valid,
        locus_region=locus_region,
        epsilon=epsilon,
    )
    if features.n_observed_edges == 0:
        raise ValueError(f"sample {phenotype.sample_id!r} has zero observed region-mapped CpGs")
    return HierSampleRecord(
        sample_id=phenotype.sample_id,
        donor_id=phenotype.donor_id,
        class_index=phenotype.class_index,
        features=features,
    )


def pack_hier_records_to_batch(
    records: list[HierSampleRecord],
    *,
    locus_region: LocusRegionGeneIndex,
    age_values: list[float | None],
    age_enabled: list[bool],
    tissue_enabled: list[bool],
    sex_enabled: list[bool] | None = None,
    sex_class_indices: list[int] | None = None,
    allowed_region_type_ids: set[int] | None = None,
) -> HierBatch:
    """Pack ragged samples with region/gene offsets for one hierarchical forward.

    Optional ``allowed_region_type_ids`` zeros features for disallowed region
    types (eval ablations). Topology (region_type / region_to_gene) stays fixed.
    """
    if not records:
        raise ValueError("pack_hier_records_to_batch requires at least one record")
    if not (len(records) == len(age_values) == len(age_enabled) == len(tissue_enabled)):
        raise ValueError("records and label lists must have equal length")
    sex_flags = sex_enabled if sex_enabled is not None else [False] * len(records)
    sex_idxs = sex_class_indices if sex_class_indices is not None else [0] * len(records)
    if len(sex_flags) != len(records) or len(sex_idxs) != len(records):
        raise ValueError("sex label lists must match records length")

    n_regions = locus_region.n_regions
    n_genes = locus_region.n_genes
    region_type = torch.from_numpy(locus_region.region_type_id.astype(np.int64, copy=False))
    region_to_gene = torch.from_numpy(locus_region.region_to_gene.astype(np.int64, copy=False))

    feat_parts: list[Tensor] = []
    region_parts: list[Tensor] = []
    type_parts: list[Tensor] = []
    gene_parts: list[Tensor] = []
    tissue_targets: list[int] = []
    tissue_masks: list[bool] = []
    age_targets: list[float] = []
    age_masks: list[bool] = []
    sex_targets: list[int] = []
    sex_masks: list[bool] = []
    sample_ids: list[str] = []

    allow = allowed_region_type_ids
    for i, record in enumerate(records):
        feats = record.features
        cpg_feat = torch.from_numpy(feats.cpg_features.copy())
        cpg_to_region = torch.from_numpy(feats.cpg_to_region.astype(np.int64, copy=False))
        if allow is not None:
            edge_types = region_type[cpg_to_region]
            keep = torch.tensor(
                [int(t) in allow for t in edge_types.tolist()],
                dtype=torch.bool,
            )
            cpg_feat = cpg_feat[keep]
            cpg_to_region = cpg_to_region[keep]
        if cpg_feat.shape[0] == 0:
            raise ValueError(f"sample {record.sample_id!r} has zero edges after region-type filter")
        feat_parts.append(cpg_feat)
        region_parts.append(cpg_to_region + i * n_regions)
        type_parts.append(region_type.clone())
        gene_parts.append(region_to_gene + i * n_genes)

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

    return HierBatch(
        sample_ids=sample_ids,
        cpg_features=torch.cat(feat_parts, dim=0),
        cpg_to_region=torch.cat(region_parts, dim=0),
        region_type=torch.cat(type_parts, dim=0),
        region_to_gene=torch.cat(gene_parts, dim=0),
        n_regions=n_regions,
        n_genes=n_genes,
        tissue_target=torch.tensor(tissue_targets, dtype=torch.long),
        tissue_mask=torch.tensor(tissue_masks, dtype=torch.bool),
        age_target=torch.tensor(age_targets, dtype=torch.float32),
        age_mask=torch.tensor(age_masks, dtype=torch.bool),
        sex_target=torch.tensor(sex_targets, dtype=torch.long),
        sex_mask=torch.tensor(sex_masks, dtype=torch.bool),
    )


def make_synthetic_hier_overfit_bundle(
    *,
    n_samples: int = 12,
    n_cpgs: int = 12,
    n_genes: int = 3,
    n_classes: int = 3,
    static_dim: int = 4,
    seed: int = 0,
    include_unassigned: bool = True,
) -> dict[str, Any]:
    """Tiny hierarchical fixture: one typed region per gene + optional orphans."""
    if n_genes < n_classes:
        raise ValueError("n_genes must be >= n_classes for the overfit fixture")
    rng = np.random.default_rng(seed)

    # First n_genes CpGs: typed gene_body regions; remaining: unassigned singletons.
    n_typed = n_genes
    n_orphan = max(0, n_cpgs - n_typed) if include_unassigned else 0
    n_cpgs_eff = n_typed + n_orphan

    gene_ids = [f"GENE{i}" for i in range(n_genes)]
    if n_orphan > 0:
        gene_ids = [*gene_ids, UNASSIGNED_GENE_ID]

    region_ids: list[str] = []
    region_type_id: list[int] = []
    region_to_gene: list[int] = []
    edge_col: list[int] = []
    edge_region: list[int] = []

    body_id = REGION_TYPE_TO_ID["gene_body"]
    promoter_id = REGION_TYPE_TO_ID["promoter_core"]
    for g in range(n_genes):
        body_rid = f"GENE{g}:gene_body"
        prom_rid = f"GENE{g}:promoter_core"
        body_ridx = len(region_ids)
        region_ids.append(body_rid)
        region_type_id.append(body_id)
        region_to_gene.append(g)
        edge_col.append(g)
        edge_region.append(body_ridx)
        prom_ridx = len(region_ids)
        region_ids.append(prom_rid)
        region_type_id.append(promoter_id)
        region_to_gene.append(g)
        edge_col.append(g)
        edge_region.append(prom_ridx)

    n_typed_edges = len(edge_col)
    unassigned_type = REGION_TYPE_TO_ID[UNASSIGNED_REGION_TYPE]
    unassigned_gene_idx = gene_ids.index(UNASSIGNED_GENE_ID) if n_orphan > 0 else None
    for o in range(n_orphan):
        col = n_genes + o
        ridx = len(region_ids)
        region_ids.append(f"unassigned:L{col}")
        region_type_id.append(unassigned_type)
        if unassigned_gene_idx is None:
            raise RuntimeError("unassigned gene index missing")
        region_to_gene.append(unassigned_gene_idx)
        edge_col.append(col)
        edge_region.append(ridx)

    locus_region = LocusRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=np.asarray(edge_col, dtype=np.int64),
        edge_region_index=np.asarray(edge_region, dtype=np.int64),
        region_type_id=np.asarray(region_type_id, dtype=np.int64),
        region_to_gene=np.asarray(region_to_gene, dtype=np.int64),
        region_ids=region_ids,
        n_study_loci=n_cpgs_eff,
        n_typed_edges=n_typed_edges,
        n_unassigned_regions=n_orphan,
    )

    class_names = [f"class_{i}" for i in range(n_classes)]
    records: list[HierSampleRecord] = []
    ages: list[float] = []
    for s in range(n_samples):
        cls = s % n_classes
        betas = np.full(n_cpgs_eff, 0.05, dtype=np.float32)
        betas[cls] = 0.95
        betas = np.clip(betas + rng.normal(0, 0.005, size=n_cpgs_eff), 0.01, 0.99).astype(
            np.float32
        )
        static = rng.normal(0, 0.01, size=(n_cpgs_eff, static_dim)).astype(np.float32)
        static[cls, 0] = 1.0 + float(cls)
        static_valid = np.ones(n_cpgs_eff, dtype=bool)
        features = gather_hier_sample_features(
            beta_row=betas,
            static_by_col=static,
            static_valid=static_valid,
            locus_region=locus_region,
        )
        records.append(
            HierSampleRecord(
                sample_id=f"SYN{s:03d}",
                donor_id=str((s % 4) + 1),
                class_index=cls,
                features=features,
            )
        )
        ages.append(20.0 + 10.0 * cls)

    return {
        "records": records,
        "locus_region": locus_region,
        "gene_ids": gene_ids,
        "class_names": class_names,
        "ages": ages,
        "n_genes": len(gene_ids),
        "n_regions": locus_region.n_regions,
        "input_dim": 2 + static_dim,
        "n_classes": n_classes,
        "region_types": list(HIER_REGION_TYPES),
    }
