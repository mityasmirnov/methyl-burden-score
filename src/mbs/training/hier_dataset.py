"""Hierarchical sample records and packed batches for HierarchicalDeepSet."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import torch
from torch import Tensor

from mbs.batch import (
    ANNOTATION_STATUS_AMBIGUOUS,
    ANNOTATION_STATUS_MAPPED,
    ANNOTATION_STATUS_MULTI_MAPPED,
    ANNOTATION_STATUS_UNMAPPED,
    annotation_status_masks,
)
from mbs.training.features import beta_to_m_value
from mbs.training.locus_region_gene import (
    HIER_REGION_TYPES,
    REGION_TYPE_TO_ID,
    RESIDUAL_PANEL_ID,
    LocusRegionGeneIndex,
)
from mbs.training.phenotypes import SamplePhenotype


@dataclass(frozen=True, slots=True)
class HierSampleFeatureBundle:
    cpg_features: np.ndarray  # float32 [n_edges_obs, feat_dim] — mapped path
    cpg_to_region: np.ndarray  # int64 [n_edges_obs]
    edge_col_index: np.ndarray  # int64 [n_edges_obs] — study columns for status
    residual_features: np.ndarray  # float32 [n_residual_obs, feat_dim]
    n_observed_edges: int
    n_observed_residual: int
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
    """Packed multi-sample batch for ``HierarchicalDeepSet`` + residual path."""

    sample_ids: list[str]
    cpg_features: Tensor
    cpg_to_region: Tensor
    region_type: Tensor
    region_to_gene: Tensor
    residual_features: Tensor
    residual_sample_index: Tensor
    annotation_status: Tensor  # per mapped edge row (study-column status)
    n_regions: int
    n_genes: int
    tissue_target: Tensor
    tissue_mask: Tensor
    age_target: Tensor | None
    age_mask: Tensor
    sex_target: Tensor | None = None
    sex_mask: Tensor | None = None

    def annotation_masks(self) -> dict[str, Tensor]:
        return annotation_status_masks(self.annotation_status)

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
            residual_features=self.residual_features.to(device, non_blocking=non_blocking),
            residual_sample_index=self.residual_sample_index.to(device, non_blocking=non_blocking),
            annotation_status=self.annotation_status.to(device, non_blocking=non_blocking),
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


def _assemble_features(
    *,
    betas: np.ndarray,
    static_rows: np.ndarray,
    include_m_value: bool,
    epsilon: float,
) -> np.ndarray:
    parts: list[np.ndarray] = [betas.reshape(-1, 1)]
    if include_m_value:
        parts.append(beta_to_m_value(betas, epsilon=epsilon).reshape(-1, 1))
    parts.append(static_rows.astype(np.float32, copy=False))
    return np.concatenate(parts, axis=1).astype(np.float32, copy=False)


def gather_hier_sample_features(
    *,
    beta_row: np.ndarray,
    static_by_col: np.ndarray,
    static_valid: np.ndarray,
    locus_region: LocusRegionGeneIndex,
    epsilon: float = 0.001,
    include_m_value: bool = True,
) -> HierSampleFeatureBundle:
    """Build ragged hierarchical + residual features for one sample."""
    beta_row = np.asarray(beta_row, dtype=np.float32)
    if beta_row.ndim != 1:
        raise ValueError("beta_row must be 1-D")
    if beta_row.shape[0] < locus_region.n_study_loci:
        raise ValueError(
            f"beta_row length {beta_row.shape[0]} < n_study_loci {locus_region.n_study_loci}"
        )

    cols = locus_region.edge_col_index
    regions = locus_region.edge_region_index
    n_dropped_nan = 0
    n_dropped_static = 0

    if cols.size == 0:
        feat_dim = 2 + static_by_col.shape[1] if include_m_value else 1 + static_by_col.shape[1]
        mapped_features = np.zeros((0, feat_dim), dtype=np.float32)
        mapped_regions = np.zeros(0, dtype=np.int64)
        mapped_cols = np.zeros(0, dtype=np.int64)
    else:
        betas = beta_row[cols]
        finite = np.isfinite(betas)
        static_ok = static_valid[cols]
        keep = finite & static_ok
        n_dropped_nan += int((~finite).sum())
        n_dropped_static += int((finite & ~static_ok).sum())
        cols_k = cols[keep]
        regions_k = regions[keep]
        mapped_features = _assemble_features(
            betas=betas[keep],
            static_rows=static_by_col[cols_k],
            include_m_value=include_m_value,
            epsilon=epsilon,
        )
        mapped_regions = regions_k.astype(np.int64, copy=False)
        mapped_cols = cols_k.astype(np.int64, copy=False)

    res_cols = locus_region.residual_col_index
    if res_cols.size == 0:
        feat_dim = (
            mapped_features.shape[1]
            if mapped_features.size
            else (2 + static_by_col.shape[1] if include_m_value else 1 + static_by_col.shape[1])
        )
        residual_features = np.zeros((0, feat_dim), dtype=np.float32)
    else:
        betas_r = beta_row[res_cols]
        finite_r = np.isfinite(betas_r)
        # Residual may lack CpGPT rows; zero-fill static when missing.
        static_ok_r = static_valid[res_cols]
        n_dropped_nan += int((~finite_r).sum())
        keep_r = finite_r
        cols_r = res_cols[keep_r]
        betas_rk = betas_r[keep_r]
        static_r = static_by_col[cols_r].astype(np.float32, copy=True)
        static_r[~static_ok_r[keep_r]] = 0.0
        residual_features = _assemble_features(
            betas=betas_rk,
            static_rows=static_r,
            include_m_value=include_m_value,
            epsilon=epsilon,
        )

    return HierSampleFeatureBundle(
        cpg_features=mapped_features,
        cpg_to_region=mapped_regions,
        edge_col_index=mapped_cols,
        residual_features=residual_features,
        n_observed_edges=int(mapped_features.shape[0]),
        n_observed_residual=int(residual_features.shape[0]),
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
    if features.n_observed_edges == 0 and features.n_observed_residual == 0:
        raise ValueError(f"sample {phenotype.sample_id!r} has zero observed CpGs")
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
    include_residual: bool = True,
) -> HierBatch:
    """Pack ragged samples with region/gene offsets for one hierarchical forward.

    Optional ``allowed_region_type_ids`` filters mapped-path edges (eval slices).
    ``include_residual=False`` drops residual features (mapped_only slice).
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
    col_status = locus_region.column_annotation_status

    feat_parts: list[Tensor] = []
    region_parts: list[Tensor] = []
    status_parts: list[Tensor] = []
    type_parts: list[Tensor] = []
    gene_parts: list[Tensor] = []
    residual_parts: list[Tensor] = []
    residual_sample_parts: list[Tensor] = []
    tissue_targets: list[int] = []
    tissue_masks: list[bool] = []
    age_targets: list[float] = []
    age_masks: list[bool] = []
    sex_targets: list[int] = []
    sex_masks: list[bool] = []
    sample_ids: list[str] = []

    allow = allowed_region_type_ids
    allow_np: np.ndarray | None = None
    if allow is not None:
        allow_np = np.zeros(len(HIER_REGION_TYPES), dtype=np.bool_)
        for tid in allow:
            t = int(tid)
            if 0 <= t < allow_np.shape[0]:
                allow_np[t] = True
    feat_dim: int | None = None
    for i, record in enumerate(records):
        feats = record.features
        cpg_feat_np = feats.cpg_features
        cpg_to_region_np = feats.cpg_to_region.astype(np.int64, copy=False)
        edge_cols_np = feats.edge_col_index.astype(np.int64, copy=False)
        if allow_np is not None and cpg_to_region_np.size > 0:
            edge_types = locus_region.region_type_id[cpg_to_region_np]
            keep = allow_np[edge_types]
            cpg_feat_np = cpg_feat_np[keep]
            cpg_to_region_np = cpg_to_region_np[keep]
            edge_cols_np = edge_cols_np[keep]
        residual_feat_np = feats.residual_features
        if not include_residual:
            residual_feat_np = residual_feat_np[:0]
        if cpg_feat_np.shape[0] == 0 and residual_feat_np.shape[0] == 0:
            raise ValueError(f"sample {record.sample_id!r} has zero edges after filters")
        if feat_dim is None:
            feat_dim = int(
                cpg_feat_np.shape[1] if cpg_feat_np.shape[0] else residual_feat_np.shape[1]
            )
        if cpg_feat_np.shape[0]:
            status_np = col_status[edge_cols_np].astype(np.int64, copy=False)
            status_parts.append(torch.from_numpy(np.ascontiguousarray(status_np)))
            feat_parts.append(torch.from_numpy(np.ascontiguousarray(cpg_feat_np)))
            region_parts.append(
                torch.from_numpy(
                    np.ascontiguousarray(cpg_to_region_np + i * max(n_regions, 1))
                )
            )
        if residual_feat_np.shape[0]:
            residual_parts.append(torch.from_numpy(np.ascontiguousarray(residual_feat_np)))
            residual_sample_parts.append(
                torch.full((residual_feat_np.shape[0],), i, dtype=torch.long)
            )
        if n_regions > 0:
            type_parts.append(region_type.clone())
            gene_parts.append(region_to_gene + i * max(n_genes, 1))

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

    empty_feat = torch.zeros(0, feat_dim or 1, dtype=torch.float32)
    empty_long = torch.zeros(0, dtype=torch.long)
    return HierBatch(
        sample_ids=sample_ids,
        cpg_features=torch.cat(feat_parts, dim=0) if feat_parts else empty_feat,
        cpg_to_region=torch.cat(region_parts, dim=0) if region_parts else empty_long,
        region_type=torch.cat(type_parts, dim=0) if type_parts else empty_long,
        region_to_gene=torch.cat(gene_parts, dim=0) if gene_parts else empty_long,
        residual_features=torch.cat(residual_parts, dim=0) if residual_parts else empty_feat,
        residual_sample_index=(
            torch.cat(residual_sample_parts, dim=0) if residual_sample_parts else empty_long
        ),
        annotation_status=torch.cat(status_parts, dim=0) if status_parts else empty_long,
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
    include_residual: bool = True,
) -> dict[str, Any]:
    """Tiny hierarchical fixture: typed regions + optional residual columns."""
    if n_genes < n_classes:
        raise ValueError("n_genes must be >= n_classes for the overfit fixture")
    rng = np.random.default_rng(seed)

    n_typed = n_genes
    n_orphan = max(0, n_cpgs - n_typed) if include_residual else 0
    n_cpgs_eff = n_typed + n_orphan

    gene_ids = [f"GENE{i}" for i in range(n_genes)]
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
    residual_cols = list(range(n_genes, n_cpgs_eff))
    status = np.full(n_cpgs_eff, ANNOTATION_STATUS_UNMAPPED, dtype=np.int8)
    status[:n_genes] = ANNOTATION_STATUS_MAPPED
    # Mark one residual as ambiguous for mask coverage when present.
    if residual_cols:
        status[residual_cols[0]] = ANNOTATION_STATUS_AMBIGUOUS
        if len(residual_cols) > 1:
            status[residual_cols[1]] = ANNOTATION_STATUS_UNMAPPED

    locus_region = LocusRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=np.asarray(edge_col, dtype=np.int64),
        edge_region_index=np.asarray(edge_region, dtype=np.int64),
        region_type_id=np.asarray(region_type_id, dtype=np.int64),
        region_to_gene=np.asarray(region_to_gene, dtype=np.int64),
        region_ids=region_ids,
        residual_col_index=np.asarray(residual_cols, dtype=np.int64),
        column_annotation_status=status,
        n_study_loci=n_cpgs_eff,
        n_typed_edges=n_typed_edges,
    )

    class_names = [f"class_{i}" for i in range(n_classes)]
    records: list[HierSampleRecord] = []
    ages: list[float] = []
    for s in range(n_samples):
        cls = s % n_classes
        betas = np.full(n_cpgs_eff, 0.05, dtype=np.float32)
        betas[cls] = 0.95
        if residual_cols:
            betas[residual_cols[0]] = 0.9 if cls == 0 else 0.1
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
        "panel_ids": [*gene_ids, RESIDUAL_PANEL_ID],
        "class_names": class_names,
        "ages": ages,
        "n_genes": len(gene_ids),
        "n_panel": locus_region.n_panel,
        "n_regions": locus_region.n_regions,
        "input_dim": 2 + static_dim,
        "n_classes": n_classes,
        "region_types": list(HIER_REGION_TYPES),
        "annotation_statuses": {
            "mapped": ANNOTATION_STATUS_MAPPED,
            "unmapped": ANNOTATION_STATUS_UNMAPPED,
            "ambiguous": ANNOTATION_STATUS_AMBIGUOUS,
            "multi_mapped": ANNOTATION_STATUS_MULTI_MAPPED,
        },
    }
