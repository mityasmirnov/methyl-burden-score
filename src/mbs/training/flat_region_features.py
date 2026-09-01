"""FlatDeepSetRegion feature assembly: M-value + region-type one-hot + observed."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mbs.training.cascade_assign import CascadeAssignment, assignment_gene_linked_only
from mbs.training.features import beta_to_m_value


@dataclass(frozen=True, slots=True)
class FlatRegionGeneIndex:
    """Gene-linked typed edges for flat region pooling."""

    gene_ids: list[str]
    region_types: tuple[str, ...]
    edge_col_index: np.ndarray  # int64 [n_edges]
    edge_gene_index: np.ndarray  # int64 [n_edges]
    edge_type_id: np.ndarray  # int64 [n_edges]
    n_study_loci: int

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_edges(self) -> int:
        return int(self.edge_col_index.shape[0])


def flat_region_input_dim(n_region_types: int) -> int:
    """M-value + type one-hot + observed flag."""
    if n_region_types <= 0:
        raise ValueError("n_region_types must be positive")
    return 1 + int(n_region_types) + 1


def build_flat_region_gene_index(assignment: CascadeAssignment) -> FlatRegionGeneIndex:
    """Build index from a gene-linked cascade assignment."""
    linked = assignment_gene_linked_only(assignment)
    if linked.edge_col_index.size == 0:
        raise ValueError("no gene-linked edges for flat region index")
    gene_edge = linked.region_to_gene[linked.edge_region_index] >= 0
    if not np.any(gene_edge):
        raise ValueError("no gene-linked edges after filter")
    cols = linked.edge_col_index[gene_edge]
    regs = linked.edge_region_index[gene_edge]
    genes = linked.region_to_gene[regs]
    types = linked.region_type_id[regs]
    return FlatRegionGeneIndex(
        gene_ids=list(linked.gene_ids),
        region_types=linked.region_types,
        edge_col_index=cols.astype(np.int64, copy=False),
        edge_gene_index=genes.astype(np.int64, copy=False),
        edge_type_id=types.astype(np.int64, copy=False),
        n_study_loci=linked.n_study_loci,
    )


def gather_flat_region_features(
    *,
    beta_row: np.ndarray,
    index: FlatRegionGeneIndex,
    epsilon: float = 0.001,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(features [n_edges, dim], cpg_to_gene [n_edges])``."""
    betas = np.asarray(beta_row, dtype=np.float32).reshape(-1)
    n_types = len(index.region_types)
    dim = flat_region_input_dim(n_types)
    n_edges = index.n_edges
    if n_edges == 0:
        return np.zeros((0, dim), dtype=np.float32), np.zeros(0, dtype=np.int64)
    feats = np.zeros((n_edges, dim), dtype=np.float32)
    cols = index.edge_col_index
    obs = np.isfinite(betas[cols])
    safe = np.where(obs, betas[cols], 0.5)
    m_vals = beta_to_m_value(safe, epsilon=epsilon)
    m_vals = np.where(obs, m_vals, 0.0).astype(np.float32)
    feats[:, 0] = m_vals
    type_ids = index.edge_type_id.astype(np.int64)
    for i in range(n_edges):
        tid = int(type_ids[i])
        if 0 <= tid < n_types:
            feats[i, 1 + tid] = 1.0
    feats[:, -1] = obs.astype(np.float32)
    return feats, index.edge_gene_index.astype(np.int64, copy=False)
