"""FlatDeepSetRegion feature assembly: annotated CpG → gene one-hop."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from mbs.annotation.gencode_regions import REGION_TYPES as GENE_ROLE_TYPES
from mbs.training.cascade_assign import CascadeAssignment, assignment_gene_linked_only
from mbs.training.features import beta_to_m_value

GENE_ROLES: tuple[str, ...] = (*GENE_ROLE_TYPES, "other_gene")
CPG_CONTEXTS: tuple[str, ...] = (
    "island",
    "north_shore",
    "south_shore",
    "north_shelf",
    "south_shelf",
    "open_sea",
    "unknown",
)
# Reserved multi-hot slots; graph-v2 has no SCREEN/cCRE yet → stay zero.
REGULATORY_CHANNELS: tuple[str, ...] = (
    "promoter_like_ccre",
    "enhancer_like_ccre",
    "ctcf_only_ccre",
    "dhs_only_ccre",
    "chromhmm",
    "dhs_flag",
)
PRESENCE_FLAGS: tuple[str, ...] = (
    "gene_role_present",
    "cpg_context_present",
    "regulatory_annotation_present",
)


@dataclass(frozen=True, slots=True)
class FlatRegionGeneIndex:
    """Gene-linked annotated edges for flat gene pooling."""

    gene_ids: list[str]
    edge_col_index: np.ndarray  # int64 [n_edges]
    edge_gene_index: np.ndarray  # int64 [n_edges]
    edge_role_id: np.ndarray  # int64 [n_edges] into GENE_ROLES
    edge_context_id: np.ndarray  # int64 [n_edges] into CPG_CONTEXTS
    edge_role_present: np.ndarray  # bool [n_edges]
    edge_context_present: np.ndarray  # bool [n_edges]
    edge_regulatory_present: np.ndarray  # bool [n_edges]
    edge_regulatory_multi_hot: np.ndarray  # float32 [n_edges, n_reg]
    n_study_loci: int
    n_other_gene_edges: int

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_edges(self) -> int:
        return int(self.edge_col_index.shape[0])

    @property
    def region_types(self) -> tuple[str, ...]:
        """Backward-compatible alias used by older callers."""
        return GENE_ROLES


def flat_region_input_dim(
    n_region_types: int | None = None,
    *,
    n_gene_roles: int = len(GENE_ROLES),
    n_cpg_contexts: int = len(CPG_CONTEXTS),
    n_regulatory: int = len(REGULATORY_CHANNELS),
) -> int:
    """M-value + gene-role + CGI context + regulatory multi-hot + flags + observed."""
    del n_region_types  # legacy kwarg ignored; dim is fixed by channel tables
    if n_gene_roles <= 0 or n_cpg_contexts <= 0:
        raise ValueError("gene-role and cpg-context dims must be positive")
    # M + roles + contexts + regulatory + 3 presence flags + observed
    return 1 + int(n_gene_roles) + int(n_cpg_contexts) + int(n_regulatory) + 3 + 1


def count_other_gene_edges(
    region_type_ids: np.ndarray,
    region_types: tuple[str, ...],
) -> int:
    """Count edges whose region type is not one of the five GENCODE roles."""
    role_set = set(GENE_ROLE_TYPES)
    n = 0
    for tid in np.asarray(region_type_ids, dtype=np.int64).tolist():
        if 0 <= int(tid) < len(region_types) and region_types[int(tid)] not in role_set:
            n += 1
    return int(n)


def assert_other_gene_count(
    *,
    n_other_gene_edges: int,
    allow_nonzero: bool = False,
) -> None:
    """Stage A graph check: five-role graphs should report zero ``other_gene``."""
    if n_other_gene_edges and not allow_nonzero:
        raise AssertionError(
            f"expected zero other_gene edges on five-role graph; found {n_other_gene_edges}"
        )


def _role_id_for_type(region_type: str) -> tuple[int, bool]:
    if region_type in GENE_ROLE_TYPES:
        return GENE_ROLES.index(region_type), True
    return GENE_ROLES.index("other_gene"), True


def _context_id(label: str | None) -> tuple[int, bool]:
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return CPG_CONTEXTS.index("unknown"), False
    text = str(label).strip().lower()
    if not text or text in {"nan", "none", "null"}:
        return CPG_CONTEXTS.index("unknown"), False
    if text in CPG_CONTEXTS:
        return CPG_CONTEXTS.index(text), True
    return CPG_CONTEXTS.index("unknown"), False


def build_flat_region_gene_index(
    assignment: CascadeAssignment,
    *,
    locus_index: pd.DataFrame | None = None,
    cpg_context_by_locus: dict[str, str] | None = None,
    allow_other_gene: bool = False,
) -> FlatRegionGeneIndex:
    """Build annotated index from a gene-linked cascade assignment."""
    linked = assignment_gene_linked_only(assignment)
    if linked.edge_col_index.size == 0:
        raise ValueError("no gene-linked edges for flat region index")
    gene_edge = linked.region_to_gene[linked.edge_region_index] >= 0
    if not np.any(gene_edge):
        raise ValueError("no gene-linked edges after filter")
    cols = linked.edge_col_index[gene_edge].astype(np.int64, copy=False)
    regs = linked.edge_region_index[gene_edge].astype(np.int64, copy=False)
    genes = linked.region_to_gene[regs].astype(np.int64, copy=False)
    type_ids = linked.region_type_id[regs].astype(np.int64, copy=False)

    n_edges = int(cols.shape[0])
    role_ids = np.zeros(n_edges, dtype=np.int64)
    role_present = np.zeros(n_edges, dtype=bool)
    context_ids = np.full(n_edges, CPG_CONTEXTS.index("unknown"), dtype=np.int64)
    context_present = np.zeros(n_edges, dtype=bool)
    regulatory = np.zeros((n_edges, len(REGULATORY_CHANNELS)), dtype=np.float32)
    regulatory_present = np.zeros(n_edges, dtype=bool)

    locus_ids: list[str] | None = None
    if locus_index is not None and "locus_id" in locus_index.columns:
        locus_ids = [str(x) for x in locus_index["locus_id"].tolist()]

    for i in range(n_edges):
        tid = int(type_ids[i])
        rtype = (
            linked.region_types[tid]
            if 0 <= tid < len(linked.region_types)
            else "other_gene"
        )
        rid, rpres = _role_id_for_type(rtype)
        role_ids[i] = rid
        role_present[i] = rpres
        ctx_label = None
        if cpg_context_by_locus is not None and locus_ids is not None:
            col = int(cols[i])
            if 0 <= col < len(locus_ids):
                ctx_label = cpg_context_by_locus.get(locus_ids[col])
        elif locus_index is not None and "cpg_context" in locus_index.columns:
            col = int(cols[i])
            if 0 <= col < len(locus_index):
                ctx_label = locus_index.iloc[col]["cpg_context"]
        cid, cpres = _context_id(ctx_label)
        context_ids[i] = cid
        context_present[i] = cpres

    n_other = int(np.sum(role_ids == GENE_ROLES.index("other_gene")))
    assert_other_gene_count(n_other_gene_edges=n_other, allow_nonzero=allow_other_gene)

    return FlatRegionGeneIndex(
        gene_ids=list(linked.gene_ids),
        edge_col_index=cols,
        edge_gene_index=genes,
        edge_role_id=role_ids,
        edge_context_id=context_ids,
        edge_role_present=role_present,
        edge_context_present=context_present,
        edge_regulatory_present=regulatory_present,
        edge_regulatory_multi_hot=regulatory,
        n_study_loci=linked.n_study_loci,
        n_other_gene_edges=n_other,
    )


def gather_flat_region_features(
    *,
    beta_row: np.ndarray,
    index: FlatRegionGeneIndex,
    epsilon: float = 0.001,
    base_features: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(features [n_edges, dim], cpg_to_gene [n_edges])``.

    When ``base_features`` is provided (annotation channels cached once per fold),
    only the M-value and observed flag are refreshed from ``beta_row``.
    """
    betas = np.asarray(beta_row, dtype=np.float32).reshape(-1)
    dim = flat_region_input_dim()
    n_edges = index.n_edges
    if n_edges == 0:
        return np.zeros((0, dim), dtype=np.float32), np.zeros(0, dtype=np.int64)
    cols = index.edge_col_index
    obs = np.isfinite(betas[cols])
    safe = np.where(obs, betas[cols], 0.5)
    m_vals = beta_to_m_value(safe, epsilon=epsilon)
    m_vals = np.where(obs, m_vals, 0.0).astype(np.float32)
    if base_features is not None:
        feats = np.array(base_features, dtype=np.float32, copy=True)
        if feats.shape != (n_edges, dim):
            raise ValueError(
                f"base_features shape {feats.shape} != expected {(n_edges, dim)}"
            )
        feats[:, 0] = m_vals
        feats[:, -1] = obs.astype(np.float32)
        return feats, index.edge_gene_index.astype(np.int64, copy=False)
    feats = np.zeros((n_edges, dim), dtype=np.float32)
    feats[:, 0] = m_vals
    offset = 1
    role_ids = np.asarray(index.edge_role_id, dtype=np.int64)
    valid_role = (role_ids >= 0) & (role_ids < len(GENE_ROLES))
    if np.any(valid_role):
        feats[np.nonzero(valid_role)[0], offset + role_ids[valid_role]] = 1.0
    offset += len(GENE_ROLES)
    context_ids = np.asarray(index.edge_context_id, dtype=np.int64)
    valid_ctx = (context_ids >= 0) & (context_ids < len(CPG_CONTEXTS))
    if np.any(valid_ctx):
        feats[np.nonzero(valid_ctx)[0], offset + context_ids[valid_ctx]] = 1.0
    offset += len(CPG_CONTEXTS)
    feats[:, offset : offset + len(REGULATORY_CHANNELS)] = index.edge_regulatory_multi_hot
    offset += len(REGULATORY_CHANNELS)
    feats[:, offset] = index.edge_role_present.astype(np.float32)
    feats[:, offset + 1] = index.edge_context_present.astype(np.float32)
    feats[:, offset + 2] = index.edge_regulatory_present.astype(np.float32)
    feats[:, -1] = obs.astype(np.float32)
    return feats, index.edge_gene_index.astype(np.int64, copy=False)


def build_flat_region_base_features(index: FlatRegionGeneIndex) -> np.ndarray:
    """Annotation-only feature template (M-value/observed filled per sample)."""
    dim = flat_region_input_dim()
    n_edges = index.n_edges
    feats = np.zeros((n_edges, dim), dtype=np.float32)
    if n_edges == 0:
        return feats
    offset = 1
    role_ids = np.asarray(index.edge_role_id, dtype=np.int64)
    valid_role = (role_ids >= 0) & (role_ids < len(GENE_ROLES))
    if np.any(valid_role):
        feats[np.nonzero(valid_role)[0], offset + role_ids[valid_role]] = 1.0
    offset += len(GENE_ROLES)
    context_ids = np.asarray(index.edge_context_id, dtype=np.int64)
    valid_ctx = (context_ids >= 0) & (context_ids < len(CPG_CONTEXTS))
    if np.any(valid_ctx):
        feats[np.nonzero(valid_ctx)[0], offset + context_ids[valid_ctx]] = 1.0
    offset += len(CPG_CONTEXTS)
    feats[:, offset : offset + len(REGULATORY_CHANNELS)] = index.edge_regulatory_multi_hot
    offset += len(REGULATORY_CHANNELS)
    feats[:, offset] = index.edge_role_present.astype(np.float32)
    feats[:, offset + 1] = index.edge_context_present.astype(np.float32)
    feats[:, offset + 2] = index.edge_regulatory_present.astype(np.float32)
    return feats
