"""Epoch-varying within-gene CpG sampler (Milestone 12b §3.5).

Train may expose ≤K CpGs per gene, resampling each epoch so large genes are
not starved. Validation and product scoring must use the full edge index
(never call this sampler).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from mbs.training.flat_region_features import GENE_ROLES, FlatRegionGeneIndex


def _subset_edges(index: FlatRegionGeneIndex, keep: np.ndarray) -> FlatRegionGeneIndex:
    keep = np.asarray(keep, dtype=np.int64)
    role_ids = np.asarray(index.edge_role_id, dtype=np.int64)[keep]
    other_id = GENE_ROLES.index("other_gene")
    return FlatRegionGeneIndex(
        gene_ids=list(index.gene_ids),
        edge_col_index=np.asarray(index.edge_col_index, dtype=np.int64)[keep],
        edge_gene_index=np.asarray(index.edge_gene_index, dtype=np.int64)[keep],
        edge_role_id=role_ids,
        edge_context_id=np.asarray(index.edge_context_id, dtype=np.int64)[keep],
        edge_role_present=np.asarray(index.edge_role_present, dtype=bool)[keep],
        edge_context_present=np.asarray(index.edge_context_present, dtype=bool)[keep],
        edge_regulatory_present=np.asarray(index.edge_regulatory_present, dtype=bool)[
            keep
        ],
        edge_regulatory_multi_hot=np.asarray(
            index.edge_regulatory_multi_hot, dtype=np.float32
        )[keep],
        n_study_loci=int(index.n_study_loci),
        n_other_gene_edges=int(np.sum(role_ids == other_id)),
    )


@dataclass
class WithinGeneEpochSampler:
    """Resample ≤K edges per gene each epoch; track locus/edge exposure."""

    full_index: FlatRegionGeneIndex
    max_cpgs_per_gene: int
    seed: int = 0
    times_available: np.ndarray = field(init=False)
    times_selected: np.ndarray = field(init=False)
    _gene_edge_idxs: list[np.ndarray] = field(init=False, repr=False)
    _epochs_seen: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        k = int(self.max_cpgs_per_gene)
        if k < 1:
            raise ValueError(f"max_cpgs_per_gene must be >= 1, got {self.max_cpgs_per_gene}")
        n = int(self.full_index.n_edges)
        self.times_available = np.zeros(n, dtype=np.int64)
        self.times_selected = np.zeros(n, dtype=np.int64)
        genes = np.asarray(self.full_index.edge_gene_index, dtype=np.int64)
        self._gene_edge_idxs = [
            np.flatnonzero(genes == g).astype(np.int64, copy=False)
            for g in range(self.full_index.n_genes)
        ]

    def sample_epoch(self, epoch: int) -> FlatRegionGeneIndex:
        """Deterministic sample for ``epoch`` (1-based training epochs)."""
        k = int(self.max_cpgs_per_gene)
        rng = np.random.default_rng(int(self.seed) + int(epoch) * 1_000_003)
        keep_parts: list[np.ndarray] = []
        for idxs in self._gene_edge_idxs:
            if idxs.size == 0:
                continue
            self.times_available[idxs] += 1
            chosen = idxs if idxs.size <= k else rng.choice(idxs, size=k, replace=False)
            keep_parts.append(np.sort(chosen))
            self.times_selected[chosen] += 1
        self._epochs_seen += 1
        if not keep_parts:
            raise ValueError("within-gene sampler produced zero edges")
        keep = np.concatenate(keep_parts)
        return _subset_edges(self.full_index, keep)

    def exposure_report(self) -> dict[str, Any]:
        avail = np.maximum(self.times_available, 1)
        prob = self.times_selected.astype(np.float64) / avail.astype(np.float64)
        cols = np.asarray(self.full_index.edge_col_index, dtype=np.int64)
        return {
            "max_cpgs_per_gene": int(self.max_cpgs_per_gene),
            "seed": int(self.seed),
            "n_epochs": int(self._epochs_seen),
            "n_edges_full": int(self.full_index.n_edges),
            "n_genes": int(self.full_index.n_genes),
            "times_selected_sum": int(self.times_selected.sum()),
            "mean_selection_probability": float(prob.mean()) if prob.size else 0.0,
            "per_edge": {
                "edge_col_index": cols.tolist(),
                "times_available": self.times_available.tolist(),
                "times_selected": self.times_selected.tolist(),
                "selection_probability": prob.tolist(),
            },
        }
