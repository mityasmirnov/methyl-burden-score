"""Gene-holdout partition for DeepRVAT-style encoder transfer (GATE G1).

Train φ/ρ on one gene set, score a disjoint heldout set with the frozen
encoder, then fit the post-hoc nested elastic-net on heldout MBS columns the
encoder never saw during training.

Methods:
- ``random`` — per-gene shuffle (locality leakage possible across the boundary).
- ``chromosome`` — hold out whole chromosomes (bounds neighbour leakage).
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from mbs.training.cascade_assign import ORPHAN_GENE_INDEX, CascadeAssignment
from mbs.training.flat_region_features import GENE_ROLES, FlatRegionGeneIndex

GeneHoldoutMethod = Literal["random", "chromosome"]

_RANDOM_CAVEAT = (
    "random split: neighbouring/co-regulated genes may straddle the "
    "boundary (locality leakage). Chromosome-held-out arm bounds that."
)
_CHROM_CAVEAT = (
    "chromosome holdout: entire chromosomes are held out so neighbouring "
    "genes cannot straddle the train/heldout boundary. Sex chromosomes may "
    "be included unless excluded via config."
)


@dataclass(frozen=True, slots=True)
class GeneHoldoutPartition:
    """Disjoint train / heldout gene id lists (string gene ids from the index)."""

    method: GeneHoldoutMethod
    seed: int
    heldout_fraction: float
    train_gene_ids: list[str]
    heldout_gene_ids: list[str]
    heldout_chromosomes: list[str] = field(default_factory=list)
    train_chromosomes: list[str] = field(default_factory=list)

    @property
    def n_train(self) -> int:
        return len(self.train_gene_ids)

    @property
    def n_heldout(self) -> int:
        return len(self.heldout_gene_ids)

    @property
    def actual_heldout_fraction(self) -> float:
        n = self.n_train + self.n_heldout
        return float(self.n_heldout) / float(n) if n else 0.0

    def caveat(self) -> str:
        return _CHROM_CAVEAT if self.method == "chromosome" else _RANDOM_CAVEAT

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["actual_heldout_fraction"] = self.actual_heldout_fraction
        out["caveat"] = self.caveat()
        return out


def gene_chromosome_map(genes: pd.DataFrame) -> dict[str, str]:
    """``gene_id → chromosome`` from a genes.parquet-like table."""
    if genes.empty or "gene_id" not in genes.columns or "chromosome" not in genes.columns:
        raise ValueError("genes table must have gene_id and chromosome columns")
    out: dict[str, str] = {}
    for row in genes.itertuples(index=False):
        gid = getattr(row, "gene_id", None)
        chrom = getattr(row, "chromosome", None)
        if gid is None or chrom is None or (isinstance(chrom, float) and np.isnan(chrom)):
            continue
        out[str(gid)] = str(chrom)
    if not out:
        raise ValueError("no gene_id/chromosome pairs in genes table")
    return out


def partition_genes_random(
    gene_ids: Sequence[str],
    *,
    heldout_fraction: float = 0.2,
    seed: int = 42,
) -> GeneHoldoutPartition:
    """Random gene split. Neighbouring genes may straddle the boundary (locality leakage)."""
    if not (0.0 < float(heldout_fraction) < 1.0):
        raise ValueError(f"heldout_fraction must be in (0, 1); got {heldout_fraction}")
    ids = [str(g) for g in gene_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("gene_ids must be unique")
    if len(ids) < 2:
        raise ValueError("need at least 2 genes to partition")
    rng = np.random.default_rng(int(seed))
    order = rng.permutation(len(ids))
    n_heldout = max(1, round(len(ids) * float(heldout_fraction)))
    n_heldout = min(n_heldout, len(ids) - 1)
    held_idx = {int(i) for i in order[:n_heldout].tolist()}
    heldout = [ids[i] for i in range(len(ids)) if i in held_idx]
    train = [ids[i] for i in range(len(ids)) if i not in held_idx]
    if not train or not heldout:
        raise RuntimeError("gene partition produced an empty train or heldout set")
    return GeneHoldoutPartition(
        method="random",
        seed=int(seed),
        heldout_fraction=float(heldout_fraction),
        train_gene_ids=train,
        heldout_gene_ids=heldout,
    )


def partition_genes_by_chromosome(
    gene_ids: Sequence[str],
    gene_chromosome: Mapping[str, str],
    *,
    heldout_fraction: float = 0.2,
    seed: int = 42,
    heldout_chromosomes: Sequence[str] | None = None,
    exclude_chromosomes: Sequence[str] | None = None,
) -> GeneHoldoutPartition:
    """Hold out whole chromosomes so neighbours cannot straddle the split.

    When ``heldout_chromosomes`` is omitted, chromosomes are shuffled with
    ``seed`` and greedily accumulated until the heldout gene count reaches
    ``heldout_fraction`` (at least one chromosome, leaving ≥1 train gene).
    Genes missing from ``gene_chromosome`` stay in train.
    """
    if not (0.0 < float(heldout_fraction) < 1.0):
        raise ValueError(f"heldout_fraction must be in (0, 1); got {heldout_fraction}")
    ids = [str(g) for g in gene_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("gene_ids must be unique")
    if len(ids) < 2:
        raise ValueError("need at least 2 genes to partition")

    exclude = {str(c) for c in (exclude_chromosomes or ())}
    chrom_of = {g: str(gene_chromosome[g]) for g in ids if g in gene_chromosome}
    by_chrom: dict[str, list[str]] = {}
    for g, c in chrom_of.items():
        if c in exclude:
            continue
        by_chrom.setdefault(c, []).append(g)
    if not by_chrom:
        raise ValueError("no chromosomes available after excludes / missing map")

    if heldout_chromosomes is not None:
        held_chroms = [str(c) for c in heldout_chromosomes]
        missing = [c for c in held_chroms if c not in by_chrom]
        if missing:
            raise KeyError(f"heldout_chromosomes not present: {missing[:3]}")
    else:
        rng = np.random.default_rng(int(seed))
        chroms = list(by_chrom.keys())
        rng.shuffle(chroms)
        target = max(1, round(len(ids) * float(heldout_fraction)))
        target = min(target, len(ids) - 1)
        held_chroms = []
        n_held = 0
        for c in chroms:
            if n_held >= target and held_chroms:
                break
            # Leave at least one train gene.
            if n_held + len(by_chrom[c]) >= len(ids) and held_chroms:
                break
            if n_held + len(by_chrom[c]) >= len(ids):
                continue
            held_chroms.append(c)
            n_held += len(by_chrom[c])
        if not held_chroms:
            # Degenerate: pick the smallest chromosome that leaves ≥1 train gene.
            sized = sorted(chroms, key=lambda c: len(by_chrom[c]))
            for c in sized:
                if len(by_chrom[c]) < len(ids):
                    held_chroms = [c]
                    break
        if not held_chroms:
            raise RuntimeError("chromosome partition could not select a heldout set")

    held_set = {g for c in held_chroms for g in by_chrom[c]}
    heldout = [g for g in ids if g in held_set]
    train = [g for g in ids if g not in held_set]
    if not train or not heldout:
        raise RuntimeError("chromosome partition produced an empty train or heldout set")
    train_chroms = sorted({chrom_of[g] for g in train if g in chrom_of})
    return GeneHoldoutPartition(
        method="chromosome",
        seed=int(seed),
        heldout_fraction=float(heldout_fraction),
        train_gene_ids=train,
        heldout_gene_ids=heldout,
        heldout_chromosomes=list(held_chroms),
        train_chromosomes=train_chroms,
    )


def resolve_gene_holdout_partition(
    gene_ids: Sequence[str],
    *,
    method: str = "random",
    heldout_fraction: float = 0.2,
    seed: int = 42,
    gene_chromosome: Mapping[str, str] | None = None,
    heldout_chromosomes: Sequence[str] | None = None,
    exclude_chromosomes: Sequence[str] | None = None,
) -> GeneHoldoutPartition:
    """Dispatch ``random`` / ``chromosome`` from config-like kwargs."""
    m = str(method).lower().strip()
    if m == "random":
        return partition_genes_random(
            gene_ids, heldout_fraction=heldout_fraction, seed=seed
        )
    if m == "chromosome":
        if gene_chromosome is None:
            raise ValueError("chromosome gene-holdout requires gene_chromosome map")
        return partition_genes_by_chromosome(
            gene_ids,
            gene_chromosome,
            heldout_fraction=heldout_fraction,
            seed=seed,
            heldout_chromosomes=heldout_chromosomes,
            exclude_chromosomes=exclude_chromosomes,
        )
    raise ValueError(
        f"unsupported gene_holdout.method={method!r} (expected 'random' or 'chromosome')"
    )


def subset_flat_region_gene_index(
    index: FlatRegionGeneIndex,
    keep_gene_ids: Sequence[str],
) -> FlatRegionGeneIndex:
    """Keep edges whose gene is in ``keep_gene_ids``; remap gene indices to ``0..K-1``.

    Gene order in the returned index follows ``keep_gene_ids`` (caller order),
    not the original index order. Genes with zero edges after filtering are
    dropped so ``n_genes`` matches columns that can receive MBS mass.
    """
    keep = [str(g) for g in keep_gene_ids]
    if not keep:
        raise ValueError("keep_gene_ids is empty")
    if len(keep) != len(set(keep)):
        raise ValueError("keep_gene_ids must be unique")
    old_ids = list(index.gene_ids)
    old_pos = {g: i for i, g in enumerate(old_ids)}
    missing = [g for g in keep if g not in old_pos]
    if missing:
        raise KeyError(f"{len(missing)} keep genes absent from index (e.g. {missing[0]!r})")

    old_gene = np.asarray(index.edge_gene_index, dtype=np.int64)
    keep_old_indices = np.asarray([old_pos[g] for g in keep], dtype=np.int64)
    edge_mask = np.isin(old_gene, keep_old_indices)
    if not np.any(edge_mask):
        raise ValueError("no edges remain after gene subset")

    present_old = sorted({int(g) for g in old_gene[edge_mask].tolist()})
    present_ids = [old_ids[i] for i in present_old]
    present_set = set(present_ids)
    new_gene_ids = [g for g in keep if g in present_set]
    old_to_new = {old_pos[g]: j for j, g in enumerate(new_gene_ids)}

    def _take(arr: np.ndarray) -> np.ndarray:
        return np.asarray(arr, copy=False)[edge_mask]

    remapped = [old_to_new[int(g)] for g in old_gene[edge_mask].tolist()]
    new_gene = np.asarray(remapped, dtype=np.int64)
    role_ids = _take(index.edge_role_id).astype(np.int64, copy=False)
    other_id = GENE_ROLES.index("other_gene")
    return FlatRegionGeneIndex(
        gene_ids=new_gene_ids,
        edge_col_index=_take(index.edge_col_index).astype(np.int64, copy=False),
        edge_gene_index=new_gene,
        edge_role_id=role_ids,
        edge_context_id=_take(index.edge_context_id).astype(np.int64, copy=False),
        edge_role_present=_take(index.edge_role_present).astype(bool, copy=False),
        edge_context_present=_take(index.edge_context_present).astype(bool, copy=False),
        edge_regulatory_present=_take(index.edge_regulatory_present).astype(bool, copy=False),
        edge_regulatory_multi_hot=np.asarray(
            index.edge_regulatory_multi_hot, dtype=np.float32, copy=False
        )[edge_mask],
        n_study_loci=int(index.n_study_loci),
        n_other_gene_edges=int(np.sum(role_ids == other_id)),
    )


def subset_cascade_assignment(
    assignment: CascadeAssignment,
    keep_gene_ids: Sequence[str],
    *,
    keep_orphans: bool = True,
    keep_direct: bool = True,
) -> CascadeAssignment:
    """Keep gene-linked regions for ``keep_gene_ids``; remap ``region_to_gene``.

    Orphan RBS and direct leftover paths are optional (drop both for a pure
    heldout-gene MBS score matrix used by nested enet).
    """
    keep = [str(g) for g in keep_gene_ids]
    if not keep:
        raise ValueError("keep_gene_ids is empty")
    if len(keep) != len(set(keep)):
        raise ValueError("keep_gene_ids must be unique")
    old_ids = list(assignment.gene_ids)
    old_pos = {g: i for i, g in enumerate(old_ids)}
    missing = [g for g in keep if g not in old_pos]
    if missing:
        raise KeyError(f"{len(missing)} keep genes absent from assignment (e.g. {missing[0]!r})")
    keep_old = {old_pos[g] for g in keep}

    region_keep: list[int] = []
    for i in range(assignment.n_regions):
        g = int(assignment.region_to_gene[i])
        if g < 0:
            if keep_orphans:
                region_keep.append(i)
        elif g in keep_old:
            region_keep.append(i)
    if not region_keep and not (keep_direct and assignment.n_direct):
        raise ValueError("no regions/direct columns remain after gene subset")

    used_regions = np.asarray(region_keep, dtype=np.int64)
    old_reg_to_new = {int(old): j for j, old in enumerate(used_regions.tolist())}

    # Genes that still have ≥1 region, in caller keep order.
    present_old_genes: set[int] = set()
    if used_regions.size:
        for g in assignment.region_to_gene[used_regions].tolist():
            if int(g) >= 0:
                present_old_genes.add(int(g))
    new_gene_ids = [g for g in keep if old_pos[g] in present_old_genes]
    old_gene_to_new = {old_pos[g]: j for j, g in enumerate(new_gene_ids)}

    if used_regions.size:
        new_region_to_gene = np.asarray(
            [
                old_gene_to_new[int(g)] if int(g) >= 0 else ORPHAN_GENE_INDEX
                for g in assignment.region_to_gene[used_regions].tolist()
            ],
            dtype=np.int64,
        )
        new_region_ids = [assignment.region_ids[int(i)] for i in used_regions.tolist()]
        new_region_type_id = assignment.region_type_id[used_regions]
        new_orphan = new_region_to_gene < 0
        new_allocated = [assignment.allocated_gene_id[int(i)] for i in used_regions.tolist()]
        if assignment.edge_region_index.size:
            edge_mask = np.asarray(
                [int(r) in old_reg_to_new for r in assignment.edge_region_index.tolist()],
                dtype=bool,
            )
            edge_col = assignment.edge_col_index[edge_mask]
            edge_reg = np.asarray(
                [old_reg_to_new[int(r)] for r in assignment.edge_region_index[edge_mask].tolist()],
                dtype=np.int64,
            )
        else:
            edge_col = np.zeros(0, dtype=np.int64)
            edge_reg = np.zeros(0, dtype=np.int64)
    else:
        new_region_ids = []
        new_region_type_id = np.zeros(0, dtype=np.int64)
        new_region_to_gene = np.zeros(0, dtype=np.int64)
        new_orphan = np.zeros(0, dtype=bool)
        new_allocated = []
        edge_col = np.zeros(0, dtype=np.int64)
        edge_reg = np.zeros(0, dtype=np.int64)

    if not new_gene_ids:
        raise ValueError("gene subset left zero genes with edges")

    direct = (
        np.asarray(assignment.direct_col_index, dtype=np.int64).copy()
        if keep_direct
        else np.zeros(0, dtype=np.int64)
    )
    return CascadeAssignment(
        gene_ids=new_gene_ids,
        region_ids=new_region_ids,
        region_type_id=np.asarray(new_region_type_id, dtype=np.int64),
        region_to_gene=new_region_to_gene,
        orphan_region_mask=new_orphan,
        edge_col_index=edge_col.astype(np.int64, copy=False),
        edge_region_index=edge_reg,
        direct_col_index=direct,
        region_types=assignment.region_types,
        n_study_loci=int(assignment.n_study_loci),
        allocated_gene_id=new_allocated,
    )


def write_gene_holdout_artifacts(
    score_dir: Path,
    *,
    partition: GeneHoldoutPartition,
    train_gene_ids: Sequence[str],
    heldout_gene_ids: Sequence[str],
    mbs_heldout: np.ndarray,
    present_heldout: np.ndarray,
    extra: dict[str, Any] | None = None,
) -> Path:
    """Persist holdout partition + heldout MBS next to Stage A / cascade score artifacts."""
    score_dir.mkdir(parents=True, exist_ok=True)
    np.save(score_dir / "mbs_heldout.npy", np.asarray(mbs_heldout, dtype=np.float32))
    np.save(
        score_dir / "mbs_heldout_present.npy",
        np.asarray(present_heldout, dtype=bool).astype(np.uint8),
    )
    (score_dir / "gene_ids_train.json").write_text(
        json.dumps(list(train_gene_ids), indent=2) + "\n", encoding="utf-8"
    )
    (score_dir / "gene_ids_heldout.json").write_text(
        json.dumps(list(heldout_gene_ids), indent=2) + "\n", encoding="utf-8"
    )
    payload: dict[str, Any] = {
        **partition.to_dict(),
        "n_train_genes_scored": len(train_gene_ids),
        "n_heldout_genes_scored": len(heldout_gene_ids),
        "mbs_heldout_shape": list(np.asarray(mbs_heldout).shape),
    }
    if extra:
        payload.update(extra)
    path = score_dir / "gene_holdout.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path
