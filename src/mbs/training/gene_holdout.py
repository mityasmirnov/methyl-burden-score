"""Gene-holdout partition for DeepRVAT-style encoder transfer (GATE G1).

Train φ/ρ on one gene set, score a disjoint heldout set with the frozen
encoder, then fit the post-hoc nested elastic-net on heldout MBS columns the
encoder never saw during training.

Methods:
- ``random`` — per-gene shuffle (locality leakage possible across the boundary).
- ``chromosome`` — hold out whole chromosomes (bounds neighbour leakage).
- ``seed`` — DeepRVAT-aligned: train on fold seed-gene union (CpG-first
  best-association → gene), score the complement. Not 9c head masks.
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

GeneHoldoutMethod = Literal["random", "chromosome", "seed"]
TraitTask = Literal["age", "tissue", "sex"]

_RANDOM_CAVEAT = (
    "random split: neighbouring/co-regulated genes may straddle the "
    "boundary (locality leakage). Chromosome-held-out arm bounds that."
)
_CHROM_CAVEAT = (
    "chromosome holdout: entire chromosomes are held out so neighbouring "
    "genes cannot straddle the train/heldout boundary. Sex chromosomes may "
    "be included unless excluded via config."
)
_SEED_CAVEAT = (
    "seed holdout (DeepRVAT-aligned): phi/rho train on phenotype-informed seed "
    "genes (CpG-first best association -> gene; multi-trait union); score the "
    "complement. Discovery uses outer-train labels only. Not a head mask."
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
    discovery: dict[str, Any] = field(default_factory=dict)

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
        if self.method == "chromosome":
            return _CHROM_CAVEAT
        if self.method == "seed":
            return _SEED_CAVEAT
        return _RANDOM_CAVEAT

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["actual_heldout_fraction"] = self.actual_heldout_fraction
        out["caveat"] = self.caveat()
        return out


def cpg_association_strength(
    x: np.ndarray,
    y: np.ndarray,
    *,
    task: TraitTask,
) -> np.ndarray:
    """Per-column association strength (higher = stronger). Outer-train only.

    Age: absolute Pearson-like correlation. Tissue/sex: max absolute class-mean
    difference. Same spirit as ``fold_safe_panel._univariate_prefilter`` scores,
    but returns the full column vector for CpG→gene aggregation.
    """
    x64 = np.asarray(x, dtype=np.float64)
    if x64.ndim != 2:
        raise ValueError(f"x must be 2-D; got shape {x64.shape}")
    n_cols = int(x64.shape[1])
    if x64.shape[0] < 2:
        return np.zeros(n_cols, dtype=np.float64)
    col_mean = np.nanmean(x64, axis=0)
    filled = np.where(np.isfinite(x64), x64, col_mean)
    if task == "age":
        y64 = np.asarray(y, dtype=np.float64)
        y_c = y64 - float(np.mean(y64))
        x_c = filled - filled.mean(axis=0, keepdims=True)
        denom = np.sqrt((x_c * x_c).sum(axis=0) * float((y_c * y_c).sum())) + 1e-12
        return np.abs((x_c * y_c[:, None]).sum(axis=0) / denom)
    y_i = np.asarray(y).astype(np.int64, copy=False)
    classes = np.unique(y_i)
    score = np.zeros(n_cols, dtype=np.float64)
    for c in classes:
        mask_c = y_i == c
        if not mask_c.any() or bool(mask_c.all()):
            continue
        diff = filled[mask_c].mean(axis=0) - filled[~mask_c].mean(axis=0)
        score = np.maximum(score, np.abs(diff))
    return score


def discover_seed_genes_best_cpg(
    gene_ids: Sequence[str],
    *,
    edge_col_index: np.ndarray,
    edge_gene_index: np.ndarray,
    x_train: np.ndarray,
    trait_labels: Mapping[str, tuple[np.ndarray, np.ndarray]],
    n_genes_per_trait: int = 256,
    extra_seed_gene_ids: Sequence[str] | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """CpG-first seed discovery: best linked-CpG association → gene; union traits.

    For each trait, score every CpG on labeled outer-train rows, assign each gene
    the **max** linked CpG strength (strongest association / lowest-p analogue),
    keep the top ``n_genes_per_trait``, then union across traits. Optional
    ``extra_seed_gene_ids`` (e.g. Atlas ``external_clean``) are merged in.
    """
    ids = [str(g) for g in gene_ids]
    if not ids:
        raise ValueError("gene_ids is empty")
    n_genes = len(ids)
    cols = np.asarray(edge_col_index, dtype=np.int64)
    genes = np.asarray(edge_gene_index, dtype=np.int64)
    if cols.size == 0 or genes.size == 0:
        raise ValueError("no edges for seed discovery")
    if cols.shape != genes.shape:
        raise ValueError("edge_col_index and edge_gene_index shape mismatch")
    x = np.asarray(x_train)
    if x.ndim != 2:
        raise ValueError(f"x_train must be 2-D; got {x.shape}")
    n_cols = int(x.shape[1])
    if int(cols.max()) >= n_cols:
        raise ValueError(
            f"edge col {int(cols.max())} out of range for x_train ncols={n_cols}"
        )
    if int(genes.max()) >= n_genes or int(genes.min()) < 0:
        raise ValueError("edge_gene_index out of range for gene_ids")

    per_trait: dict[str, Any] = {}
    union: set[str] = set()
    for trait, (y_raw, mask_raw) in trait_labels.items():
        task = str(trait).lower().strip()
        if task not in {"age", "tissue", "sex"}:
            raise ValueError(
                f"unsupported seed trait {trait!r}; expected age/tissue/sex"
            )
        mask = np.asarray(mask_raw, dtype=bool)
        if mask.shape[0] != x.shape[0]:
            raise ValueError(
                f"trait {trait!r} mask length {mask.shape[0]} != n_train {x.shape[0]}"
            )
        if int(mask.sum()) < 4:
            per_trait[task] = {"n_labeled": int(mask.sum()), "skipped": True}
            continue
        y = np.asarray(y_raw)[mask]
        col_scores = cpg_association_strength(x[mask], y, task=task)  # type: ignore[arg-type]
        gene_best = np.full(n_genes, -np.inf, dtype=np.float64)
        np.maximum.at(gene_best, genes, col_scores[cols])
        finite = np.isfinite(gene_best)
        if not finite.any():
            per_trait[task] = {"n_labeled": int(mask.sum()), "n_seed_genes": 0}
            continue
        order = np.argsort(-gene_best, kind="stable")
        k = max(1, int(n_genes_per_trait))
        picked: list[str] = []
        strengths: list[float] = []
        for ix in order.tolist():
            if not finite[ix]:
                continue
            picked.append(ids[int(ix)])
            strengths.append(float(gene_best[int(ix)]))
            if len(picked) >= k:
                break
        union.update(picked)
        per_trait[task] = {
            "n_labeled": int(mask.sum()),
            "n_seed_genes": len(picked),
            "top_gene_ids": picked[:16],
            "top_strength": strengths[:8],
        }

    if extra_seed_gene_ids:
        panel = set(ids)
        extra = [str(g) for g in extra_seed_gene_ids if str(g) in panel]
        union.update(extra)
        n_extra = len(extra)
    else:
        n_extra = 0

    if not union:
        raise RuntimeError("seed discovery produced an empty gene union")
    # Stable order: panel order for reproducibility.
    seed_list = [g for g in ids if g in union]
    meta: dict[str, Any] = {
        "discovery": "best_cpg_p",
        "n_genes_per_trait": int(n_genes_per_trait),
        "traits": per_trait,
        "n_seed_genes": len(seed_list),
        "n_extra_seed_genes": n_extra,
    }
    return seed_list, meta


def partition_genes_by_seed(
    gene_ids: Sequence[str],
    seed_gene_ids: Sequence[str],
    *,
    seed: int = 42,
    discovery: Mapping[str, Any] | None = None,
) -> GeneHoldoutPartition:
    """Train = seed ∩ panel; heldout = panel \\ seed (DeepRVAT-aligned)."""
    ids = [str(g) for g in gene_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("gene_ids must be unique")
    if len(ids) < 2:
        raise ValueError("need at least 2 genes to partition")
    seed_set = {str(g) for g in seed_gene_ids}
    train = [g for g in ids if g in seed_set]
    heldout = [g for g in ids if g not in seed_set]
    if not train:
        raise RuntimeError("seed partition: no seed genes overlap the panel")
    if not heldout:
        raise RuntimeError(
            "seed partition: seed union covers the whole panel; nothing to score"
        )
    return GeneHoldoutPartition(
        method="seed",
        seed=int(seed),
        heldout_fraction=float(len(heldout)) / float(len(ids)),
        train_gene_ids=train,
        heldout_gene_ids=heldout,
        discovery=dict(discovery or {}),
    )

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
    seed_gene_ids: Sequence[str] | None = None,
    discovery: Mapping[str, Any] | None = None,
) -> GeneHoldoutPartition:
    """Dispatch ``random`` / ``chromosome`` / ``seed`` from config-like kwargs."""
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
    if m == "seed":
        if seed_gene_ids is None:
            raise ValueError(
                "seed gene-holdout requires seed_gene_ids "
                "(run discover_seed_genes_best_cpg on outer-train first)"
            )
        return partition_genes_by_seed(
            gene_ids,
            seed_gene_ids,
            seed=seed,
            discovery=discovery,
        )
    raise ValueError(
        "unsupported gene_holdout.method="
        f"{method!r} (expected 'random', 'chromosome', or 'seed')"
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
