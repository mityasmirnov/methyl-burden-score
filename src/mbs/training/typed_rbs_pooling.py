"""CPU typed-RBS pooling: presence-aware gene/role feature builders + eval helpers.

Milestone 7G amendment. Given the ``all_gene_rbs`` region matrix (``[n_samples,
n_regions]``) plus a matching ``all_gene_rbs_present`` mask and a region index
(``region_id, gene_id, region_type, column_index``), build a family of pooling
arms (R0..R5 + a role-shuffle control) and score them with the transparent
linear multitask heads from :mod:`mbs.training.transparent_baselines`.

Per (gene, role) the base feature is presence-aware and centered::

    x_{g,r} = (RBS_{g,r} - 0.5) * present_{g,r}

Pooling is over *present* regions only (neutral 0 when a role is absent), so
missingness never reads as low burden (scientific invariant 5).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import partial
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.gencode_regions import REGION_TYPES
from mbs.training.transparent_baselines import run_mean_baseline

# R4 collapses the five roles into promoter- vs body-side channels.
R4_ROLE_MAP: dict[str, str] = {
    "promoter_core": "promoter",
    "promoter_proximal": "promoter",
    "five_prime": "promoter",
    "gene_body": "body",
    "three_prime": "body",
}
R4_ROLES: tuple[str, ...] = ("promoter", "body")
UNTYPED_ROLE = "gene"


@dataclass(frozen=True)
class GeneRoleLayout:
    """Maps region-matrix columns onto ``(gene, role)`` pooling groups."""

    genes: tuple[str, ...]
    role_names: tuple[str, ...]
    # groups[(gene_idx, role_idx)] -> column positions in the region matrix.
    groups: Mapping[tuple[int, int], np.ndarray]

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_roles(self) -> int:
        return len(self.role_names)


def _ordered_region_index(region_df: pd.DataFrame) -> pd.DataFrame:
    """Return ``region_df`` ordered so row ``j`` is region-matrix column ``j``."""
    required = {"gene_id", "region_type"}
    missing = required - set(region_df.columns)
    if missing:
        raise ValueError(f"region index missing columns: {sorted(missing)}")
    if "column_index" in region_df.columns:
        ordered = region_df.sort_values("column_index", kind="mergesort").reset_index(drop=True)
        expected = np.arange(len(ordered), dtype=np.int64)
        if not np.array_equal(ordered["column_index"].to_numpy(dtype=np.int64), expected):
            raise ValueError("column_index must be a dense 0..n-1 range")
        return ordered
    return region_df.reset_index(drop=True)


def build_layout(
    region_df: pd.DataFrame,
    *,
    role_of: Callable[[str], str | None] | None = None,
    role_names: Sequence[str] = REGION_TYPES,
) -> GeneRoleLayout:
    """Build a pooling layout from a region index.

    ``role_of`` maps a raw ``region_type`` to a pooling-role name (or ``None`` to
    drop the region). Defaults to the identity over ``role_names``.
    """
    ordered = _ordered_region_index(region_df)
    names = tuple(role_names)
    role_to_idx = {r: i for i, r in enumerate(names)}

    def _identity_role(region_type: str) -> str | None:
        return region_type if region_type in role_to_idx else None

    resolve_role = _identity_role if role_of is None else role_of

    gene_ids = ordered["gene_id"].astype(str).to_numpy()
    region_types = ordered["region_type"].astype(str).to_numpy()
    gene_order = tuple(dict.fromkeys(gene_ids.tolist()))
    gene_to_idx = {g: i for i, g in enumerate(gene_order)}

    groups: dict[tuple[int, int], list[int]] = {}
    for col, (gene, rtype) in enumerate(zip(gene_ids, region_types, strict=True)):
        role = resolve_role(str(rtype))
        if role is None:
            continue
        key = (gene_to_idx[str(gene)], role_to_idx[role])
        groups.setdefault(key, []).append(col)

    frozen = {k: np.asarray(v, dtype=np.int64) for k, v in groups.items()}
    return GeneRoleLayout(genes=gene_order, role_names=names, groups=frozen)


def pool_gene_role(
    rbs: np.ndarray,
    present: np.ndarray,
    layout: GeneRoleLayout,
) -> dict[str, np.ndarray]:
    """Presence-aware pooling per ``(gene, role)``.

    Returns four ``[n_samples, n_genes, n_roles]`` arrays: ``max`` and ``mean``
    of centered RBS over *present* regions (0 when none present), a ``present``
    flag, and the ``count`` of contributing regions.
    """
    vals = np.asarray(rbs, dtype=np.float64)
    obs = np.asarray(present).astype(bool)
    if vals.shape != obs.shape:
        raise ValueError("rbs and present shape mismatch")
    centered = vals - 0.5
    ns = vals.shape[0]
    ng, nr = layout.n_genes, layout.n_roles
    mx = np.zeros((ns, ng, nr), dtype=np.float32)
    mn = np.zeros((ns, ng, nr), dtype=np.float32)
    pf = np.zeros((ns, ng, nr), dtype=np.float32)
    cnt = np.zeros((ns, ng, nr), dtype=np.float32)
    for (gi, ri), cols in layout.groups.items():
        block = centered[:, cols]
        mask = obs[:, cols]
        k = mask.sum(axis=1)
        has = k > 0
        cnt[:, gi, ri] = k
        pf[has, gi, ri] = 1.0
        sums = np.where(mask, block, 0.0).sum(axis=1)
        mn[has, gi, ri] = sums[has] / k[has]
        maxed = np.where(mask, block, -np.inf).max(axis=1)
        mx[has, gi, ri] = maxed[has]
    return {"max": mx, "mean": mn, "present": pf, "count": cnt}


def _flatten(
    pooled: dict[str, np.ndarray],
    layout: GeneRoleLayout,
    keys: Sequence[str],
) -> tuple[np.ndarray, list[str]]:
    ns = next(iter(pooled.values())).shape[0]
    blocks: list[np.ndarray] = []
    names: list[str] = []
    for key in keys:
        arr = pooled[key].reshape(ns, layout.n_genes * layout.n_roles)
        blocks.append(arr)
        for gi in range(layout.n_genes):
            for ri in range(layout.n_roles):
                names.append(f"{layout.genes[gi]}:{layout.role_names[ri]}:{key}")
    return np.concatenate(blocks, axis=1).astype(np.float32, copy=False), names


# --- Arm feature builders -------------------------------------------------


def features_untyped(
    rbs: np.ndarray, present: np.ndarray, region_df: pd.DataFrame
) -> tuple[np.ndarray, list[str]]:
    """R0: untyped pooling — one channel per gene (max + mean + present + count)."""
    layout = build_layout(region_df, role_of=lambda _t: UNTYPED_ROLE, role_names=(UNTYPED_ROLE,))
    pooled = pool_gene_role(rbs, present, layout)
    return _flatten(pooled, layout, ("max", "mean", "present", "count"))


def features_typed(
    rbs: np.ndarray,
    present: np.ndarray,
    region_df: pd.DataFrame,
    *,
    stats: Sequence[str],
    role_of: Callable[[str], str | None] | None = None,
    role_names: Sequence[str] = REGION_TYPES,
) -> tuple[np.ndarray, list[str]]:
    """Typed pooling with the requested ``stats`` plus present flags + counts."""
    layout = build_layout(region_df, role_of=role_of, role_names=role_names)
    pooled = pool_gene_role(rbs, present, layout)
    keys = [*stats, "present", "count"]
    return _flatten(pooled, layout, keys)


def features_passthrough(
    rbs: np.ndarray, present: np.ndarray, region_df: pd.DataFrame
) -> tuple[np.ndarray, list[str]]:
    """R5 ceiling: presence-aware region matrix, no pooling."""
    ordered = _ordered_region_index(region_df)
    vals = np.asarray(rbs, dtype=np.float64) - 0.5
    obs = np.asarray(present).astype(bool)
    feats = np.where(obs, vals, 0.0).astype(np.float32, copy=False)
    names = [f"{g}:{t}" for g, t in zip(ordered["region_id"], ordered["region_type"], strict=False)]
    if len(names) != feats.shape[1]:
        names = [f"col_{j}" for j in range(feats.shape[1])]
    return feats, names


ArmBuilder = Callable[[np.ndarray, np.ndarray, pd.DataFrame], tuple[np.ndarray, list[str]]]


def arm_builders(*, include_r5: bool = False) -> dict[str, ArmBuilder]:
    """Return the ``arm_name -> builder`` map (R5 optional, it is the slow ceiling)."""
    builders: dict[str, ArmBuilder] = {
        "R0": features_untyped,
        "R1": partial(features_typed, stats=("max",)),
        "R2": partial(features_typed, stats=("mean",)),
        "R3": partial(features_typed, stats=("max", "mean")),
        "R4": partial(
            features_typed,
            stats=("max", "mean"),
            role_of=R4_ROLE_MAP.get,
            role_names=R4_ROLES,
        ),
    }
    if include_r5:
        builders["R5"] = features_passthrough
    return builders


# --- Role-shuffle control -------------------------------------------------


def shuffle_region_types(
    region_df: pd.DataFrame, *, seed: int
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Permute ``region_type`` within each gene; report how much actually changed.

    Returns a new region index and stats: ``frac_genes_altered`` (genes with at
    least one region whose role changed) and ``frac_columns_changed`` (fraction
    of region columns whose assigned role differs from the original).
    """
    ordered = _ordered_region_index(region_df).copy()
    rng = np.random.default_rng(seed)
    original = ordered["region_type"].astype(str).to_numpy()
    shuffled = original.copy()
    genes_altered = 0
    total_genes = 0
    for idx in ordered.groupby("gene_id", sort=False).groups.values():
        pos = np.asarray(idx, dtype=np.int64)
        total_genes += 1
        if pos.size < 2:
            continue
        perm = rng.permutation(pos.size)
        shuffled[pos] = original[pos][perm]
        if np.any(shuffled[pos] != original[pos]):
            genes_altered += 1
    ordered["region_type"] = shuffled
    stats = {
        "seed": float(seed),
        "frac_genes_altered": genes_altered / total_genes if total_genes else 0.0,
        "frac_columns_changed": float(np.mean(shuffled != original)) if original.size else 0.0,
    }
    return ordered, stats


# --- Normalization (train-fold only) --------------------------------------


def fit_standardizer(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column mean/std from the train fold; zero-variance columns get std 1."""
    arr = np.asarray(x, dtype=np.float64)
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    std = np.where(std < 1e-8, 1.0, std)
    return mean.astype(np.float64), std.astype(np.float64)


def apply_standardizer(x: np.ndarray, stats: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    mean, std = stats
    return ((np.asarray(x, dtype=np.float64) - mean) / std).astype(np.float32, copy=False)


# --- Arm evaluation -------------------------------------------------------


def evaluate_arm(
    *,
    name: str,
    x_train: np.ndarray,
    x_test: np.ndarray,
    pheno_train: Mapping[str, np.ndarray],
    pheno_test: Mapping[str, np.ndarray],
    class_names: list[str] | None = None,
    standardize: bool = True,
) -> dict[str, Any]:
    """Fold-fit linear multitask heads on one arm's features and score holdout."""
    x_tr = np.asarray(x_train, dtype=np.float32)
    x_te = np.asarray(x_test, dtype=np.float32)
    if standardize:
        stats = fit_standardizer(x_tr)
        x_tr = apply_standardizer(x_tr, stats)
        x_te = apply_standardizer(x_te, stats)
    result = run_mean_baseline(
        x_train=x_tr,
        x_test=x_te,
        age_train=pheno_train.get("age"),
        age_mask_train=pheno_train.get("age_mask"),
        tissue_train=pheno_train.get("tissue"),
        tissue_mask_train=pheno_train.get("tissue_mask"),
        sex_train=pheno_train.get("sex"),
        sex_mask_train=pheno_train.get("sex_mask"),
        age_test=pheno_test.get("age"),
        age_mask_test=pheno_test.get("age_mask"),
        tissue_test=pheno_test.get("tissue"),
        tissue_mask_test=pheno_test.get("tissue_mask"),
        sex_test=pheno_test.get("sex"),
        sex_mask_test=pheno_test.get("sex_mask"),
        study_ids_test=pheno_test.get("study_ids"),
        tissue_class_names=class_names,
    )
    result["arm"] = name
    return result


# --- Promotion gate -------------------------------------------------------


def typed_pool_promotion_gate(
    *,
    age_mae_r0: float,
    age_mae_typed: float,
    age_r2_r0: float,
    age_r2_typed: float,
    tissue_f1_r0: float,
    tissue_f1_typed: float,
    sex_auroc_r0: float | None,
    sex_auroc_typed: float | None,
    shuffle_age_mae: float,
) -> dict[str, Any]:
    """Decide whether typed pooling beats the untyped R0 baseline.

    Promote when: age MAE improves by >=1 year OR age R^2 improves by >=0.05;
    AND tissue macro-F1 loss <=0.03; AND (if both present) sex AUROC loss
    <=0.03; AND the role-shuffle control is clearly worse than typed on age MAE.
    """
    age_mae_gain = float(age_mae_r0) - float(age_mae_typed)
    age_r2_gain = float(age_r2_typed) - float(age_r2_r0)
    tissue_f1_loss = float(tissue_f1_r0) - float(tissue_f1_typed)
    age_improves = age_mae_gain >= 1.0 or age_r2_gain >= 0.05
    tissue_ok = tissue_f1_loss <= 0.03
    if sex_auroc_r0 is not None and sex_auroc_typed is not None:
        sex_loss: float | None = float(sex_auroc_r0) - float(sex_auroc_typed)
        sex_ok = sex_loss <= 0.03
    else:
        sex_loss = None
        sex_ok = True
    # Require a material collapse under role permutation (not a 0.1 y wobble).
    shuffle_delta = float(shuffle_age_mae) - float(age_mae_typed)
    shuffle_worse = shuffle_delta >= 1.0
    promote = bool(age_improves and tissue_ok and sex_ok and shuffle_worse)
    return {
        "promote": promote,
        "age_mae_gain": age_mae_gain,
        "age_r2_gain": age_r2_gain,
        "tissue_f1_loss": tissue_f1_loss,
        "sex_auroc_loss": sex_loss,
        "age_improves": bool(age_improves),
        "tissue_ok": bool(tissue_ok),
        "sex_ok": bool(sex_ok),
        "shuffle_worse_than_typed": bool(shuffle_worse),
        "shuffle_age_mae_delta": shuffle_delta,
    }
