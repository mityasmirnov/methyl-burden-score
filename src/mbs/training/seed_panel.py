"""``internal_fold`` seed-gene panel constructor (7G′ Stage B, ADR 0010/0011/0012).

Gene-first, fold-fitted seed selection: reuse ``stability_select_columns`` on
outer-train samples only, map selected CpG columns to genes via *explicit*
region→gene edges, rank genes by capped burden strength, then enrich each
selected gene with all its gene-linked CpGs (siblings). DeepRVAT-faithful:
selection is a burden over the gene's CpG set, not a single min-P locus.

Discovery CpGs (prefilter / stability) rank genes only — they are not the G2/C2
input panel (ADR 0012). Callers should pass an ``explicit_only``
:class:`CascadeAssignment` (ADR 0010) so neural and classical arms share
evidence-backed gene columns.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from mbs.annotation.manifest import sha256_file
from mbs.training.cascade_assign import CascadeAssignment
from mbs.training.fold_safe_panel import stability_select_columns

TraitRole = Literal["primary", "secondary", "auxiliary"]
_PROMOTER_ROLES = frozenset({"promoter_core", "promoter_proximal", "five_prime"})
_BODY_ROLES = frozenset({"gene_body", "three_prime"})
_SEX_CHROMS = frozenset({"chrx", "chry", "x", "y"})
_GRAPH_HASH_FILES = (
    "genes.parquet",
    "regions.parquet",
    "locus_region_edges.parquet",
    "region_gene_edges.parquet",
)
DEFAULT_STRENGTH_CAP_QUANTILE = 0.99
DEFAULT_MIN_FREQUENCY = 0.34
# ponytail: full 51k-col × 5×3×9 enet grid is multi-hour; univariate prefilter
# + lighter CV keeps DeepRVAT-faithful fold-fitting while finishing a screen.
DEFAULT_PREFILTER_MAX_COLS = 4096
DEFAULT_N_INNER_FOLDS = 2
DEFAULT_N_REPEATS = 2
SELECTION_METHOD = "study_grouped_enet_stability_gene_first"

# Traits with label arrays wired into :func:`build_internal_fold_seed_panel`.
_SUPPORTED_PANEL_TRAITS = frozenset({"age", "tissue", "sex"})
# Future traits: config-driven, but not yet labeled in this constructor.
_BLOCKED_TRAIT_HINTS: dict[str, str] = {
    "bmi": (
        "BMI is blocked on ATS; eligible only on matrix-hub-bmi-full-v1 "
        "(≥1000 samples, ≥5 studies). Do not join BMI onto the ATS matrix."
    ),
    "disease": (
        "disease requires documented multi-study cases+controls; "
        "unknown labels must remain unknown (not controls)"
    ),
    "cancer": (
        "cancer requires documented multi-study cases+controls; "
        "unknown labels must remain unknown (not controls)"
    ),
    "blood": "blood subtraits blocked until ontology / label quality confirmed",
    "brain": "brain subtraits blocked until ontology / label quality confirmed",
}


@dataclass(frozen=True, slots=True)
class SeedTraitSpec:
    """One seed-panel trait from experiment YAML (not Hub phenotype_registry)."""

    id: str
    role: TraitRole = "primary"
    autosome_control: bool = False


DEFAULT_SEED_TRAITS: tuple[SeedTraitSpec, ...] = (
    SeedTraitSpec("age", "primary"),
    SeedTraitSpec("tissue", "secondary"),
    SeedTraitSpec("sex", "auxiliary", autosome_control=True),
)

# Back-compat alias for callers / tests that still import TRAITS.
TRAITS: tuple[str, ...] = tuple(t.id for t in DEFAULT_SEED_TRAITS)


def resolve_seed_panel_traits(
    traits: Sequence[SeedTraitSpec | Mapping[str, Any] | str] | None = None,
) -> list[SeedTraitSpec]:
    """Parse experiment ``seed_panel.traits`` into :class:`SeedTraitSpec` list.

    ``None`` / empty → ATS default (age / tissue / sex + sex autosome control).
    Unknown or blocked ids raise with an eligibility hint.
    """
    if not traits:
        return list(DEFAULT_SEED_TRAITS)
    out: list[SeedTraitSpec] = []
    seen: set[str] = set()
    for raw in traits:
        if isinstance(raw, SeedTraitSpec):
            spec = raw
        elif isinstance(raw, str):
            spec = SeedTraitSpec(id=raw)
        elif isinstance(raw, Mapping):
            tid = str(raw.get("id") or raw.get("trait") or "").strip()
            if not tid:
                raise ValueError(f"seed_panel trait entry missing id: {raw!r}")
            role_raw = str(raw.get("role") or "primary").strip().lower()
            if role_raw not in {"primary", "secondary", "auxiliary"}:
                raise ValueError(
                    f"seed_panel trait {tid!r} has invalid role {role_raw!r}; "
                    "expected primary|secondary|auxiliary"
                )
            spec = SeedTraitSpec(
                id=tid,
                role=role_raw,  # type: ignore[arg-type]
                autosome_control=bool(raw.get("autosome_control", False)),
            )
        else:
            raise TypeError(f"unsupported seed_panel trait entry: {type(raw)!r}")
        tid = spec.id.strip().lower()
        if tid in seen:
            raise ValueError(f"duplicate seed_panel trait {tid!r}")
        seen.add(tid)
        if tid in _BLOCKED_TRAIT_HINTS:
            raise ValueError(f"seed_panel trait {tid!r}: {_BLOCKED_TRAIT_HINTS[tid]}")
        if tid not in _SUPPORTED_PANEL_TRAITS:
            raise ValueError(
                f"unknown seed_panel trait {tid!r}; supported={sorted(_SUPPORTED_PANEL_TRAITS)}; "
                f"blocked={sorted(_BLOCKED_TRAIT_HINTS)}"
            )
        if spec.autosome_control and tid != "sex":
            raise ValueError(
                f"autosome_control is only valid for sex; got trait={tid!r}"
            )
        out.append(
            SeedTraitSpec(id=tid, role=spec.role, autosome_control=spec.autosome_control)
        )
    return out


def hash_graph_tables(graph_root: Path) -> str:
    """Combined SHA-256 of the graph tables used for seed-panel construction."""
    root = Path(graph_root)
    digest = hashlib.sha256()
    for name in _GRAPH_HASH_FILES:
        path = root / name
        if not path.is_file():
            raise FileNotFoundError(f"graph table missing for content hash: {path}")
        digest.update(name.encode("utf-8"))
        digest.update(b"\x00")
        digest.update(sha256_file(path).encode("utf-8"))
        digest.update(b"\x00")
    return digest.hexdigest()


def _is_sex_chrom(chrom: str) -> bool:
    c = str(chrom).strip().lower()
    if c.startswith("chr"):
        c = c[3:]
    return c in {"x", "y"}


def _autosomal_col_mask(
    locus_chrom: np.ndarray | list[str] | None, n_cols: int
) -> np.ndarray:
    """Boolean mask over local columns; True = autosomal.

    Raises if ``locus_chrom`` is missing — sex_autosome control must not silently
    treat unknown chroms as autosomal.
    """
    if locus_chrom is None:
        raise ValueError(
            "locus_chrom is required for autosomal filtering / sex_autosome control"
        )
    chroms = np.asarray(locus_chrom, dtype=object)
    if chroms.shape[0] != n_cols:
        raise ValueError(f"locus_chrom length {chroms.shape[0]} != n_cols {n_cols}")
    return np.asarray([not _is_sex_chrom(str(c)) for c in chroms.tolist()], dtype=bool)


@dataclass(frozen=True, slots=True)
class SeedPanelArtifacts:
    """Materialised seed panel: JSON summary + gene / locus tables + content hash."""

    panel_json: dict[str, Any]
    genes: pd.DataFrame  # seed_panel_gene
    loci: pd.DataFrame  # seed_panel_locus
    panel_hash: str


def _gene_edge_maps(
    assignment: CascadeAssignment,
) -> tuple[dict[int, set[int]], dict[int, set[int]], dict[tuple[int, int], str]]:
    """Return (col→genes, gene→cols, (col,gene)→role) for explicit gene edges only."""
    col_to_genes: dict[int, set[int]] = {}
    gene_to_cols: dict[int, set[int]] = {}
    edge_role: dict[tuple[int, int], str] = {}
    if assignment.edge_col_index.size == 0:
        return col_to_genes, gene_to_cols, edge_role
    region_to_gene = assignment.region_to_gene
    region_type_id = assignment.region_type_id
    region_types = assignment.region_types
    for col, reg in zip(
        assignment.edge_col_index.tolist(),
        assignment.edge_region_index.tolist(),
        strict=True,
    ):
        gene = int(region_to_gene[int(reg)])
        if gene < 0:
            continue
        col_i = int(col)
        col_to_genes.setdefault(col_i, set()).add(gene)
        gene_to_cols.setdefault(gene, set()).add(col_i)
        key = (col_i, gene)
        if key not in edge_role:
            tid = int(region_type_id[int(reg)])
            edge_role[key] = region_types[tid] if 0 <= tid < len(region_types) else "unknown"
    return col_to_genes, gene_to_cols, edge_role


def _tissue_multi_study_mask(
    tissue: np.ndarray,
    mask: np.ndarray,
    study_ids: np.ndarray,
) -> np.ndarray:
    """Keep only tissue classes observed in >= 2 studies (drop single-study markers)."""
    keep = np.zeros(mask.shape[0], dtype=bool)
    studies_per_class: dict[int, set[str]] = defaultdict(set)
    labeled = np.flatnonzero(np.asarray(mask, dtype=bool))
    for i in labeled.tolist():
        studies_per_class[int(tissue[i])].add(str(study_ids[i]))
    multi = {cls for cls, studies in studies_per_class.items() if len(studies) >= 2}
    for i in labeled.tolist():
        if int(tissue[i]) in multi:
            keep[i] = True
    return keep


def _promoter_body_coverage(roles: set[str]) -> str | None:
    if not roles:
        return None
    has_promoter = bool(roles & _PROMOTER_ROLES)
    has_body = bool(roles & _BODY_ROLES)
    if has_promoter and has_body:
        return "both"
    if has_promoter:
        return "promoter"
    if has_body:
        return "body"
    return "other"


def _gene_chrom_is_autosomal(
    cols: set[int],
    locus_chrom: np.ndarray | list[str] | None,
) -> bool | None:
    """None if unknown; else False when any linked CpG sits on a sex chromosome."""
    if locus_chrom is None:
        return None
    for col in cols:
        if 0 <= col < len(locus_chrom):
            chrom = str(locus_chrom[col]).strip().lower()
            if chrom in _SEX_CHROMS:
                return False
    return True


def _univariate_prefilter(
    x: np.ndarray,
    y: np.ndarray,
    *,
    task: Literal["age", "sex", "tissue"],
    max_cols: int,
) -> np.ndarray:
    """Keep the top ``max_cols`` by univariate association (outer-train only)."""
    n_cols = int(x.shape[1])
    if n_cols <= max_cols:
        return np.arange(n_cols, dtype=np.int64)
    x64 = np.asarray(x, dtype=np.float64)
    # Column-wise nan-safe centering for correlation / mean diffs.
    col_mean = np.nanmean(x64, axis=0)
    filled = np.where(np.isfinite(x64), x64, col_mean)
    if task == "age":
        y64 = np.asarray(y, dtype=np.float64)
        y_c = y64 - float(np.mean(y64))
        x_c = filled - filled.mean(axis=0, keepdims=True)
        denom = np.sqrt((x_c * x_c).sum(axis=0) * float((y_c * y_c).sum())) + 1e-12
        score = np.abs((x_c * y_c[:, None]).sum(axis=0) / denom)
    else:
        y_i = np.asarray(y).astype(np.int64, copy=False)
        classes = np.unique(y_i)
        score = np.zeros(n_cols, dtype=np.float64)
        for c in classes:
            mask_c = y_i == c
            if not mask_c.any() or mask_c.all():
                continue
            diff = filled[mask_c].mean(axis=0) - filled[~mask_c].mean(axis=0)
            score = np.maximum(score, np.abs(diff))
    order = np.argsort(-score, kind="stable")
    return order[:max_cols].astype(np.int64)


def _select_trait_genes(
    *,
    trait: str,
    y: np.ndarray,
    mask: np.ndarray,
    x_train: np.ndarray,
    study_ids: np.ndarray,
    col_to_genes: dict[int, set[int]],
    seed: int,
    strength_cap_quantile: float,
    min_frequency: float,
    prefilter_max_cols: int = DEFAULT_PREFILTER_MAX_COLS,
    n_inner_folds: int = DEFAULT_N_INNER_FOLDS,
    n_repeats: int = DEFAULT_N_REPEATS,
    allowed_cols: np.ndarray | None = None,
) -> dict[str, Any] | None:
    """Run stability selection for one trait; score + rank supported genes.

    ``allowed_cols`` is an optional int index array into ``x_train`` columns
    (e.g. autosomal-only). Selection runs only on those columns; returned
    ``seed_cols`` are indices into the full ``x_train``.
    """
    trait_mask = np.asarray(mask, dtype=bool)
    if trait == "tissue":
        trait_mask = trait_mask & _tissue_multi_study_mask(y, trait_mask, study_ids)
    if not trait_mask.any():
        return None
    task: Literal["age", "sex", "tissue"] = (
        "age" if trait in {"age"} else "sex" if trait in {"sex", "sex_autosome"} else "tissue"
    )
    x_full = x_train[trait_mask]
    y_labeled = np.asarray(y)[trait_mask]
    studies_labeled = np.asarray(study_ids, dtype=object)[trait_mask]
    if allowed_cols is not None:
        allow = np.asarray(allowed_cols, dtype=np.int64)
        if allow.size == 0:
            return None
        x_labeled = x_full[:, allow]
        # Map local (allowed) indices → full matrix columns after selection.
        local_to_full = allow
    else:
        x_labeled = x_full
        local_to_full = None
    n_cols_input = int(x_labeled.shape[1])
    keep = _univariate_prefilter(
        x_labeled, y_labeled, task=task, max_cols=prefilter_max_cols
    )
    # Cap so a failed frequency gate cannot dump the whole prefilter set as "seeds".
    max_stability_seeds = max(512, min(int(prefilter_max_cols) // 4, 1024))
    cols_local, meta = stability_select_columns(
        x_labeled[:, keep],
        y_labeled,
        study_ids=studies_labeled,
        task=task,
        min_frequency=min_frequency,
        seed=seed,
        n_inner_folds=n_inner_folds,
        n_repeats=n_repeats,
        max_seeds=max_stability_seeds,
    )
    # Remap: prefilter-local → allowed-local → full matrix columns.
    seed_in_allowed = [
        int(keep[int(c)]) for c in np.asarray(cols_local, dtype=np.int64).tolist()
    ]
    if local_to_full is not None:
        seed_cols = [int(local_to_full[i]) for i in seed_in_allowed]
        coef_meta = {
            int(local_to_full[int(keep[int(k)])]): float(v)
            for k, v in meta.get("mean_abs_coef", {}).items()
        }
        freq_meta = {
            int(local_to_full[int(keep[int(k)])]): float(v)
            for k, v in meta.get("frequency", {}).items()
        }
    else:
        seed_cols = seed_in_allowed
        coef_meta = {
            int(keep[int(k)]): float(v) for k, v in meta.get("mean_abs_coef", {}).items()
        }
        freq_meta = {
            int(keep[int(k)]): float(v) for k, v in meta.get("frequency", {}).items()
        }
    coefs = np.asarray([coef_meta.get(c, 0.0) for c in seed_cols], dtype=np.float64)
    # Reject insane SGD scales (age-years × near-singular cols → ~1e10 caps).
    coefs_usable = bool(coefs.size) and float(np.nanmax(np.abs(coefs))) < 1e3
    ranking_fallback = None
    if not coefs_usable:
        # Stability frequencies are uninformative when every col is selected every
        # run with explosive coefs; fall back to univariate-ranked prefilter head.
        ranking_fallback = "univariate_prefilter_top_k"
        seed_in_allowed = [int(c) for c in keep[:max_stability_seeds].tolist()]
        if local_to_full is not None:
            seed_cols = [int(local_to_full[i]) for i in seed_in_allowed]
        else:
            seed_cols = seed_in_allowed
        # Synthetic unit frequencies so gene scoring stays frequency-primary.
        freq_meta = {int(c): 1.0 for c in seed_cols}
        coef_meta = {int(c): 0.0 for c in seed_cols}
        coefs = np.zeros(len(seed_cols), dtype=np.float64)
    cap = (
        float(np.quantile(coefs, strength_cap_quantile)) if coefs_usable else float("nan")
    )

    gene_freq: dict[int, float] = defaultdict(float)
    gene_coefs: dict[int, list[float]] = defaultdict(list)
    gene_seed_cols: dict[int, set[int]] = defaultdict(set)
    for col in seed_cols:
        f = freq_meta.get(col, 0.0)
        raw = coef_meta.get(col, 0.0)
        for gene in col_to_genes.get(col, ()):  # explicit edges only
            gene_freq[gene] += f
            gene_coefs[gene].append(raw)
            gene_seed_cols[gene].add(col)

    genes: list[dict[str, Any]] = []
    for gene, cols_g in gene_seed_cols.items():
        support = sorted(cols_g)
        if support:
            observed = np.isfinite(x_full[:, support]).any(axis=1)
            n_studies = len({str(s) for s in studies_labeled[observed].tolist()})
        else:
            n_studies = 0
        coefs_g = gene_coefs[gene]
        score = float(gene_freq[gene])
        if coefs_usable:
            score = score + 1e-6 * sum(min(c, cap) for c in coefs_g)
        genes.append(
            {
                "gene_index": gene,
                "score": score,
                "freq_score": float(gene_freq[gene]),
                "n_cpgs": len(support),
                "n_studies": int(n_studies),
                "mean_abs_coef": float(np.mean(coefs_g)) if coefs_g else 0.0,
            }
        )
    genes.sort(key=lambda g: (-g["score"], -g["n_cpgs"], g["gene_index"]))
    n_passing = int(meta.get("n_passing_min_frequency", 0))
    return {
        "genes": genes,
        "seed_cols": seed_cols,
        "selection_meta": {
            "n_runs": int(meta.get("n_runs", 0)),
            # ADR 0012: discovery set size (not G2/C2 input width).
            "n_discovery_cpgs": len(seed_cols),
            "n_seed_cpgs_after_stability": len(seed_cols),
            "n_selected_cols": len(seed_cols),
            "n_passing_min_frequency": n_passing,
            "fallback_to_top_freq": bool(n_passing == 0 and len(seed_cols) > 0),
            "max_stability_seeds": int(max_stability_seeds),
            "strength_cap": None if not coefs_usable else cap,
            "strength_cap_warning": bool(not coefs_usable),
            "coefs_numerically_usable": coefs_usable,
            "ranking_fallback": ranking_fallback,
            "prefilter_max_cols": int(prefilter_max_cols),
            "n_cols_input": n_cols_input,
            "n_cols_prefiltered": int(keep.size),
            "n_inner_folds": int(n_inner_folds),
            "n_repeats": int(n_repeats),
            "standardization": meta.get("standardization"),
            "n_zero_variance_dropped": int(meta.get("n_zero_variance_dropped", 0)),
            "n_cols_nonzero_variance": int(meta.get("n_cols_nonzero_variance", 0)),
            "frequency_quantiles": meta.get("frequency_quantiles") or {},
            "coef_abs_quantiles_selected": meta.get("coef_abs_quantiles_selected") or {},
            "n_fits_attempted": int(meta.get("n_fits_attempted", 0)),
            "n_fits_converged": int(meta.get("n_fits_converged", 0)),
            "n_fits_nonconverged": int(meta.get("n_fits_nonconverged", 0)),
            "sparsity_ok": bool(meta.get("sparsity_ok", False)),
            "top_k_per_fit": meta.get("top_k_per_fit"),
            "selection_frequency_by_seed_col": {
                str(k): freq_meta[k] for k in seed_cols if k in freq_meta
            },
        },
    }


def _trait_label_arrays(
    trait_id: str,
    *,
    age: np.ndarray,
    age_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if trait_id == "age":
        return np.asarray(age), np.asarray(age_mask, dtype=bool)
    if trait_id == "tissue":
        return np.asarray(tissue), np.asarray(tissue_mask, dtype=bool)
    if trait_id == "sex":
        return np.asarray(sex), np.asarray(sex_mask, dtype=bool)
    raise ValueError(f"no label arrays wired for trait {trait_id!r}")


def _trait_expansion_stats(
    *,
    trait: str,
    seed_cols_set: set[int],
    locus_records: list[dict[str, Any]],
) -> dict[str, Any]:
    """ADR 0012 discovery-vs-expanded counts for one trait summary."""
    trait_loci = [r for r in locus_records if r["trait"] == trait]
    expanded_cols = {int(r["locus_col"]) for r in trait_loci}
    discovery_in_expanded = {
        int(r["locus_col"]) for r in trait_loci if bool(r["is_seed_cpg"])
    }
    n_edges = len(trait_loci)
    n_unique = len(expanded_cols)
    n_discovery = len(seed_cols_set)
    # Multi-gene: same locus_col linked to >1 gene_id within this trait.
    by_col: dict[int, set[str]] = defaultdict(set)
    for r in trait_loci:
        by_col[int(r["locus_col"])].add(str(r["gene_id"]))
    n_multigene = sum(1 for genes in by_col.values() if len(genes) > 1)
    seed_frac = (
        float(len(discovery_in_expanded) / n_unique) if n_unique > 0 else 0.0
    )
    return {
        "n_discovery_cpgs": int(n_discovery),
        "n_seed_genes": None,  # filled by caller
        "n_expanded_gene_cpg_edges": int(n_edges),
        "n_unique_expanded_gene_cpgs": int(n_unique),
        "n_multigene_cpgs": int(n_multigene),
        "seed_fraction_of_expanded": seed_frac,
        # Aliases for older audits / reports.
        "n_seed_cpgs": int(n_discovery),
        "n_enriched_locus_rows": int(n_edges),
    }


def build_internal_fold_seed_panel(
    *,
    x_train: np.ndarray,
    age: np.ndarray,
    age_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
    study_ids: np.ndarray,
    assignment: CascadeAssignment,
    locus_chrom: np.ndarray | list[str] | None = None,
    gene_length_bp: dict[str, int] | None = None,
    n_genes: int = 256,
    selection_source: Literal["internal_fold"] = "internal_fold",
    fold_id: int | str = 0,
    excluded_study_ids: list[str] | None = None,
    graph_id: str | None = None,
    graph_content_hash: str | None = None,
    matrix_id: str | None = None,
    seed: int = 42,
    min_genes: int = 32,
    traits: Sequence[SeedTraitSpec | Mapping[str, Any] | str] | None = None,
) -> SeedPanelArtifacts:
    """Build fold-fitted per-trait seed-gene panels from outer-train samples only.

    Raises ``ValueError`` if any configured trait with labeled samples yields
    fewer than ``min_genes``, or if a requested trait is unknown / blocked.
    """
    if selection_source != "internal_fold":
        raise ValueError(f"seed_panel only builds internal_fold; got {selection_source!r}")
    if n_genes < 1:
        raise ValueError("n_genes must be >= 1")
    if not graph_content_hash:
        raise ValueError(
            "graph_content_hash is required (hash_graph_tables(graph_root)); "
            "refusing to write an unauditable seed panel"
        )
    trait_specs = resolve_seed_panel_traits(traits)
    col_to_genes, gene_to_cols, edge_role = _gene_edge_maps(assignment)
    gene_ids = list(assignment.gene_ids)
    n_cols = int(np.asarray(x_train).shape[1])

    need_autosome = any(s.id == "sex" and s.autosome_control for s in trait_specs)
    auto_mask: np.ndarray | None = None
    autosomal_cols: np.ndarray | None = None
    if need_autosome:
        auto_mask = _autosomal_col_mask(locus_chrom, n_cols)
        autosomal_cols = np.flatnonzero(auto_mask).astype(np.int64)

    trait_jobs: list[tuple[str, np.ndarray, np.ndarray, np.ndarray | None]] = []
    for spec in trait_specs:
        y, mask = _trait_label_arrays(
            spec.id,
            age=age,
            age_mask=age_mask,
            sex=sex,
            sex_mask=sex_mask,
            tissue=tissue,
            tissue_mask=tissue_mask,
        )
        if not np.asarray(mask, dtype=bool).any():
            raise ValueError(
                f"seed_panel trait {spec.id!r} has no labeled samples on this matrix"
            )
        trait_jobs.append((spec.id, y, mask, None))
        if spec.id == "sex" and spec.autosome_control:
            if autosomal_cols is None:
                raise RuntimeError("sex autosome_control requires autosomal_cols")
            trait_jobs.append(
                ("sex_autosome", y, mask, autosomal_cols)
            )

    gene_records: list[dict[str, Any]] = []
    locus_records: list[dict[str, Any]] = []
    trait_summaries: dict[str, Any] = {}
    overlap_traits = [s.id for s in trait_specs]

    for trait, y, mask, allowed in trait_jobs:
        print(f"[seed_panel] selecting trait={trait} …", flush=True)
        result = _select_trait_genes(
            trait=trait,
            y=y,
            mask=mask,
            x_train=x_train,
            study_ids=study_ids,
            col_to_genes=col_to_genes,
            seed=seed,
            strength_cap_quantile=DEFAULT_STRENGTH_CAP_QUANTILE,
            min_frequency=DEFAULT_MIN_FREQUENCY,
            allowed_cols=allowed,
        )
        if result is None:
            raise ValueError(
                f"seed_panel trait {trait!r} produced no selection "
                "(empty after study filters?)"
            )
        ranked = result["genes"][:n_genes]
        seed_cols_set = set(result["seed_cols"])
        if len(ranked) < min_genes:
            raise ValueError(
                f"seed panel for {trait!r} has {len(ranked)} genes < min_genes={min_genes}"
            )
        for rank, gene in enumerate(ranked):
            gene_idx = int(gene["gene_index"])
            gene_id = gene_ids[gene_idx] if 0 <= gene_idx < len(gene_ids) else str(gene_idx)
            all_cols = sorted(gene_to_cols.get(gene_idx, set()))
            if trait == "sex_autosome":
                if auto_mask is None:
                    raise RuntimeError("sex_autosome requires auto_mask")
                all_cols = [c for c in all_cols if bool(auto_mask[int(c)])]
                gene_seed = seed_cols_set & set(gene_to_cols.get(gene_idx, set()))
                if any(not bool(auto_mask[int(c)]) for c in gene_seed):
                    raise AssertionError(
                        f"sex_autosome seed CpG on sex chrom for {gene_id}"
                    )
            roles = {edge_role.get((c, gene_idx), "unknown") for c in all_cols}
            autosomal = _gene_chrom_is_autosomal(set(all_cols), locus_chrom)
            if trait == "sex_autosome" and autosomal is not True:
                raise AssertionError(
                    f"sex_autosome gene {gene_id} is not autosomal-only (autosome_only={autosomal})"
                )
            gene_records.append(
                {
                    "gene_id": gene_id,
                    "trait": trait,
                    "rank": rank,
                    "score": gene["score"],
                    # Discovery support (stability seed CpGs for this gene).
                    "n_cpgs": gene["n_cpgs"],
                    # Full sibling-enriched width used by G2/C2.
                    "n_expanded_cpgs": len(all_cols),
                    "n_studies": gene["n_studies"],
                    "mean_abs_coef": gene["mean_abs_coef"],
                    "promoter_body_coverage": _promoter_body_coverage(roles),
                    "inclusion_reason": "stability_selection",
                    "autosome_only": autosomal,
                }
            )
            locus_records.extend(
                {
                    "gene_id": gene_id,
                    "trait": trait,
                    "locus_col": int(col),
                    "locus_id": str(col),
                    "gene_role": edge_role.get((col, gene_idx), "unknown"),
                    "is_seed_cpg": col in seed_cols_set,
                }
                for col in all_cols
            )
        sex_chrom_seed = 0
        if locus_chrom is not None:
            sex_chrom_seed = sum(
                1 for c in seed_cols_set if _is_sex_chrom(str(np.asarray(locus_chrom)[int(c)]))
            )
        n_autosome_genes = sum(
            1
            for r in gene_records
            if r["trait"] == trait and r["autosome_only"] is True
        )
        expansion = _trait_expansion_stats(
            trait=trait,
            seed_cols_set=seed_cols_set,
            locus_records=locus_records,
        )
        expansion["n_seed_genes"] = len(ranked)
        trait_summaries[trait] = {
            "n_genes_requested": n_genes,
            "n_genes_actual": len(ranked),
            "n_sex_chrom_seed_cpgs": int(sex_chrom_seed),
            "n_autosome_only_genes": int(n_autosome_genes),
            **result["selection_meta"],
            **expansion,
        }
        if trait == "sex_autosome":
            if sex_chrom_seed != 0:
                raise AssertionError(
                    f"sex_autosome control has {sex_chrom_seed} sex-chrom seed CpGs"
                )
            trait_summaries[trait]["autosome_only_control"] = True

    if not trait_summaries:
        raise ValueError("no labeled samples for any trait; cannot build seed panel")

    genes_df = pd.DataFrame(
        gene_records,
        columns=[
            "gene_id",
            "trait",
            "rank",
            "score",
            "n_cpgs",
            "n_expanded_cpgs",
            "n_studies",
            "mean_abs_coef",
            "promoter_body_coverage",
            "inclusion_reason",
            "autosome_only",
        ],
    )
    loci_df = pd.DataFrame(
        locus_records,
        columns=["gene_id", "trait", "locus_col", "locus_id", "gene_role", "is_seed_cpg"],
    )

    panel_hash = _panel_hash(genes_df, loci_df)
    overlap = summarize_seed_panel_overlap(
        genes_df, loci_df, traits=overlap_traits
    )
    panel_json: dict[str, Any] = {
        "selection_source": selection_source,
        "method": SELECTION_METHOD,
        "fold_id": fold_id,
        "seed": seed,
        "n_genes_requested": n_genes,
        "min_genes": min_genes,
        "gene_allocation": "explicit_only",
        "configured_traits": [
            {
                "id": s.id,
                "role": s.role,
                "autosome_control": s.autosome_control,
            }
            for s in trait_specs
        ],
        "thresholds": {
            "min_frequency": DEFAULT_MIN_FREQUENCY,
            "strength_cap_quantile": DEFAULT_STRENGTH_CAP_QUANTILE,
        },
        "excluded_study_ids": list(excluded_study_ids or []),
        "graph_id": graph_id,
        "graph_content_hash": graph_content_hash,
        "matrix_id": matrix_id,
        "n_train_samples": int(np.asarray(x_train).shape[0]),
        "gene_length_bp_provided": gene_length_bp is not None,
        "traits": trait_summaries,
        "overlap": overlap,
        "panel_hash": panel_hash,
    }
    return SeedPanelArtifacts(
        panel_json=panel_json,
        genes=genes_df,
        loci=loci_df,
        panel_hash=panel_hash,
    )


def summarize_seed_panel_overlap(
    genes: pd.DataFrame,
    loci: pd.DataFrame,
    *,
    traits: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Gene / CpG overlap and coverage stats across configured traits.

    ``sex_autosome`` is excluded unless explicitly listed in ``traits``.
    """
    if traits is None:
        present = [
            str(t)
            for t in genes["trait"].astype(str).unique().tolist()
            if str(t) != "sex_autosome"
        ]
        primary = tuple(present) if present else ("age", "tissue", "sex")
    else:
        primary = tuple(str(t) for t in traits if str(t) != "sex_autosome")
    if not primary:
        primary = ("age", "tissue", "sex")
    gene_sets = {
        t: set(genes.loc[genes["trait"] == t, "gene_id"].astype(str)) for t in primary
    }
    cpg_sets = {
        t: set(
            loci.loc[loci["trait"] == t, "locus_col"].astype(int).tolist()
        )
        for t in primary
    }
    seed_cpg_sets = {
        t: set(
            loci.loc[
                (loci["trait"] == t) & (loci["is_seed_cpg"].astype(bool)),
                "locus_col",
            ]
            .astype(int)
            .tolist()
        )
        for t in primary
    }
    role_cov: dict[str, dict[str, int]] = {}
    for t in primary:
        sub = loci.loc[loci["trait"] == t]
        role_cov[t] = {
            str(k): int(v) for k, v in sub["gene_role"].value_counts().to_dict().items()
        }
    genes_one_cpg = {
        t: int((genes.loc[genes["trait"] == t, "n_cpgs"] == 1).sum()) for t in primary
    }
    seed_fraction = {
        t: (
            float(len(seed_cpg_sets[t]) / len(cpg_sets[t])) if cpg_sets[t] else 0.0
        )
        for t in primary
    }
    # Multi-gene CpGs: same locus_col linked to >1 gene_id within a trait.
    multi_gene_cpg: dict[str, int] = {}
    for t in primary:
        sub = loci.loc[loci["trait"] == t, ["locus_col", "gene_id"]]
        if sub.empty:
            multi_gene_cpg[t] = 0
            continue
        counts = sub.groupby("locus_col")["gene_id"].nunique()
        multi_gene_cpg[t] = int((counts > 1).sum())
    pairwise = [
        (primary[i], primary[j])
        for i in range(len(primary))
        for j in range(i + 1, len(primary))
    ]
    nonempty_gene = [gene_sets[t] for t in primary if gene_sets[t]]
    nonempty_cpg = [cpg_sets[t] for t in primary if cpg_sets[t]]
    return {
        "traits": list(primary),
        "gene_set_sizes": {t: len(gene_sets[t]) for t in primary},
        "gene_union_size": len(set.union(*nonempty_gene)) if nonempty_gene else 0,
        "gene_pairwise_overlap": {
            f"{a}_∩_{b}": len(gene_sets[a] & gene_sets[b]) for a, b in pairwise
        },
        "cpg_set_sizes": {t: len(cpg_sets[t]) for t in primary},
        "cpg_union_size": len(set.union(*nonempty_cpg)) if nonempty_cpg else 0,
        "cpg_pairwise_overlap": {
            f"{a}_∩_{b}": len(cpg_sets[a] & cpg_sets[b]) for a, b in pairwise
        },
        "stability_seed_cpg_set_sizes": {t: len(seed_cpg_sets[t]) for t in primary},
        "n_unique_expanded_gene_cpgs": {t: len(cpg_sets[t]) for t in primary},
        "seed_fraction_of_expanded": seed_fraction,
        "gene_role_coverage": role_cov,
        "genes_with_only_one_seed_cpg": genes_one_cpg,
        "multi_gene_cpg_count": multi_gene_cpg,
    }


def _panel_hash(genes: pd.DataFrame, loci: pd.DataFrame) -> str:
    """Content hash over ordered gene / locus records (order-independent)."""
    digest = hashlib.sha256()
    gene_rows: list[dict[str, Any]] = genes.to_dict(orient="records")
    locus_rows: list[dict[str, Any]] = loci.to_dict(orient="records")
    gene_rows.sort(key=lambda r: (str(r["trait"]), str(r["gene_id"])))
    locus_rows.sort(key=lambda r: (str(r["trait"]), str(r["gene_id"]), int(r["locus_col"])))
    digest.update(json.dumps(gene_rows, sort_keys=True, default=str).encode("utf-8"))
    digest.update(b"\x00")
    digest.update(json.dumps(locus_rows, sort_keys=True, default=str).encode("utf-8"))
    return digest.hexdigest()[:32]


def matched_random_gene_panel(
    seed_gene_ids: list[str] | list[int],
    *,
    candidate_gene_ids: Sequence[str] | Sequence[int],
    gene_cpg_counts: dict[Any, int],
    gene_length_bp: dict[Any, int] | None = None,
    gene_role_coverage: dict[Any, Sequence[float]] | None = None,
    rng: np.random.Generator,
) -> list[Any]:
    """G3 arm: match each seed gene to a distinct non-seed candidate.

    Match on CpG count (primary), gene length and promoter/body coverage when
    available. Greedy nearest-neighbour without replacement; ties broken with
    ``rng``. Never returns a seed gene.
    """
    matched, _quality = matched_random_gene_panel_with_quality(
        seed_gene_ids,
        candidate_gene_ids=candidate_gene_ids,
        gene_cpg_counts=gene_cpg_counts,
        gene_length_bp=gene_length_bp,
        gene_role_coverage=gene_role_coverage,
        rng=rng,
    )
    return matched


def matched_random_gene_panel_with_quality(
    seed_gene_ids: list[str] | list[int],
    *,
    candidate_gene_ids: Sequence[str] | Sequence[int],
    gene_cpg_counts: dict[Any, int],
    gene_length_bp: dict[Any, int] | None = None,
    gene_role_coverage: dict[Any, Sequence[float]] | None = None,
    rng: np.random.Generator,
) -> tuple[list[Any], dict[str, Any]]:
    """Like :func:`matched_random_gene_panel` plus CpG-count matching diagnostics."""
    seed_set = set(seed_gene_ids)
    available = [g for g in candidate_gene_ids if g not in seed_set]
    if len(available) < len(seed_gene_ids):
        raise ValueError(
            f"need {len(seed_gene_ids)} matched genes but only {len(available)} candidates"
        )

    def _distance(seed: Any, cand: Any) -> tuple[float, float, float]:
        cpg_diff = abs(gene_cpg_counts.get(seed, 0) - gene_cpg_counts.get(cand, 0))
        len_diff = 0.0
        if gene_length_bp is not None:
            len_diff = abs(gene_length_bp.get(seed, 0) - gene_length_bp.get(cand, 0))
        role_diff = 0.0
        if gene_role_coverage is not None:
            a = np.asarray(gene_role_coverage.get(seed, ()), dtype=np.float64)
            b = np.asarray(gene_role_coverage.get(cand, ()), dtype=np.float64)
            if a.shape == b.shape and a.size:
                role_diff = float(np.abs(a - b).sum())
        return float(cpg_diff), float(len_diff), role_diff

    used: set[Any] = set()
    matched: dict[Any, Any] = {}
    cpg_abs_err: list[float] = []
    for seed in seed_gene_ids:
        pool = [c for c in available if c not in used]
        best = min(pool, key=lambda c: (*_distance(seed, c), rng.random()))
        used.add(best)
        matched[seed] = best
        cpg_abs_err.append(
            abs(float(gene_cpg_counts.get(seed, 0) - gene_cpg_counts.get(best, 0)))
        )
    err = np.asarray(cpg_abs_err, dtype=np.float64)
    quality = {
        "n_seed_genes": len(seed_gene_ids),
        "n_matched": len(matched),
        "seed_genes_disjoint_from_matched": bool(
            set(seed_gene_ids).isdisjoint(matched.values())
        ),
        "cpg_count_abs_err_mean": float(err.mean()) if err.size else 0.0,
        "cpg_count_abs_err_median": float(np.median(err)) if err.size else 0.0,
        "cpg_count_abs_err_p90": float(np.quantile(err, 0.9)) if err.size else 0.0,
        "cpg_count_abs_err_max": float(err.max()) if err.size else 0.0,
        "fraction_exact_cpg_match": (
            float(np.mean(err == 0.0)) if err.size else 0.0
        ),
        "gene_length_bp_used": gene_length_bp is not None,
        "gene_role_coverage_used": gene_role_coverage is not None,
    }
    return [matched[s] for s in seed_gene_ids], quality


def gene_mask_tensor(
    selected_gene_indices: Sequence[int],
    n_genes: int,
    *,
    n_outputs: int = 1,
) -> np.ndarray:
    """Return float32 ``[n_outputs, n_genes]`` mask (1.0 at selected genes)."""
    if n_genes < 0:
        raise ValueError("n_genes must be >= 0")
    if n_outputs < 1:
        raise ValueError("n_outputs must be >= 1")
    mask = np.zeros((n_outputs, n_genes), dtype=np.float32)
    idx = np.asarray(list(selected_gene_indices), dtype=np.int64)
    if idx.size:
        if int(idx.min()) < 0 or int(idx.max()) >= n_genes:
            raise ValueError("selected_gene_indices out of range [0, n_genes)")
        mask[:, idx] = 1.0
    return mask


def write_seed_panel(out_dir: Path, artifacts: SeedPanelArtifacts) -> dict[str, str]:
    """Write ``seed_panel.json`` + gene / locus parquet; return paths + panel_hash."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "seed_panel.json"
    gene_path = out_dir / "seed_panel_gene.parquet"
    locus_path = out_dir / "seed_panel_locus.parquet"
    panel_path.write_text(
        json.dumps(artifacts.panel_json, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    artifacts.genes.to_parquet(gene_path, index=False)
    artifacts.loci.to_parquet(locus_path, index=False)
    return {
        "seed_panel_json": str(panel_path),
        "seed_panel_gene_parquet": str(gene_path),
        "seed_panel_locus_parquet": str(locus_path),
        "panel_hash": artifacts.panel_hash,
    }


def load_seed_panel(out_dir: Path) -> SeedPanelArtifacts:
    """Reload hashed panel artifacts written by :func:`write_seed_panel`."""
    out_dir = Path(out_dir)
    panel_json = json.loads((out_dir / "seed_panel.json").read_text(encoding="utf-8"))
    genes = pd.read_parquet(out_dir / "seed_panel_gene.parquet")
    loci = pd.read_parquet(out_dir / "seed_panel_locus.parquet")
    panel_hash = str(panel_json.get("panel_hash") or _panel_hash(genes, loci))
    return SeedPanelArtifacts(
        panel_json=panel_json, genes=genes, loci=loci, panel_hash=panel_hash
    )
