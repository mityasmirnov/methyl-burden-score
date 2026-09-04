"""``internal_fold`` seed-gene panel constructor (7G′ Stage B, ADR 0010/0011).

Gene-first, fold-fitted seed selection: reuse ``stability_select_columns`` on
outer-train samples only, map selected CpG columns to genes via *explicit*
region→gene edges, rank genes by capped burden strength, then enrich each
selected gene with all its gene-linked CpGs (siblings). DeepRVAT-faithful:
selection is a burden over the gene's CpG set, not a single min-P locus.

Callers should pass an ``explicit_only`` :class:`CascadeAssignment` (ADR 0010)
so neural and classical arms share evidence-backed gene columns.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from mbs.training.cascade_assign import CascadeAssignment
from mbs.training.fold_safe_panel import stability_select_columns

TRAITS: tuple[str, ...] = ("age", "tissue", "sex")
_PROMOTER_ROLES = frozenset({"promoter_core", "promoter_proximal", "five_prime"})
_BODY_ROLES = frozenset({"gene_body", "three_prime"})
_SEX_CHROMS = frozenset({"chrx", "chry", "x", "y"})
DEFAULT_STRENGTH_CAP_QUANTILE = 0.99
DEFAULT_MIN_FREQUENCY = 0.34
# ponytail: full 51k-col × 5×3×9 enet grid is multi-hour; univariate prefilter
# + lighter CV keeps DeepRVAT-faithful fold-fitting while finishing a screen.
DEFAULT_PREFILTER_MAX_COLS = 4096
DEFAULT_N_INNER_FOLDS = 2
DEFAULT_N_REPEATS = 2
SELECTION_METHOD = "study_grouped_enet_stability_gene_first"


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
) -> dict[str, Any] | None:
    """Run stability selection for one trait; score + rank supported genes."""
    trait_mask = np.asarray(mask, dtype=bool)
    if trait == "tissue":
        trait_mask = trait_mask & _tissue_multi_study_mask(y, trait_mask, study_ids)
    if not trait_mask.any():
        return None
    task: Literal["age", "sex", "tissue"] = (
        "age" if trait == "age" else "sex" if trait == "sex" else "tissue"
    )
    x_labeled = x_train[trait_mask]
    y_labeled = np.asarray(y)[trait_mask]
    studies_labeled = np.asarray(study_ids, dtype=object)[trait_mask]
    keep = _univariate_prefilter(
        x_labeled, y_labeled, task=task, max_cols=prefilter_max_cols
    )
    cols_local, meta = stability_select_columns(
        x_labeled[:, keep],
        y_labeled,
        study_ids=studies_labeled,
        task=task,
        min_frequency=min_frequency,
        seed=seed,
        n_inner_folds=n_inner_folds,
        n_repeats=n_repeats,
    )
    # Remap local prefilter indices back to matrix columns.
    seed_cols = [int(keep[int(c)]) for c in np.asarray(cols_local, dtype=np.int64).tolist()]
    coef_meta: dict[int, float] = {
        int(keep[int(k)]): float(v) for k, v in meta.get("mean_abs_coef", {}).items()
    }
    coefs = np.asarray([coef_meta.get(c, 0.0) for c in seed_cols], dtype=np.float64)
    # ponytail: data-adaptive burden cap (upper quantile) so one runaway CpG
    # cannot dominate a gene's rank; a fixed global cap is the upgrade path.
    cap = float(np.quantile(coefs, strength_cap_quantile)) if coefs.size else 0.0

    gene_score: dict[int, float] = defaultdict(float)
    gene_coefs: dict[int, list[float]] = defaultdict(list)
    gene_seed_cols: dict[int, set[int]] = defaultdict(set)
    for col in seed_cols:
        raw = coef_meta.get(col, 0.0)
        capped = min(raw, cap) if cap > 0 else raw
        for gene in col_to_genes.get(col, ()):  # explicit edges only
            gene_score[gene] += capped
            gene_coefs[gene].append(raw)
            gene_seed_cols[gene].add(col)

    genes: list[dict[str, Any]] = []
    for gene, cols_g in gene_seed_cols.items():
        support = sorted(cols_g)
        if support:
            observed = np.isfinite(x_labeled[:, support]).any(axis=1)
            n_studies = len({str(s) for s in studies_labeled[observed].tolist()})
        else:
            n_studies = 0
        coefs_g = gene_coefs[gene]
        genes.append(
            {
                "gene_index": gene,
                "score": float(gene_score[gene]),
                "n_cpgs": len(support),
                "n_studies": int(n_studies),
                "mean_abs_coef": float(np.mean(coefs_g)) if coefs_g else 0.0,
            }
        )
    # Rank by burden strength, then supporting-CpG count, then index (stable).
    genes.sort(key=lambda g: (-g["score"], -g["n_cpgs"], g["gene_index"]))
    return {
        "genes": genes,
        "seed_cols": seed_cols,
        "selection_meta": {
            "n_runs": int(meta.get("n_runs", 0)),
            "n_selected_cols": len(seed_cols),
            "strength_cap": cap,
            "prefilter_max_cols": int(prefilter_max_cols),
            "n_cols_prefiltered": int(keep.size),
            "n_inner_folds": int(n_inner_folds),
            "n_repeats": int(n_repeats),
        },
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
) -> SeedPanelArtifacts:
    """Build fold-fitted per-trait seed-gene panels from outer-train samples only.

    Raises ``ValueError`` if any trait with labeled samples yields < ``min_genes``.
    """
    if selection_source != "internal_fold":
        raise ValueError(f"seed_panel only builds internal_fold; got {selection_source!r}")
    if n_genes < 1:
        raise ValueError("n_genes must be >= 1")
    col_to_genes, gene_to_cols, edge_role = _gene_edge_maps(assignment)
    gene_ids = list(assignment.gene_ids)

    trait_inputs = {
        "age": (np.asarray(age), np.asarray(age_mask, dtype=bool)),
        "tissue": (np.asarray(tissue), np.asarray(tissue_mask, dtype=bool)),
        "sex": (np.asarray(sex), np.asarray(sex_mask, dtype=bool)),
    }

    gene_records: list[dict[str, Any]] = []
    locus_records: list[dict[str, Any]] = []
    trait_summaries: dict[str, Any] = {}

    for trait in TRAITS:
        y, mask = trait_inputs[trait]
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
        )
        if result is None:
            continue
        ranked = result["genes"][:n_genes]
        seed_cols_set = set(result["seed_cols"])
        if len(ranked) < min_genes:
            raise ValueError(
                f"seed panel for {trait!r} has {len(ranked)} genes < min_genes={min_genes}"
            )
        for rank, gene in enumerate(ranked):
            gene_idx = int(gene["gene_index"])
            gene_id = gene_ids[gene_idx] if 0 <= gene_idx < len(gene_ids) else str(gene_idx)
            all_cols = sorted(gene_to_cols.get(gene_idx, set()))  # enrichment: siblings too
            roles = {edge_role.get((c, gene_idx), "unknown") for c in all_cols}
            autosomal = (
                _gene_chrom_is_autosomal(set(all_cols), locus_chrom) if trait == "sex" else None
            )
            gene_records.append(
                {
                    "gene_id": gene_id,
                    "trait": trait,
                    "rank": rank,
                    "score": gene["score"],
                    "n_cpgs": gene["n_cpgs"],
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
        n_autosome_genes = (
            sum(1 for r in gene_records if r["trait"] == "sex" and r["autosome_only"] is True)
            if trait == "sex"
            else None
        )
        trait_summaries[trait] = {
            "n_genes_requested": n_genes,
            "n_genes_actual": len(ranked),
            "n_seed_cpgs": len(seed_cols_set),
            **result["selection_meta"],
            **({"n_autosome_only_genes": n_autosome_genes} if n_autosome_genes is not None else {}),
        }

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
    panel_json: dict[str, Any] = {
        "selection_source": selection_source,
        "method": SELECTION_METHOD,
        "fold_id": fold_id,
        "seed": seed,
        "n_genes_requested": n_genes,
        "min_genes": min_genes,
        "gene_allocation": "explicit_only",
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
        "panel_hash": panel_hash,
    }
    return SeedPanelArtifacts(
        panel_json=panel_json,
        genes=genes_df,
        loci=loci_df,
        panel_hash=panel_hash,
    )


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
    for seed in seed_gene_ids:
        pool = [c for c in available if c not in used]
        best = min(pool, key=lambda c: (*_distance(seed, c), rng.random()))
        used.add(best)
        matched[seed] = best
    return [matched[s] for s in seed_gene_ids]


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
