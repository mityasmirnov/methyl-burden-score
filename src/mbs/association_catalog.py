"""Atlas association + seed-panel knowledge track (CPU-only, no training).

Parallel, non-blocking track that ingests published CpG↔trait associations and
derives per-trait *seed gene panels* used for external, leakage-free feature
selection. Nothing here reads sample×CpG observations (ADR 0005). Atlas gene
symbols are kept as source metadata only and never drive gene allocation
(ADR 0006, ADR 0010); seed panels are built from explicit, non-nearest-gene
edges only (ADR 0004, ADR 0010 ``explicit_only``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from mbs.atlas_study_enrichment import ATLAS_STUDY_ID_RE

ATLAS_ASSOCIATIONS_FILE = "EWAS_Atlas_associations.tsv"

# Candidate source column names → normalized names used across this module.
_PROBE_COLS = ("probe_id", "probe_ID", "probeID", "cpg", "CpG")
_STUDY_COLS = ("association_study_id", "study_ID", "study_id", "atlas_study_id")
_GENE_SYMBOL_COLS = ("atlas_gene_symbol", "gene_symbol", "gene", "Gene")
_DIRECTION_COLS = ("effect_direction", "correlation", "direction")


def _first_present(frame: pd.DataFrame, candidates: Sequence[str]) -> str | None:
    """Return the first candidate column present in ``frame`` (or ``None``)."""
    return next((c for c in candidates if c in frame.columns), None)


def load_atlas_associations_head(
    atlas_root: Path,
    *,
    nrows: int | None = 1000,
) -> pd.DataFrame:
    """Read ``EWAS_Atlas_associations.tsv`` (may be huge — ``nrows`` bounds it).

    Pass ``nrows=None`` to read the whole file. Raises ``FileNotFoundError`` if
    the associations TSV is absent so callers fail loudly rather than silently
    returning an empty panel.
    """
    path = atlas_root / ATLAS_ASSOCIATIONS_FILE
    if not path.is_file():
        raise FileNotFoundError(f"Atlas associations file not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, encoding="latin-1", nrows=nrows)


def _explicit_non_nearest_edges(locus_gene_edges: pd.DataFrame) -> pd.DataFrame:
    """Keep only explicit, non-nearest-gene edges (seed-panel safety gate)."""
    missing = {"is_explicit", "is_nearest_gene"} - set(locus_gene_edges.columns)
    if missing:
        raise ValueError(
            "locus_gene_edges must carry is_explicit and is_nearest_gene "
            f"(missing: {sorted(missing)})"
        )
    explicit = locus_gene_edges["is_explicit"].astype(bool)
    nearest = locus_gene_edges["is_nearest_gene"].astype(bool)
    return locus_gene_edges.loc[explicit & ~nearest].copy()


def remap_associations_to_graph(
    associations: pd.DataFrame,
    *,
    locus_gene_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Remap associations onto the annotation graph via explicit edges only.

    The Atlas gene symbol (if any) is retained as ``atlas_gene_symbol`` metadata
    and is **never** used as the join key or allocation. Associations are joined
    to ``locus_gene_edges`` on ``locus_id`` (preferred) or ``probe_id`` and only
    ``is_explicit`` / non-``is_nearest_gene`` edges survive, so nearest-gene
    edges cannot seed a panel (ADR 0004/0010).
    """
    edges = _explicit_non_nearest_edges(locus_gene_edges)

    assoc = associations.copy()
    probe_col = _first_present(assoc, _PROBE_COLS)
    if probe_col and probe_col != "probe_id":
        assoc = assoc.rename(columns={probe_col: "probe_id"})
    gene_symbol_col = _first_present(assoc, _GENE_SYMBOL_COLS)
    if gene_symbol_col and gene_symbol_col != "atlas_gene_symbol":
        assoc = assoc.rename(columns={gene_symbol_col: "atlas_gene_symbol"})
    if "atlas_gene_symbol" not in assoc.columns:
        assoc["atlas_gene_symbol"] = pd.NA

    # Join key: never ES/gene-symbol; use locus_id if both sides have it, else probe_id.
    if "locus_id" in assoc.columns and "locus_id" in edges.columns:
        key = "locus_id"
    elif "probe_id" in assoc.columns and "probe_id" in edges.columns:
        key = "probe_id"
    else:
        raise ValueError("no shared locus_id/probe_id join key between associations and edges")

    edge_cols = [
        c
        for c in (
            "locus_id",
            "probe_id",
            "gene_id",
            "gene_role",
            "mapping_evidence",
            "mapping_source",
            "cpg_context",
        )
        if c in edges.columns
    ]
    return assoc.merge(edges[edge_cols], on=key, how="inner", suffixes=("", "_edge"))


def _split_atlas_ids(value: object) -> list[str]:
    """Split a ``;``-joined atlas_study_id cell into ES* tokens."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [tok.strip() for tok in str(value).split(";") if tok.strip()]


def external_clean_study_ids(
    atlas_study_ids: Sequence[str],
    *,
    benchmark_gse_ids: set[str],
    atlas_gse_map: pd.DataFrame,
    pmid_to_gse: Mapping[str, list[str]] | None = None,
) -> list[str]:
    """Atlas study IDs with no overlap to the benchmark GSE set.

    An Atlas study (ES*) is contaminated if the curated map or a shared PMID
    links it to any benchmark GSE. We never compare ES* to GSE* by string
    equality — only via the curated ``atlas_gse_map`` / PMID bridge (see
    ``mbs.atlas_study_enrichment`` and ``docs/EWAS_METADATA.md``).
    """
    benchmark = {g.strip() for g in benchmark_gse_ids if g and g.strip()}
    contaminated_es: set[str] = set()
    contaminated_pmids: set[str] = set()

    has_map = not atlas_gse_map.empty and "gse_id" in atlas_gse_map.columns
    if has_map:
        for _, row in atlas_gse_map.iterrows():
            gse = str(row.get("gse_id", "")).strip()
            es_ids = _split_atlas_ids(row.get("atlas_study_id"))
            pmid = str(row.get("pmid", "")).strip()
            if gse in benchmark:
                contaminated_es.update(es_ids)
                if pmid:
                    contaminated_pmids.update(p.strip() for p in pmid.split(";") if p.strip())

    if pmid_to_gse:
        for pmid, gses in pmid_to_gse.items():
            if any(str(g).strip() in benchmark for g in gses):
                contaminated_pmids.add(str(pmid).strip())

    # Second pass: any Atlas study sharing a contaminated PMID is also excluded.
    if has_map and contaminated_pmids and "pmid" in atlas_gse_map.columns:
        for _, row in atlas_gse_map.iterrows():
            row_pmids = {p.strip() for p in str(row.get("pmid", "")).split(";") if p.strip()}
            if row_pmids & contaminated_pmids:
                contaminated_es.update(_split_atlas_ids(row.get("atlas_study_id")))

    clean: list[str] = []
    seen: set[str] = set()
    for raw in atlas_study_ids:
        es = str(raw).strip()
        if not es or es in seen:
            continue
        seen.add(es)
        if not ATLAS_STUDY_ID_RE.match(es):
            continue
        if es.upper() in {c.upper() for c in contaminated_es}:
            continue
        clean.append(es)
    return clean


def build_external_clean_gene_list(
    remapped_associations: pd.DataFrame,
    *,
    clean_study_ids: Sequence[str],
    p_value_threshold: float = 1e-5,
) -> pd.DataFrame:
    """Aggregate explicit, leakage-free associations into a ranked gene panel.

    Operates on the output of :func:`remap_associations_to_graph` restricted to
    externally clean studies (see :func:`external_clean_study_ids`). Returns one
    row per gene with ``seed_panel_gene``-aligned columns. Works on a small
    fixture frame; the full Atlas file is bounded upstream via ``nrows``.
    """
    if "gene_id" not in remapped_associations.columns:
        raise ValueError("remapped_associations must carry gene_id (explicit edges)")
    frame = remapped_associations.copy()

    study_col = _first_present(frame, _STUDY_COLS)
    if study_col is None:
        raise ValueError(f"no study-id column found among {list(_STUDY_COLS)}")
    clean = {str(s).strip() for s in clean_study_ids}
    frame = frame.loc[frame[study_col].astype(str).str.strip().isin(clean)].copy()

    probe_col = _first_present(frame, _PROBE_COLS) or "probe_id"
    if "p_value" in frame.columns:
        frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
        frame = frame.loc[frame["p_value"] <= p_value_threshold]
    fdr_col = (
        "fdr" if "fdr" in frame.columns else ("p_value" if "p_value" in frame.columns else None)
    )
    direction_col = _first_present(frame, _DIRECTION_COLS)

    rows: list[dict[str, object]] = []
    for gene_id, grp in frame.groupby("gene_id", sort=True):
        n_cpgs = int(grp[probe_col].nunique()) if probe_col in grp.columns else int(len(grp))
        n_studies = int(grp[study_col].nunique())
        min_fdr = float(grp[fdr_col].min()) if fdr_col and grp[fdr_col].notna().any() else None
        if direction_col and grp[direction_col].notna().any():
            counts = grp[direction_col].astype(str).str.strip().value_counts()
            direction_consistency = float(counts.iloc[0] / counts.sum())
        else:
            direction_consistency = None
        rows.append(
            {
                "gene_id": str(gene_id),
                "score": float(n_studies * n_cpgs),
                "n_associated_cpgs": n_cpgs,
                "n_independent_studies": n_studies,
                "direction_consistency": direction_consistency,
                "min_fdr": min_fdr,
                "inclusion_reason": "external_clean",
            }
        )

    panel = pd.DataFrame(
        rows,
        columns=[
            "gene_id",
            "score",
            "n_associated_cpgs",
            "n_independent_studies",
            "direction_consistency",
            "min_fdr",
            "inclusion_reason",
        ],
    )
    panel = panel.sort_values(
        ["n_independent_studies", "n_associated_cpgs", "score"],
        ascending=[False, False, False],
        kind="stable",
    ).reset_index(drop=True)
    panel.insert(1, "rank", range(1, len(panel) + 1))
    return panel


def build_hybrid_fold_gene_list(
    *,
    external_clean: pd.DataFrame,
    internal_fold: pd.DataFrame,
    n_genes: int | None = None,
) -> pd.DataFrame:
    """Combine external_clean + internal_fold gene lists inside the train fold.

    Union by ``gene_id``; prefer the higher score when both sources nominate a
    gene. Marks ``inclusion_reason`` as ``hybrid_fold`` (or
    ``hybrid_fold+external`` / ``hybrid_fold+internal`` when only one source
    contributes). Ranking is by score descending. Callers must ensure
    ``internal_fold`` used only outer-train samples (ADR 0011).
    """
    need = {"gene_id", "score"}
    for name, frame in (("external_clean", external_clean), ("internal_fold", internal_fold)):
        missing = need - set(frame.columns)
        if missing:
            raise ValueError(f"{name} missing columns: {sorted(missing)}")

    ext = external_clean.copy()
    inn = internal_fold.copy()
    ext["gene_id"] = ext["gene_id"].astype(str)
    inn["gene_id"] = inn["gene_id"].astype(str)
    ext["_src"] = "external"
    inn["_src"] = "internal"
    keep_cols = ["gene_id", "score", "_src"]
    for optional in (
        "n_associated_cpgs",
        "n_independent_studies",
        "direction_consistency",
        "min_fdr",
    ):
        if optional in ext.columns or optional in inn.columns:
            keep_cols.append(optional)
            if optional not in ext.columns:
                ext[optional] = None
            if optional not in inn.columns:
                inn[optional] = None

    stacked = pd.concat([ext[keep_cols], inn[keep_cols]], ignore_index=True)
    stacked["score"] = pd.to_numeric(stacked["score"], errors="coerce").fillna(0.0)
    rows: list[dict[str, object]] = []
    for gene_id, grp in stacked.groupby("gene_id", sort=False):
        sources = set(grp["_src"].astype(str))
        best = grp.loc[grp["score"].idxmax()]
        if sources == {"external", "internal"}:
            reason = "hybrid_fold"
        elif "external" in sources:
            reason = "hybrid_fold+external"
        else:
            reason = "hybrid_fold+internal"
        row: dict[str, object] = {
            "gene_id": str(gene_id),
            "score": float(best["score"]),
            "inclusion_reason": reason,
        }
        for optional in (
            "n_associated_cpgs",
            "n_independent_studies",
            "direction_consistency",
            "min_fdr",
        ):
            if optional in best.index and pd.notna(best[optional]):
                row[optional] = best[optional]
        rows.append(row)

    panel = pd.DataFrame(rows)
    panel = panel.sort_values("score", ascending=False, kind="stable").reset_index(drop=True)
    if n_genes is not None:
        panel = panel.iloc[: int(n_genes)].copy()
    panel.insert(1, "rank", range(1, len(panel) + 1))
    return panel
