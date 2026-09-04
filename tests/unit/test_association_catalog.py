from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd
import pytest

from mbs.association_catalog import (
    ATLAS_ASSOCIATIONS_FILE,
    build_external_clean_gene_list,
    build_hybrid_fold_gene_list,
    external_clean_study_ids,
    load_atlas_associations_head,
    remap_associations_to_graph,
)


def _write_atlas_associations(atlas_root: Path, n_rows: int = 5) -> Path:
    atlas_root.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        {
            "Association_ID": [f"EA{i:08d}" for i in range(n_rows)],
            "probe_ID": [f"cg{i:08d}" for i in range(n_rows)],
            "trait": ["obesity"] * n_rows,
            "p_value": ["1e-12"] * n_rows,
            "correlation": ["pos"] * n_rows,
            "study_ID": [f"ES{i:05d}" for i in range(n_rows)],
            "PMID": ["12345"] * n_rows,
        }
    )
    path = atlas_root / ATLAS_ASSOCIATIONS_FILE
    frame.to_csv(path, sep="\t", index=False)
    return path


def test_load_atlas_associations_head_respects_nrows(tmp_path: Path) -> None:
    atlas_root = tmp_path / "ewas_atlas"
    _write_atlas_associations(atlas_root, n_rows=10)

    head = load_atlas_associations_head(atlas_root, nrows=3)
    assert len(head) == 3
    assert "probe_ID" in head.columns

    full = load_atlas_associations_head(atlas_root, nrows=None)
    assert len(full) == 10


def test_load_atlas_associations_head_missing_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_atlas_associations_head(tmp_path / "absent")


def _explicit_edges() -> pd.DataFrame:
    # cg01 -> explicit body edge (keep); cg02 -> nearest-gene edge (reject);
    # cg03 -> non-explicit edge (reject).
    return pd.DataFrame(
        {
            "probe_id": ["cg01", "cg02", "cg03"],
            "gene_id": ["ENSG_A", "ENSG_B", "ENSG_C"],
            "gene_role": ["body", "body", "promoter"],
            "is_explicit": [True, True, False],
            "is_nearest_gene": [False, True, False],
            "mapping_source": ["graph", "nearest", "graph"],
        }
    )


def test_remap_drops_nearest_and_nonexplicit_edges() -> None:
    associations = pd.DataFrame(
        {
            "probe_id": ["cg01", "cg02", "cg03"],
            "study_ID": ["ES1", "ES1", "ES1"],
            # Atlas gene symbol is metadata only; must NOT be used to allocate.
            "gene": ["WRONG_A", "WRONG_B", "WRONG_C"],
            "p_value": ["1e-9", "1e-9", "1e-9"],
        }
    )
    remapped = remap_associations_to_graph(associations, locus_gene_edges=_explicit_edges())

    # Only the explicit, non-nearest cg01 edge survives.
    assert remapped["probe_id"].tolist() == ["cg01"]
    assert remapped["gene_id"].tolist() == ["ENSG_A"]
    # Atlas symbol retained as metadata, but the graph gene_id wins allocation.
    assert remapped["atlas_gene_symbol"].tolist() == ["WRONG_A"]
    assert remapped["gene_id"].iloc[0] != remapped["atlas_gene_symbol"].iloc[0]


def test_remap_requires_edge_safety_columns() -> None:
    associations = pd.DataFrame({"probe_id": ["cg01"], "study_ID": ["ES1"]})
    bad_edges = pd.DataFrame({"probe_id": ["cg01"], "gene_id": ["ENSG_A"]})
    with pytest.raises(ValueError, match="is_explicit"):
        remap_associations_to_graph(associations, locus_gene_edges=bad_edges)


def _atlas_gse_map() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gse_id": ["GSE100", "GSE200", "GSE300"],
            "atlas_study_id": ["ES00100", "ES00200;ES00201", "ES00300"],
            "pmid": ["111", "222", "333"],
        }
    )


def test_external_clean_excludes_benchmark_overlap_via_map() -> None:
    candidates = ["ES00100", "ES00200", "ES00201", "ES00300", "ES00999"]
    clean = external_clean_study_ids(
        candidates,
        benchmark_gse_ids={"GSE200"},  # contaminates ES00200 + ES00201 (same row)
        atlas_gse_map=_atlas_gse_map(),
    )
    assert "ES00200" not in clean
    assert "ES00201" not in clean
    assert clean == ["ES00100", "ES00300", "ES00999"]


def test_external_clean_never_matches_es_to_gse_string_equality() -> None:
    # Benchmark literally names an ES id; string equality must NOT drop it,
    # and an unrelated GSE benchmark must not remove anything.
    candidates = ["ES00100", "ES00300"]
    clean = external_clean_study_ids(
        candidates,
        benchmark_gse_ids={"ES00100", "GSE999"},
        atlas_gse_map=_atlas_gse_map(),
    )
    assert clean == ["ES00100", "ES00300"]


def test_external_clean_excludes_via_shared_pmid_bridge() -> None:
    # pmid_to_gse links PMID 333 to a benchmark GSE; ES00300 shares that PMID.
    clean = external_clean_study_ids(
        ["ES00100", "ES00300"],
        benchmark_gse_ids={"GSE_BENCH"},
        atlas_gse_map=_atlas_gse_map(),
        pmid_to_gse={"333": ["GSE_BENCH"]},
    )
    assert clean == ["ES00100"]


def test_build_external_clean_gene_list_ranks_and_filters() -> None:
    remapped = pd.DataFrame(
        {
            "probe_id": ["cg1", "cg2", "cg3", "cg4", "cg5"],
            "gene_id": ["ENSG_A", "ENSG_A", "ENSG_B", "ENSG_A", "ENSG_B"],
            "study_ID": ["ESa", "ESb", "ESa", "ESa", "ES_DIRTY"],
            "p_value": ["1e-9", "1e-8", "1e-7", "1e-6", "1e-9"],
            "correlation": ["pos", "pos", "neg", "pos", "pos"],
        }
    )
    panel = build_external_clean_gene_list(
        remapped,
        clean_study_ids=["ESa", "ESb"],  # ES_DIRTY excluded
    )

    # ENSG_A: 3 cpgs across 2 studies; ENSG_B: 1 cpg (cg3) after dropping dirty cg5.
    assert panel["gene_id"].tolist() == ["ENSG_A", "ENSG_B"]
    assert panel["rank"].tolist() == [1, 2]
    row_a = panel.loc[panel["gene_id"] == "ENSG_A"].iloc[0]
    assert row_a["n_independent_studies"] == 2
    assert row_a["n_associated_cpgs"] == 3
    assert row_a["direction_consistency"] == pytest.approx(1.0)
    assert row_a["inclusion_reason"] == "external_clean"


def test_build_hybrid_fold_gene_list_unions_and_caps() -> None:
    external = pd.DataFrame(
        {
            "gene_id": ["ENSG_A", "ENSG_B"],
            "score": [10.0, 4.0],
            "n_associated_cpgs": [3, 1],
        }
    )
    internal = pd.DataFrame(
        {
            "gene_id": ["ENSG_B", "ENSG_C"],
            "score": [8.0, 7.0],
            "n_associated_cpgs": [2, 5],
        }
    )
    hybrid = build_hybrid_fold_gene_list(
        external_clean=external, internal_fold=internal, n_genes=2
    )
    assert hybrid["gene_id"].tolist() == ["ENSG_A", "ENSG_B"]
    assert hybrid["inclusion_reason"].tolist() == ["hybrid_fold+external", "hybrid_fold"]
    assert len(hybrid) == 2


def test_sql_013_ddl_applies_and_round_trips(tmp_path: Path) -> None:
    sql = Path(__file__).resolve().parents[2] / "sql" / "013_association_seed_panels.sql"
    con = duckdb.connect(str(tmp_path / "assoc.duckdb"))
    try:
        con.execute(sql.read_text(encoding="utf-8"))
        con.execute(
            "INSERT INTO trait_ontology (trait_id, raw_trait_name, value_type, usage_class) "
            "VALUES ('T1', 'obesity', 'binary', 'core')"
        )
        con.execute(
            "INSERT INTO seed_panel (seed_panel_id, trait_id, selection_source) "
            "VALUES ('SP1', 'T1', 'external_clean')"
        )
        con.execute(
            "INSERT INTO seed_panel_gene (seed_panel_id, gene_id, rank, score) "
            "VALUES ('SP1', 'ENSG_A', 1, 6.0)"
        )
        n = con.execute(
            "SELECT count(*) FROM seed_panel_gene WHERE seed_panel_id = 'SP1'"
        ).fetchone()
        assert n is not None and n[0] == 1
        # Intra-track FK is enforced: unknown seed_panel_id must fail.
        with pytest.raises(duckdb.ConstraintException):
            con.execute(
                "INSERT INTO seed_panel_gene (seed_panel_id, gene_id) VALUES ('MISSING', 'ENSG_Z')"
            )
    finally:
        con.close()
