from __future__ import annotations

from pathlib import Path
from uuid import uuid4

import duckdb

from mbs.catalog import build_catalog, init_catalog
from mbs.paths import DataPaths

# sql/001_schema.sql + sql/010_views.sql (keep in sync when schema changes)
EXPECTED_TABLES = 18
EXPECTED_VIEWS = 8


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def test_build_catalog_applies_project_sql(tmp_path: Path) -> None:
    repo = _repo_root()
    sql_dir = repo / "sql"
    parquet_root = tmp_path / "tables"
    database = tmp_path / "catalog.duckdb"
    parquet_root.mkdir()

    result = build_catalog(
        database=database,
        sql_dir=sql_dir,
        parquet_root=parquet_root,
    )

    assert database.exists()
    assert result["executed_sql"] == ["001_schema.sql", "010_views.sql"]
    assert result["tables"] == EXPECTED_TABLES
    assert result["views"] == EXPECTED_VIEWS

    connection = duckdb.connect(str(database), read_only=True)
    try:
        tables = {
            row[0]
            for row in connection.execute(
                """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = 'main' AND table_type = 'BASE TABLE'
                """
            ).fetchall()
        }
    finally:
        connection.close()

    assert "source_release" in tables
    assert "sample" in tables
    assert "fold_assignment" in tables


def test_init_catalog_creates_dirs_and_schema() -> None:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"catalog-init-{uuid4().hex}"
    workspace.mkdir()

    paths = DataPaths(
        project_root=repo,
        data_root=workspace / "data",
        scratch_root=workspace / "scratch",
        cache_root=workspace / "cache",
        artifact_root=workspace / "artifacts",
        docker_root=workspace / "docker",
    )

    result = init_catalog(paths=paths)

    assert Path(result["database"]).exists()
    assert (workspace / "data" / "canonical" / "catalog" / "tables").is_dir()
    assert (workspace / "data" / "raw" / "cpgcorpus").is_dir()
    assert result["tables"] == EXPECTED_TABLES
    assert result["views"] == EXPECTED_VIEWS
