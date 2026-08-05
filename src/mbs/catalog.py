"""DuckDB catalog construction utilities."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

import duckdb

if TYPE_CHECKING:
    from mbs.paths import DataPaths


class CatalogBuildResult(TypedDict):
    database: str
    read_only: bool
    executed_sql: list[str]
    tables: int
    views: int


def _render_sql(sql: str, parquet_root: Path) -> str:
    escaped_root = parquet_root.absolute().as_posix().replace("'", "''")
    return sql.replace("${PARQUET_ROOT}", escaped_root)


def build_catalog(
    *,
    database: Path,
    sql_dir: Path,
    parquet_root: Path,
    read_only: bool = False,
) -> CatalogBuildResult:
    """Create or inspect the analytical catalog.

    Numbered ``*.sql`` files are executed in lexical order. SQL files may use
    ``${PARQUET_ROOT}``, which is replaced with an escaped absolute path.
    """
    database = database.absolute()
    sql_dir = sql_dir.absolute()
    parquet_root = parquet_root.absolute()

    if not sql_dir.is_dir():
        raise FileNotFoundError(f"SQL directory does not exist: {sql_dir}")
    if not parquet_root.exists():
        raise FileNotFoundError(f"Parquet root does not exist: {parquet_root}")

    database.parent.mkdir(parents=True, exist_ok=True)
    if read_only and not database.exists():
        raise FileNotFoundError(f"Read-only catalog does not exist: {database}")

    sql_files = sorted(sql_dir.glob("[0-9][0-9][0-9]_*.sql"))
    if not sql_files:
        raise FileNotFoundError(f"No numbered SQL files found under {sql_dir}")

    executed: list[str] = []
    connection = duckdb.connect(str(database), read_only=read_only)
    try:
        if not read_only:
            connection.execute("BEGIN TRANSACTION")
            try:
                for sql_file in sql_files:
                    rendered = _render_sql(sql_file.read_text(encoding="utf-8"), parquet_root)
                    connection.execute(rendered)
                    executed.append(sql_file.name)
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise

        relation_counts = connection.execute(
            """
            SELECT table_type, count(*)
            FROM information_schema.tables
            WHERE table_schema = 'main'
            GROUP BY table_type
            """
        ).fetchall()
        counts = {str(table_type): int(count) for table_type, count in relation_counts}
    finally:
        connection.close()

    return {
        "database": str(database),
        "read_only": read_only,
        "executed_sql": executed,
        "tables": counts.get("BASE TABLE", 0),
        "views": counts.get("VIEW", 0),
    }


def init_catalog(
    *,
    paths: DataPaths,
    sql_dir: Path | None = None,
    database: Path | None = None,
    parquet_root: Path | None = None,
) -> CatalogBuildResult:
    """Ensure Stage 0 directories exist and apply numbered SQL under defaults."""
    paths.ensure_directories()
    resolved_sql_dir = (sql_dir or (paths.project_root / "sql")).absolute()
    resolved_database = (
        database or (paths.data_root / "canonical" / "catalog" / "catalog.duckdb")
    ).absolute()
    resolved_parquet_root = (
        parquet_root or (paths.data_root / "canonical" / "catalog" / "tables")
    ).absolute()
    resolved_parquet_root.mkdir(parents=True, exist_ok=True)
    return build_catalog(
        database=resolved_database,
        sql_dir=resolved_sql_dir,
        parquet_root=resolved_parquet_root,
        read_only=False,
    )
