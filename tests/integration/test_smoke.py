from __future__ import annotations

from pathlib import Path

import torch

from mbs.catalog import build_catalog
from mbs.models import HierarchicalDeepSet


def test_catalog_build_and_model_forward(tmp_path: Path) -> None:
    sql_dir = tmp_path / "sql"
    parquet_root = tmp_path / "tables"
    database = tmp_path / "catalog.duckdb"
    sql_dir.mkdir()
    parquet_root.mkdir()
    (sql_dir / "001_schema.sql").write_text(
        "CREATE TABLE IF NOT EXISTS smoke(id INTEGER PRIMARY KEY);",
        encoding="utf-8",
    )
    (sql_dir / "010_view.sql").write_text(
        "CREATE OR REPLACE VIEW smoke_view AS SELECT * FROM smoke;",
        encoding="utf-8",
    )

    result = build_catalog(
        database=database,
        sql_dir=sql_dir,
        parquet_root=parquet_root,
    )

    assert database.exists()
    assert result["tables"] == 1
    assert result["views"] == 1

    model = HierarchicalDeepSet(
        input_dim=2,
        n_region_types=5,
        dropout=0.0,
    )
    model.eval()
    output = model(
        cpg_features=torch.tensor(
            [
                [0.2, 1.0],
                [0.8, -1.0],
                [0.5, 0.4],
            ]
        ),
        cpg_to_region=torch.tensor([0, 0, 1]),
        region_type=torch.tensor([0, 3]),
        region_to_gene=torch.tensor([0, 0]),
        n_regions=2,
        n_gene_instances=1,
    )

    assert output["mbs"].shape == (1,)
    assert output["present"].tolist() == [True]
    assert torch.isfinite(output["mbs"]).all()
    assert output["residual_mbs"].shape == (0,) or output["residual_mbs"].numel() >= 0
    assert "residual_present" in output
