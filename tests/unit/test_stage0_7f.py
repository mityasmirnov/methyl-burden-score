"""Milestone 7F: leftover→direct, RBS→gene allocation, no TBS fusion."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mbs.models import CascadeDeepSet
from mbs.training.cascade_assign import (
    ORPHAN_GENE_INDEX,
    build_cascade_assignment,
    nearest_gene_on_chromosome,
)
from mbs.training.cascade_loop import (
    make_synthetic_cascade_tables,
    run_cascade_fixture,
    score_samples,
)
from mbs.training.cascade_scores import (
    fusion_feature_matrix,
    load_cascade_score_blocks,
    write_cascade_score_dir,
)


def test_nearest_gene_and_orphan_chromosome() -> None:
    genes = pd.DataFrame(
        {
            "gene_id": ["ENSG1"],
            "chromosome": ["chr1"],
            "start": [100],
            "end": [200],
        }
    )
    assert nearest_gene_on_chromosome("chr1", 150, genes) == "ENSG1"
    assert nearest_gene_on_chromosome("chr2", 150, genes) is None


def test_cascade_assign_leftover_direct_and_rbs_gene() -> None:
    tables = make_synthetic_cascade_tables(seed=0)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    # TBS locus col 4 and leftover col 3 → direct; never nearest-gene as CpGs.
    assert set(assignment.direct_col_index.tolist()) == {3, 4}
    # Gene promoter allocated to ENSG1; CGI near ENSG1 allocated; chr2 RBS orphan.
    assert "ENSG1" in assignment.gene_ids
    assert assignment.n_orphan_rbs == 1
    assert assignment.orphan_region_ids == ["RBS:orphan_island"]
    # Nearest-gene for RBS:cgi_1 → ENSG1
    cgi_i = assignment.region_ids.index("RBS:cgi_1")
    assert assignment.allocated_gene_id[cgi_i] == "ENSG1"
    assert assignment.region_to_gene[cgi_i] >= 0
    orphan_i = assignment.region_ids.index("RBS:orphan_island")
    assert assignment.region_to_gene[orphan_i] == ORPHAN_GENE_INDEX
    # No tbs region in panel
    assert all("TILE" not in rid for rid in assignment.region_ids)


def test_cascade_deepset_forward_and_empty_gene_mask() -> None:
    tables = make_synthetic_cascade_tables(seed=1)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    model = CascadeDeepSet(1, len(assignment.region_types), cpg_hidden_dim=8, region_hidden_dim=4)
    mbs, present, orphan = score_samples(
        model, assignment, tables["betas"][:2], device=torch.device("cpu")
    )
    assert mbs.shape == (2, assignment.n_genes)
    assert present.shape == mbs.shape
    assert orphan.shape == (2, assignment.n_orphan_rbs)
    # Absent genes stay neutral when no active regions (ENSG2 unused).
    ens2 = assignment.gene_ids.index("ENSG2") if "ENSG2" in assignment.gene_ids else None
    if ens2 is not None:
        assert not bool(present[0, ens2])
        assert abs(float(mbs[0, ens2]) - 0.5) < 1e-5


def test_fusion_matrix_rejects_tbs_and_writes_scores(tmp_path: Path) -> None:
    sample_ids = ["a", "b"]
    gene_ids = ["G1"]
    orphan_ids = ["R1"]
    write_cascade_score_dir(
        tmp_path,
        sample_ids=sample_ids,
        gene_ids=gene_ids,
        orphan_region_ids=orphan_ids,
        mbs=np.ones((2, 1), dtype=np.float32),
        gene_present=np.ones((2, 1), dtype=bool),
        orphan_rbs=np.full((2, 1), 0.7, dtype=np.float32),
        direct_contrib=np.zeros((2, 1), dtype=np.float32),
        direct_task_names=["age"],
    )
    assert not (tmp_path / "tbs.zarr").exists()
    blocks = load_cascade_score_blocks(tmp_path)
    x = fusion_feature_matrix(blocks)
    assert x.shape == (2, 3)
    try:
        fusion_feature_matrix({**blocks, "tbs": np.zeros((2, 1), dtype=np.float32)})
        raise AssertionError("expected TBS rejection")
    except ValueError as exc:
        assert "TBS" in str(exc)


def test_train_cascade_fixture_writes_report(tmp_path: Path) -> None:
    artifact_root = tmp_path / "artifacts"
    project_root = tmp_path / "proj"
    (project_root / "reports" / "inspection").mkdir(parents=True)
    artifact_root.mkdir(parents=True)
    result = run_cascade_fixture(
        project_root=project_root,
        artifact_root=artifact_root,
        run_id="7f-test",
        seed=0,
        max_epochs=3,
        device_str="cpu",
    )
    assert result.report_dir.is_dir()
    assert (result.report_dir / "summary.json").is_file()
    assert (result.report_dir / "analysis.md").is_file()
    assert (result.score_dir / "mbs.zarr").exists()
    assert (result.score_dir / "rbs.zarr").exists()
    assert (result.score_dir / "direct_contrib.zarr").exists()
    assert not (result.score_dir / "tbs.zarr").exists()
    summary = result.metrics
    assert summary["tbs_arm"] is False
    assert summary["assignment"]["n_direct"] >= 1
    assert summary["assignment"]["n_orphan_rbs"] >= 1
