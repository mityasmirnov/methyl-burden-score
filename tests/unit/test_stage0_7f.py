"""Milestone 7F: leftover→direct, RBS→gene allocation, no TBS fusion."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mbs.models import CascadeDeepSet
from mbs.training.cascade_assign import (
    ORPHAN_GENE_INDEX,
    assignment_gene_linked_only,
    build_cascade_assignment,
    gene_linked_col_index,
    nearest_gene_on_chromosome,
)
from mbs.training.cascade_loop import (
    make_synthetic_cascade_tables,
    run_cascade_fixture,
    score_samples,
    train_cascade_on_arrays,
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


def test_gene_allocation_explicit_only_excludes_nearest_gene_rbs() -> None:
    tables = make_synthetic_cascade_tables(seed=0)
    legacy = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="legacy_nearest",
    )
    explicit = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    legacy_cols = set(gene_linked_col_index(legacy).tolist())
    explicit_cols = set(gene_linked_col_index(explicit).tolist())
    assert legacy_cols.issuperset(explicit_cols)
    assert legacy_cols != explicit_cols
    cgi_i = legacy.region_ids.index("RBS:cgi_1")
    assert legacy.region_to_gene[cgi_i] >= 0
    explicit_cgi_i = explicit.region_ids.index("RBS:cgi_1")
    assert explicit.region_to_gene[explicit_cgi_i] == ORPHAN_GENE_INDEX


def test_gene_allocation_bounded_nearest_distance_gate() -> None:
    tables = make_synthetic_cascade_tables(seed=0)
    regions = pd.concat(
        [
            tables["regions"],
            pd.DataFrame(
                [
                    {
                        "region_id": "RBS:far_cgi",
                        "gene_id": None,
                        "region_type": "cgi_island",
                        "region_system": "rbs",
                        "chromosome": "chr1",
                        "start": 2000,
                        "end": 2100,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    edges = pd.concat(
        [
            tables["locus_region_edges"],
            pd.DataFrame([{"locus_id": 40, "region_id": "RBS:far_cgi"}]),
        ],
        ignore_index=True,
    )
    kwargs = {
        "locus_index": tables["locus_index"],
        "locus_region_edges": edges,
        "regions": regions,
        "genes": tables["genes"],
        "gene_allocation": "bounded_nearest",
    }
    near = build_cascade_assignment(**kwargs, max_nearest_gene_bp=5000)
    far = build_cascade_assignment(**kwargs, max_nearest_gene_bp=500)
    far_i = near.region_ids.index("RBS:far_cgi")
    assert near.region_to_gene[far_i] >= 0
    assert far.region_to_gene[far_i] == ORPHAN_GENE_INDEX


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
    mbs, present, orphan, all_rbs, all_present = score_samples(
        model, assignment, tables["betas"][:2], device=torch.device("cpu")
    )
    assert mbs.shape == (2, assignment.n_genes)
    assert present.shape == mbs.shape
    assert orphan.shape == (2, assignment.n_orphan_rbs)
    assert all_rbs.shape == (2, assignment.n_regions)
    assert all_present.shape == all_rbs.shape
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


def test_assignment_gene_linked_only_drops_direct() -> None:
    tables = make_synthetic_cascade_tables(seed=1)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    assert assignment.n_direct >= 1
    linked = assignment_gene_linked_only(assignment)
    assert linked.direct_col_index.size == 0
    assert linked.n_direct == 0


def test_cascade_writes_direct_cpg_when_locus_ids_provided(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=2)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    assert assignment.n_direct >= 1
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, (n * 2) // 3), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()
    locus_ids = tables["locus_index"]["locus_id"].astype(str).tolist()
    out_dir = tmp_path / "cascade_direct_cpg"
    metrics = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=out_dir,
        max_epochs=2,
        seed=2,
        device_str="cpu",
        cpg_hidden_dim=16,
        region_hidden_dim=8,
        dropout=0.0,
        locus_ids=locus_ids,
        include_mbs_enet=False,
    )
    score_dir = Path(str(metrics["score_dir"]))
    assert (score_dir / "direct_cpg.zarr").exists()
    assert (score_dir / "direct_locus_index.parquet").is_file()
    manifest = __import__("json").loads(
        (score_dir / "score_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest.get("direct_cpg") is True
    assert int(manifest.get("n_direct_loci") or 0) == int(assignment.n_direct)
