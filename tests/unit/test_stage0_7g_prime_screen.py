"""Stage A DeepRVAT screen: pooling, vector cascade, one-hop features, RBS export."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mbs.models import CascadeDeepSet, FlatDeepSetRegion
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import make_synthetic_cascade_tables, score_samples
from mbs.training.cascade_scores import load_cascade_score_blocks, write_cascade_score_dir
from mbs.training.flat_region_features import (
    GENE_ROLES,
    assert_other_gene_count,
    build_flat_region_gene_index,
    count_other_gene_edges,
    flat_region_input_dim,
    gather_flat_region_features,
)
from mbs.training.flat_region_loop import train_flat_region_on_arrays


def test_independent_pooling_and_vector_aggregation() -> None:
    tables = make_synthetic_cascade_tables(seed=3)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    n_types = max(len(assignment.region_types), 1)
    scalar = CascadeDeepSet(
        1,
        n_types,
        cpg_pool="mean",
        region_pool="max",
        gene_aggregation="scalar_rbs",
        cpg_hidden_dim=8,
        region_hidden_dim=4,
        dropout=0.0,
    )
    vector = CascadeDeepSet(
        1,
        n_types,
        cpg_pool="mean",
        region_pool="max",
        gene_aggregation="region_hidden",
        cpg_hidden_dim=8,
        region_hidden_dim=4,
        dropout=0.0,
    )
    assert scalar.cpg_pool == "mean"
    assert scalar.region_pool == "max"
    assert vector.gene_aggregation == "region_hidden"
    mbs_s, _, _, rbs_s, _ = score_samples(
        scalar, assignment, tables["betas"][:1], device=torch.device("cpu")
    )
    mbs_v, _, _, rbs_v, _ = score_samples(
        vector, assignment, tables["betas"][:1], device=torch.device("cpu")
    )
    assert mbs_s.shape[1] == assignment.n_genes
    assert mbs_v.shape[1] == assignment.n_genes
    assert rbs_s.shape[1] == assignment.n_regions
    assert rbs_v.shape[1] == assignment.n_regions


def test_vector_permutation_invariance_within_region() -> None:
    # Region-type embeddings are randomly initialized; seed so the "region_type
    # change alters output" assertion below doesn't depend on leftover global
    # RNG state from whichever test ran before this one in the same process.
    torch.manual_seed(0)
    model = CascadeDeepSet(
        1,
        2,
        cpg_pool="mean",
        region_pool="max",
        gene_aggregation="region_hidden",
        cpg_hidden_dim=8,
        region_hidden_dim=4,
        dropout=0.0,
    )
    model.eval()
    cpg_features = torch.tensor([[0.2], [0.8], [0.5]], dtype=torch.float32)
    cpg_to_region = torch.tensor([0, 0, 1], dtype=torch.long)
    region_type = torch.tensor([0, 1], dtype=torch.long)
    region_to_gene = torch.tensor([0, 0], dtype=torch.long)
    out_a = model(
        cpg_features=cpg_features,
        cpg_to_region=cpg_to_region,
        region_type=region_type,
        region_to_gene=region_to_gene,
        n_regions=2,
        n_gene_instances=1,
    )
    perm = torch.tensor([1, 0, 2], dtype=torch.long)
    out_b = model(
        cpg_features=cpg_features[perm],
        cpg_to_region=cpg_to_region[perm],
        region_type=region_type,
        region_to_gene=region_to_gene,
        n_regions=2,
        n_gene_instances=1,
    )
    assert torch.allclose(out_a["mbs"], out_b["mbs"], atol=1e-5)
    # Region type change can alter output.
    region_type_alt = torch.tensor([1, 0], dtype=torch.long)
    out_c = model(
        cpg_features=cpg_features,
        cpg_to_region=cpg_to_region,
        region_type=region_type_alt,
        region_to_gene=region_to_gene,
        n_regions=2,
        n_gene_instances=1,
    )
    assert not torch.allclose(out_a["mbs"], out_c["mbs"], atol=1e-5)


def test_flat_region_channels_and_other_gene_check() -> None:
    tables = make_synthetic_cascade_tables(seed=4)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    # Annotate CGI context on the synthetic locus index.
    locus_index = tables["locus_index"].copy()
    n_loci = len(locus_index)
    ctx = ["island", "open_sea", "north_shore", "open_sea", "unknown"]
    locus_index["cpg_context"] = [ctx[i % len(ctx)] for i in range(n_loci)]
    index = build_flat_region_gene_index(
        assignment,
        locus_index=locus_index,
        allow_other_gene=True,
    )
    dim = flat_region_input_dim()
    assert dim == 1 + len(GENE_ROLES) + 7 + 6 + 3 + 1
    feats, cpg_to_gene = gather_flat_region_features(
        beta_row=np.asarray(tables["betas"][0], dtype=np.float32),
        index=index,
    )
    assert feats.shape == (feats.shape[0], dim)
    assert cpg_to_gene.shape[0] == feats.shape[0]
    assert feats.shape[0] <= index.n_edges
    # Presence flags live near the end (before observed).
    assert feats[:, -1].min() >= 0.0
    model = FlatDeepSetRegion(dim, phi_hidden_dim=8, rho_hidden_dim=4, phi_layers=1, rho_layers=1)
    out = model(torch.from_numpy(feats), torch.from_numpy(cpg_to_gene), index.n_genes)
    assert out["mbs"].shape == (index.n_genes,)

    # other_gene count on five-role types should be zero for gene roles only.
    gene_edge = assignment.region_to_gene[assignment.edge_region_index] >= 0
    type_ids = assignment.region_type_id[assignment.edge_region_index[gene_edge]]
    # Filter to gene-system types only for the assert.
    five = {"promoter_core", "promoter_proximal", "five_prime", "three_prime", "gene_body"}
    only_five = np.asarray(
        [
            tid
            for tid in type_ids.tolist()
            if 0 <= tid < len(assignment.region_types)
            and assignment.region_types[int(tid)] in five
        ],
        dtype=np.int64,
    )
    n_other = count_other_gene_edges(only_five, assignment.region_types)
    assert n_other == 0
    assert_other_gene_count(n_other_gene_edges=0)
    try:
        assert_other_gene_count(n_other_gene_edges=1)
        raise AssertionError("expected AssertionError")
    except AssertionError as exc:
        assert "other_gene" in str(exc)


def test_one_hop_permutation_invariance() -> None:
    dim = flat_region_input_dim()
    model = FlatDeepSetRegion(dim, phi_hidden_dim=8, rho_hidden_dim=4, phi_layers=1, rho_layers=1)
    model.eval()
    feats = torch.randn(4, dim)
    genes = torch.tensor([0, 0, 1, 1], dtype=torch.long)
    out_a = model(feats, genes, 2)
    perm = torch.tensor([1, 0, 3, 2], dtype=torch.long)
    out_b = model(feats[perm], genes[perm], 2)
    assert torch.allclose(out_a["mbs"], out_b["mbs"], atol=1e-5)


def test_all_gene_rbs_export_and_index(tmp_path: Path) -> None:
    write_cascade_score_dir(
        tmp_path,
        sample_ids=["a", "b"],
        gene_ids=["G1"],
        orphan_region_ids=[],
        mbs=np.ones((2, 1), dtype=np.float32),
        gene_present=np.ones((2, 1), dtype=bool),
        orphan_rbs=np.zeros((2, 0), dtype=np.float32),
        direct_contrib=np.zeros((2, 0), dtype=np.float32),
        direct_task_names=[],
        all_gene_rbs=np.full((2, 2), 0.6, dtype=np.float32),
        all_gene_rbs_present=np.ones((2, 2), dtype=bool),
        all_gene_region_ids=["R1", "R2"],
        all_gene_region_gene_ids=["G1", "G1"],
        all_gene_region_types=["promoter_core", "gene_body"],
        allocation_policy="explicit_only",
        extra_manifest={"cpg_pool": "mean", "region_pool": "max"},
    )
    assert (tmp_path / "all_gene_rbs.zarr").exists()
    assert (tmp_path / "all_gene_rbs_present.zarr").exists()
    idx = pd.read_parquet(tmp_path / "all_gene_region_index.parquet")
    assert list(idx.columns) == [
        "region_id",
        "gene_id",
        "region_type",
        "column_index",
        "allocation_policy",
    ]
    assert len(idx) == 2
    blocks = load_cascade_score_blocks(tmp_path)
    assert "all_gene_rbs" in blocks
    assert blocks["all_gene_rbs"].shape == (2, 2)
    # Orphan production rbs.zarr remains empty / separate.
    assert blocks["orphan_rbs"].shape == (2, 0)


def test_flat_region_loop_reports_all_tasks(tmp_path: Path) -> None:
    tables = make_synthetic_cascade_tables(seed=5)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, n - 1), dtype=np.int64)
    test_idx = np.asarray([n - 1], dtype=np.int64)
    payload = train_flat_region_on_arrays(
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
        out_dir=tmp_path / "fold",
        max_epochs=1,
        device_str="cpu",
        age_mask=np.ones(n, dtype=bool),
        tissue_mask=np.ones(n, dtype=bool),
        sex_mask=np.ones(n, dtype=bool),
        allow_other_gene=True,
        pool="max",
        arm="N-light-gene-max",
        locus_index=tables["locus_index"],
    )
    metrics = payload["metrics"]
    assert "tissue" in metrics
    assert "age" in metrics
    assert "sex" in metrics
    assert payload["eval_split"] == "test"
