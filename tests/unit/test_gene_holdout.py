"""GATE G1 gene-holdout partition + index/assignment subset."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import make_synthetic_cascade_tables
from mbs.training.flat_region_features import REGULATORY_CHANNELS, FlatRegionGeneIndex
from mbs.training.gene_holdout import (
    gene_chromosome_map,
    partition_genes_by_chromosome,
    partition_genes_random,
    resolve_gene_holdout_partition,
    subset_cascade_assignment,
    subset_flat_region_gene_index,
    write_gene_holdout_artifacts,
)


def _toy_index(*, n_genes: int = 3, edges_per_gene: int = 2) -> FlatRegionGeneIndex:
    gene_ids = [f"G{i}" for i in range(n_genes)]
    n_edges = n_genes * edges_per_gene
    gene_ix = np.repeat(np.arange(n_genes, dtype=np.int64), edges_per_gene)
    return FlatRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=np.arange(n_edges, dtype=np.int64),
        edge_gene_index=gene_ix,
        edge_role_id=np.zeros(n_edges, dtype=np.int64),
        edge_context_id=np.zeros(n_edges, dtype=np.int64),
        edge_role_present=np.ones(n_edges, dtype=bool),
        edge_context_present=np.ones(n_edges, dtype=bool),
        edge_regulatory_present=np.zeros(n_edges, dtype=bool),
        edge_regulatory_multi_hot=np.zeros(
            (n_edges, len(REGULATORY_CHANNELS)), dtype=np.float32
        ),
        n_study_loci=n_edges,
        n_other_gene_edges=0,
    )


def test_partition_genes_random_disjoint_and_reproducible() -> None:
    genes = [f"G{i}" for i in range(50)]
    a = partition_genes_random(genes, heldout_fraction=0.2, seed=42)
    b = partition_genes_random(genes, heldout_fraction=0.2, seed=42)
    assert a.train_gene_ids == b.train_gene_ids
    assert a.heldout_gene_ids == b.heldout_gene_ids
    assert set(a.train_gene_ids).isdisjoint(a.heldout_gene_ids)
    assert set(a.train_gene_ids) | set(a.heldout_gene_ids) == set(genes)
    assert a.n_heldout == 10
    assert a.n_train == 40


def test_partition_rejects_bad_fraction() -> None:
    with pytest.raises(ValueError, match="heldout_fraction"):
        partition_genes_random(["A", "B"], heldout_fraction=0.0)


def test_partition_genes_by_chromosome_holds_out_whole_chroms() -> None:
    genes = [f"G{i}" for i in range(20)]
    chrom = {g: f"chr{(i % 5) + 1}" for i, g in enumerate(genes)}
    part = partition_genes_by_chromosome(
        genes, chrom, heldout_fraction=0.25, seed=0
    )
    assert part.method == "chromosome"
    assert part.heldout_chromosomes
    assert set(part.train_gene_ids).isdisjoint(part.heldout_gene_ids)
    # Every heldout gene's chromosome is in the heldout chrom set.
    for g in part.heldout_gene_ids:
        assert chrom[g] in part.heldout_chromosomes
    # No train gene shares a heldout chromosome.
    for g in part.train_gene_ids:
        if g in chrom:
            assert chrom[g] not in set(part.heldout_chromosomes)


def test_partition_chromosome_explicit_list() -> None:
    genes = ["A", "B", "C", "D"]
    chrom = {"A": "chr1", "B": "chr1", "C": "chr2", "D": "chr3"}
    part = partition_genes_by_chromosome(
        genes, chrom, heldout_fraction=0.5, seed=0, heldout_chromosomes=["chr1"]
    )
    assert set(part.heldout_gene_ids) == {"A", "B"}
    assert set(part.train_gene_ids) == {"C", "D"}


def test_resolve_dispatches() -> None:
    genes = [f"G{i}" for i in range(10)]
    chrom = {g: "chr1" if i < 5 else "chr2" for i, g in enumerate(genes)}
    r = resolve_gene_holdout_partition(genes, method="random", seed=1)
    assert r.method == "random"
    c = resolve_gene_holdout_partition(
        genes, method="chromosome", seed=1, gene_chromosome=chrom
    )
    assert c.method == "chromosome"
    with pytest.raises(ValueError, match="unsupported"):
        resolve_gene_holdout_partition(genes, method="bogus")


def test_subset_flat_region_gene_index_remaps() -> None:
    index = _toy_index(n_genes=3, edges_per_gene=2)
    part = partition_genes_random(index.gene_ids, heldout_fraction=1 / 3, seed=1)
    train_idx = subset_flat_region_gene_index(index, part.train_gene_ids)
    held_idx = subset_flat_region_gene_index(index, part.heldout_gene_ids)

    assert set(train_idx.gene_ids).isdisjoint(held_idx.gene_ids)
    assert train_idx.n_edges + held_idx.n_edges == index.n_edges
    assert int(train_idx.edge_gene_index.min()) >= 0
    assert int(train_idx.edge_gene_index.max()) < train_idx.n_genes
    assert int(held_idx.edge_gene_index.max()) < held_idx.n_genes
    assert set(train_idx.edge_gene_index.tolist()) == set(range(train_idx.n_genes))
    assert set(held_idx.edge_gene_index.tolist()) == set(range(held_idx.n_genes))


def test_subset_cascade_assignment_remaps_genes() -> None:
    tables = make_synthetic_cascade_tables(seed=2)
    # Expand genes so both ENSG1 and ENSG2 have gene-linked regions.
    genes = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2"],
            "chromosome": ["chr1", "chr2"],
            "start": [100, 5000],
            "end": [400, 5500],
            "strand": ["+", "+"],
            "gene_name": ["G1", "G2"],
            "gene_type": ["protein_coding", "protein_coding"],
            "source_version": ["test", "test"],
        }
    )
    regions = tables["regions"].copy()
    # Attach second gene region.
    regions = pd.concat(
        [
            regions,
            pd.DataFrame(
                {
                    "region_id": ["ENSG2:promoter_core"],
                    "gene_id": ["ENSG2"],
                    "region_type": ["promoter_core"],
                    "region_system": ["gene"],
                    "chromosome": ["chr1"],
                    "start": [5000],
                    "end": [5100],
                }
            ),
        ],
        ignore_index=True,
    )
    locus_index = tables["locus_index"].copy()
    # Add a locus for ENSG2.
    locus_index = pd.concat(
        [
            locus_index,
            pd.DataFrame(
                {
                    "col_index": [5],
                    "locus_id": [60],
                    "canonical_key": ["GRCh38:chr1:5050"],
                }
            ),
        ],
        ignore_index=True,
    )
    lr_edges = pd.concat(
        [
            tables["locus_region_edges"],
            pd.DataFrame({"locus_id": [60], "region_id": ["ENSG2:promoter_core"]}),
        ],
        ignore_index=True,
    )
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        gene_allocation="explicit_only",
    )
    assert assignment.n_genes >= 2
    chrom_map = gene_chromosome_map(genes)
    assert chrom_map["ENSG1"] == "chr1"
    part = partition_genes_random(assignment.gene_ids, heldout_fraction=0.5, seed=3)
    train_asg = subset_cascade_assignment(
        assignment, part.train_gene_ids, keep_orphans=False, keep_direct=False
    )
    held_asg = subset_cascade_assignment(
        assignment, part.heldout_gene_ids, keep_orphans=False, keep_direct=False
    )
    assert set(train_asg.gene_ids).isdisjoint(held_asg.gene_ids)
    assert train_asg.n_genes >= 1
    assert held_asg.n_genes >= 1
    assert train_asg.n_direct == 0
    assert held_asg.n_direct == 0
    # Remapped gene indices are compact.
    if train_asg.edge_region_index.size:
        genes_on_edges = train_asg.region_to_gene[train_asg.edge_region_index]
        genes_on_edges = genes_on_edges[genes_on_edges >= 0]
        assert int(genes_on_edges.max()) < train_asg.n_genes


def test_write_gene_holdout_artifacts(tmp_path: Path) -> None:
    part = partition_genes_random([f"G{i}" for i in range(10)], heldout_fraction=0.3, seed=0)
    mbs = np.zeros((4, part.n_heldout), dtype=np.float32)
    present = np.ones((4, part.n_heldout), dtype=bool)
    path = write_gene_holdout_artifacts(
        tmp_path,
        partition=part,
        train_gene_ids=part.train_gene_ids,
        heldout_gene_ids=part.heldout_gene_ids,
        mbs_heldout=mbs,
        present_heldout=present,
    )
    assert path.is_file()
    assert (tmp_path / "mbs_heldout.npy").is_file()
    assert (tmp_path / "gene_ids_heldout.json").is_file()
