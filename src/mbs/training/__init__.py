"""Flat DeepRVAT-style baseline training (Milestone 5)."""

from __future__ import annotations

from mbs.training.dataset import FlatBatch, FlatSampleRecord, build_flat_sample
from mbs.training.features import beta_to_m_value, gather_sample_features
from mbs.training.locus_gene import LocusGeneIndex, build_locus_gene_index
from mbs.training.loop import TrainResult, resolve_device, train_flat_baseline
from mbs.training.phenotypes import SamplePhenotype, load_gse35069_phenotypes

__all__ = [
    "FlatBatch",
    "FlatSampleRecord",
    "LocusGeneIndex",
    "SamplePhenotype",
    "TrainResult",
    "beta_to_m_value",
    "build_flat_sample",
    "build_locus_gene_index",
    "gather_sample_features",
    "load_gse35069_phenotypes",
    "resolve_device",
    "train_flat_baseline",
]
