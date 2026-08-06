"""Unit tests for Milestone 5 flat baseline training helpers."""

from __future__ import annotations

from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import torch
from typer.testing import CliRunner

from mbs.cli import app
from mbs.paths import DataPaths
from mbs.training.dataset import make_synthetic_overfit_bundle, record_to_batch
from mbs.training.features import beta_to_m_value, gather_sample_features
from mbs.training.locus_gene import LocusGeneIndex, build_locus_gene_index
from mbs.training.loop import resolve_device, train_flat_baseline
from mbs.training.phenotypes import load_gse35069_phenotypes

runner = CliRunner()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _point_env(monkeypatch: pytest.MonkeyPatch, scratch: Path) -> None:
    repo = _repo_root()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(scratch / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(scratch / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(scratch / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(scratch / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(scratch / "docker"))


@pytest.fixture
def isolated_workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    scratch_base = _repo_root() / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    workspace = scratch_base / f"train-{uuid4().hex}"
    workspace.mkdir()
    _point_env(monkeypatch, workspace)
    return workspace


def test_beta_to_m_value_roundtrip_midpoint() -> None:
    beta = np.array([0.5], dtype=np.float32)
    m_value = beta_to_m_value(beta, epsilon=1e-3)
    assert abs(float(m_value[0])) < 1e-5


def test_build_locus_gene_index_collapses_same_gene_regions() -> None:
    locus_index = pd.DataFrame(
        {
            "col_index": [0, 1, 2],
            "locus_id": np.array([10, 11, 12], dtype=np.uint64),
        }
    )
    locus_region_edges = pd.DataFrame(
        {
            "locus_id": np.array([10, 10, 11], dtype=np.uint64),
            "region_id": ["g1:promoter", "g1:body", "g2:body"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": ["g1:promoter", "g1:body", "g2:body"],
            "gene_id": ["ENSG1", "ENSG1", "ENSG2"],
        }
    )
    index = build_locus_gene_index(
        locus_index=locus_index,
        locus_region_edges=locus_region_edges,
        regions=regions,
    )
    assert index.n_genes == 2
    assert index.n_edges == 2
    assert sorted(index.gene_ids) == ["ENSG1", "ENSG2"]


def test_gather_drops_nan_beta_and_missing_static() -> None:
    locus_gene = LocusGeneIndex(
        gene_ids=["G0", "G1"],
        edge_col_index=np.array([0, 1, 2], dtype=np.int64),
        edge_gene_index=np.array([0, 1, 0], dtype=np.int64),
        n_study_loci=3,
    )
    beta_row = np.array([0.2, np.nan, 0.8], dtype=np.float32)
    static = np.zeros((3, 2), dtype=np.float32)
    valid = np.array([True, True, False])
    bundle = gather_sample_features(
        beta_row=beta_row,
        static_by_col=static,
        static_valid=valid,
        locus_gene=locus_gene,
    )
    assert bundle.n_observed_edges == 1
    assert bundle.n_dropped_nan_beta == 1
    assert bundle.n_dropped_no_static == 1
    assert bundle.cpg_to_gene.tolist() == [0]


def test_phenotype_join_requires_gsm() -> None:
    scratch = _repo_root() / "scratch" / "pytest" / f"meta-{uuid4().hex}"
    scratch.mkdir(parents=True)
    meta_path = scratch / "metadata.arrow"
    table = pa.table(
        {
            "GSM_ID": ["GSM1", "GSM2"],
            "title": ["WB_1", "PBMC_2"],
            "tissue/cell type:ch1": ["Whole blood", "PBMC"],
        }
    )
    with (
        pa.OSFile(str(meta_path), "wb") as sink,
        pa.ipc.new_file(sink, table.schema) as writer,
    ):
        writer.write_table(table)

    phenotypes, classes = load_gse35069_phenotypes(meta_path, sample_ids=["GSM1", "GSM2"])
    assert classes == ["PBMC", "Whole blood"]
    assert phenotypes[0].donor_id == "1"
    assert phenotypes[1].class_index == classes.index("PBMC")

    with pytest.raises(KeyError):
        load_gse35069_phenotypes(meta_path, sample_ids=["GSM1", "GSM_MISSING"])


def test_overfit_fixture_reaches_perfect_accuracy(isolated_workspace: Path) -> None:
    paths = DataPaths.from_environment()
    paths.ensure_directories()
    cfg = {
        "experiment": {"name": "unit_overfit", "stage": 0, "seed": 0},
        "model": {
            "phi_layers": 2,
            "phi_hidden_dimension": 32,
            "rho_layers": 2,
            "rho_hidden_dimension": 16,
            "pooling": "max",
            "neutral_score": 0.5,
            "dropout": 0.0,
        },
        "training": {
            "optimizer": "adam",
            "learning_rate": 0.05,
            "weight_decay": 0.0,
            "max_epochs": 400,
            "early_stopping_patience": 100,
            "gradient_clip_norm": 2.0,
            "precision": "fp32",
            "require_cuda": False,
        },
        "heads": {"tissue": {"enabled": True}},
        "logging": {"tensorboard": True, "auto_tensorboard": False},
    }
    result = train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        run_id="unit-overfit",
        device_str="cpu",
        overfit_fixture=True,
        max_epochs=100,
    )
    assert result.metrics.get("overfit_ok") is True
    assert result.metrics.get("model_public_name") == "deepMAT"
    assert (result.run_dir / "resolved_config.yaml").is_file()
    assert (result.run_dir / "metrics.jsonl").is_file()
    assert (result.run_dir / "tb").is_dir()
    assert (result.checkpoint_dir / "best.pt").is_file()
    assert (result.checkpoint_dir / "checkpoint_manifest.json").is_file()


def test_resolve_device_cuda_requires_availability() -> None:
    if torch.cuda.is_available():
        device = resolve_device("cuda", require_cuda=True)
        assert device.type == "cuda"
    else:
        with pytest.raises(RuntimeError, match="CUDA"):
            resolve_device("cuda", require_cuda=True)


def test_train_flat_cli_overfit(isolated_workspace: Path) -> None:
    result = runner.invoke(
        app,
        [
            "train",
            "flat",
            "--overfit-fixture",
            "--device",
            "cpu",
            "--max-epochs",
            "250",
            "--run-id",
            "cli-overfit",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "cli-overfit" in result.output
    assert (isolated_workspace / "artifacts" / "runs" / "cli-overfit" / "metrics.json").is_file()


def test_record_to_batch_shapes() -> None:
    bundle = make_synthetic_overfit_bundle(n_samples=2, n_cpgs=6, n_genes=3, n_classes=2)
    batch = record_to_batch(bundle["records"][0], n_genes=bundle["n_genes"])
    assert batch.cpg_features.shape[0] == 6
    assert batch.cpg_to_gene.shape == (6,)
    assert batch.tissue_target.shape == (1,)
