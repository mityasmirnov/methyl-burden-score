"""Milestone 7C architecture corrections (sampler, splits, graph-v2, arms)."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest
import torch

from mbs.annotation.build import attach_cgi_tile_systems
from mbs.evaluation.splits import partition_studies_constrained
from mbs.models import FlatDeepSet, HierarchicalDeepSet, SharedMLP, center_mask_scores
from mbs.scoring.orientation import (
    apply_orientation,
    orient_run_scores,
    polarity_from_correlation,
    signed_gene_mean_m,
)
from mbs.training.branch import hub_longform_ready, train_branch_arm
from mbs.training.controls import apply_feature_control, permute_labels_within_study
from mbs.training.dataset import FlatBatch
from mbs.training.direct_cpg import fit_direct_elasticnet
from mbs.training.encoder_config import resolve_encoder
from mbs.training.loop import train_flat_baseline
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss
from mbs.training.phenotypes import SamplePhenotype, load_longform_multilabel
from mbs.training.sampler import iter_epoch_batches


def test_sampler_deterministic_and_budget() -> None:
    n = 6
    tokens = [3, 3, 3, 3, 3, 10]
    studies = ["A", "A", "B", "B", "C", "C"]
    tasks = ["t"] * n
    a = list(
        iter_epoch_batches(
            n,
            n_tokens=tokens,
            study_ids=studies,
            task_keys=tasks,
            batch_token_budget=6,
            batch_size=8,
            seed=7,
            epoch=1,
        )
    )
    b = list(
        iter_epoch_batches(
            n,
            n_tokens=tokens,
            study_ids=studies,
            task_keys=tasks,
            batch_token_budget=6,
            batch_size=8,
            seed=7,
            epoch=1,
        )
    )
    assert a == b
    for batch in a:
        if len(batch) == 1:
            continue
        assert sum(tokens[i] for i in batch) <= 6


def test_center_mask_absent_gene_is_zero() -> None:
    mbs = torch.tensor([[0.9, 0.1]])
    present = torch.tensor([[True, False]])
    centered = center_mask_scores(mbs, present)
    assert float(centered[0, 1]) == 0.0
    heads = MultitaskHeads(2, 3, sex_enabled=True)
    assert heads.forward_age(mbs, present).shape[0] == 1
    assert heads.forward_sex(mbs, present).shape[0] == 1
    assert heads.forward_tissue(mbs, present).shape[0] == 1


def test_constrained_split_blocks_donor_leak() -> None:
    samples = [
        {"sample_id": "s1", "study_id": "G1", "donor_id": "D", "platform": "HM450"},
        {"sample_id": "s2", "study_id": "G2", "donor_id": "D", "platform": "HM450"},
        {"sample_id": "s3", "study_id": "G3", "donor_id": None, "platform": "EPIC"},
        {"sample_id": "s4", "study_id": "G4", "donor_id": None, "platform": "EPIC"},
        {"sample_id": "s5", "study_id": "G5", "donor_id": None, "platform": "EPIC"},
    ]
    split = partition_studies_constrained(samples, seed=1, split_id="c-test")
    role = {r["sample_id"]: r["role"] for r in split["samples"]}
    assert role["s1"] == role["s2"]
    assert "constraints" in split


def test_orientation_flips_negative_correlation() -> None:
    mbs = np.array([0.1, 0.2, 0.9], dtype=np.float64)
    signed = np.array([0.9, 0.8, 0.1], dtype=np.float64)
    assert polarity_from_correlation(mbs, signed) == "flipped"
    out = apply_orientation(mbs, signed_m=signed)
    assert out["score_polarity"] == "flipped"
    assert out["score_family"] == "predictive_mbs"
    np.testing.assert_allclose(out["mbs"], 1.0 - mbs)


def test_feature_controls_zero_channels() -> None:
    feats = np.ones((4, 5), dtype=np.float32)
    static = apply_feature_control(feats, mode="static_only")
    assert np.allclose(static[:, :2], 0.0)
    cov = apply_feature_control(feats, mode="coverage_only")
    assert np.allclose(cov[:, :-1], 0.0)
    assert np.allclose(cov[:, -1], 1.0)


def test_label_permutation_preserves_study_multiset() -> None:
    ph = [
        SamplePhenotype("a", "t", None, "t", 0, study_id="S1", age=10.0),
        SamplePhenotype("b", "u", None, "u", 1, study_id="S1", age=20.0),
        SamplePhenotype("c", "t", None, "t", 0, study_id="S2", age=30.0),
    ]
    out = permute_labels_within_study(ph, seed=3)
    s1 = sorted(float(x.age) for x in out if x.study_id == "S1" and x.age is not None)
    assert s1 == [10.0, 20.0]


def test_graph_v2_assignment_disjoint() -> None:
    loci = pd.DataFrame(
        {
            "locus_id": [1, 2, 3, 4, 5, 6],
            "chromosome": ["chr1"] * 6,
            "position": [10, 20, 30, 40, 50, 60],
            "cpg_context": ["island", "island", "open_sea", "open_sea", "open_sea", "open_sea"],
            "mapping_status": ["mapped"] * 5 + ["unmapped"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": ["g:body"],
            "gene_id": ["ENSG1"],
            "region_type": ["gene_body"],
            "chromosome": ["chr1"],
            "start": [1],
            "end": [15],
            "strand": ["+"],
            "source_version": ["gencode"],
        }
    )
    edges = pd.DataFrame(
        {
            "locus_id": [1],
            "region_id": ["g:body"],
            "edge_weight": [1.0],
            "evidence_type": ["gene"],
            "primary_gene_role": [True],
        }
    )
    out_r, out_e = attach_cgi_tile_systems(loci, regions, edges, tile_target_n_cpgs=2)
    sys_by_region = dict(zip(out_r["region_id"], out_r["region_system"], strict=True))
    by_sys: dict[str, set[int]] = {"gene": set(), "rbs": set(), "tbs": set()}
    for lid, rid in zip(out_e["locus_id"], out_e["region_id"], strict=True):
        by_sys[str(sys_by_region[rid])].add(int(lid))
    assert by_sys["gene"].isdisjoint(by_sys["rbs"])
    assert 6 not in by_sys["tbs"]
    assert 6 not in by_sys["rbs"]
    assert 2 in by_sys["rbs"]
    assert {3, 4, 5} <= by_sys["tbs"]


def test_direct_elasticnet_and_arms(monkeypatch: pytest.MonkeyPatch) -> None:
    repo = Path(__file__).resolve().parents[2]
    workspace = repo / "scratch" / "pytest" / f"7c-branch-{uuid4().hex}"
    workspace.mkdir(parents=True)
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    rng = np.random.default_rng(0)
    z = rng.normal(size=(24, 6))
    obs = np.ones_like(z, dtype=bool)
    y = z[:, 0]
    studies = np.array(["A"] * 12 + ["B"] * 12)
    fitted = fit_direct_elasticnet(z, obs, y, studies, min_studies=2, alpha=0.01)
    assert fitted["n_loci"] == 6
    cfg = {
        "experiment": {"seed": 1, "name": "7c"},
        "training": {"max_epochs": 1, "batch_size": 4},
        "logging": {"tensorboard": False},
    }
    data_root = workspace / "data"
    art_root = workspace / "artifacts"
    data_root.mkdir()
    art_root.mkdir()
    for arm in ("direct", "gene", "rbs", "tbs"):
        out = train_branch_arm(
            arm=arm,
            project_root=repo,
            data_root=data_root,
            artifact_root=art_root,
            config=cfg,
            run_id=f"7c-{arm}",
            device="cpu",
            overfit_fixture=True,
        )
        assert out["arm"] == arm


def test_matched_encoder_sharedmlp() -> None:
    cfg = {
        "encoder": {
            "activation": "gelu",
            "dropout": 0.1,
            "layer_norm": True,
            "cpg_hidden_dim": 64,
        }
    }
    enc = resolve_encoder(cfg)
    flat = FlatDeepSet(
        4,
        phi_hidden_dim=enc["cpg_hidden_dim"],
        dropout=enc["dropout"],
        activation=enc["activation"],
        layer_norm=enc["layer_norm"],
    )
    hier = HierarchicalDeepSet(
        4,
        5,
        cpg_hidden_dim=enc["cpg_hidden_dim"],
        dropout=enc["dropout"],
        activation=enc["activation"],
        layer_norm=enc["layer_norm"],
    )
    assert type(flat.phi) is SharedMLP
    assert type(hier.cpg_encoder) is SharedMLP
    assert enc["activation"] == "gelu"
    assert enc["layer_norm"] is True
    assert enc["cpg_hidden_dim"] == 64


def test_disease_masked_bce_unknown_not_control() -> None:
    heads = MultitaskHeads(3, 2, n_disease_labels=2)
    mbs = torch.rand(2, 3)
    present = torch.ones(2, 3, dtype=torch.bool)
    batch = FlatBatch(
        sample_ids=["a", "b"],
        cpg_features=torch.zeros(1, 4),
        cpg_to_gene=torch.zeros(1, dtype=torch.int64),
        n_genes=3,
        tissue_target=torch.zeros(2, dtype=torch.int64),
        tissue_mask=torch.zeros(2, dtype=torch.bool),
        age_target=None,
        age_mask=torch.zeros(2, dtype=torch.bool),
        disease_target=torch.tensor([[1.0, 0.0], [0.0, 0.0]]),
        disease_mask=torch.tensor([[True, True], [False, False]]),
    )
    result = masked_multitask_loss(mbs=mbs, present=present, heads=heads, batch=batch)
    assert result.metrics["disease_n"] == 2.0
    assert hub_longform_ready(Path("/data/no-such-root"), "matrix-hub-disease-full-v1") is False


def test_signed_gene_mean_m_and_orient_helper() -> None:
    signed = signed_gene_mean_m(
        np.array([1.0, 3.0, 5.0]),
        np.array([0, 0, 1]),
        n_genes=2,
    )
    np.testing.assert_allclose(signed, [2.0, 5.0])
    mbs = np.array([[0.1, 0.9], [0.2, 0.8]])
    out = orient_run_scores(mbs, signed_m=np.array([0.9, 0.1]))
    assert out["score_polarity"] == "flipped"


def test_longform_multilabel_unknown_not_control(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2"],
            "phenotype_value": ["T2D", "CAD", "T2D"],
        }
    )
    path = tmp_path / "sample_phenotypes.parquet"
    frame.to_parquet(path, index=False)
    maps = load_longform_multilabel(path, sample_ids=["s1", "s2", "s3"])
    assert maps.label_names == ("CAD", "T2D")
    np.testing.assert_array_equal(maps.targets["s1"], [1.0, 1.0])
    assert maps.masks["s1"].tolist() == [True, True]
    assert maps.masks["s3"].tolist() == [False, False]
    assert float(maps.targets["s3"].sum()) == 0.0


def test_overfit_writes_score_manifest(monkeypatch: pytest.MonkeyPatch) -> None:
    repo = Path(__file__).resolve().parents[2]
    workspace = repo / "scratch" / "pytest" / f"7c-orient-{uuid4().hex}"
    workspace.mkdir(parents=True)
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(workspace / "data"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(workspace / "artifacts"))
    (workspace / "data").mkdir()
    (workspace / "artifacts").mkdir()
    result = train_flat_baseline(
        project_root=repo,
        data_root=workspace / "data",
        artifact_root=workspace / "artifacts",
        config={
            "experiment": {"seed": 2, "name": "orient"},
            "training": {"max_epochs": 1, "batch_size": 4},
            "logging": {"tensorboard": False},
            "features": {"methylation": {"m_value": True}},
        },
        run_id="7c-orient-fixture",
        device_str="cpu",
        overfit_fixture=True,
    )
    manifest_path = result.run_dir / "score_manifest.json"
    assert manifest_path.is_file()
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["score_family"] == "predictive_mbs"
    assert payload["score_polarity"] in {"hyper_aligned", "flipped"}
    assert "anchor_recipe" in payload


def test_hub_disease_cancer_longform_sidecar_join() -> None:
    """Real 7B matrices: multi-hot join when Hub full packs are on disk."""
    data_root = Path(__file__).resolve().parents[2] / "data"
    disease_id = "matrix-hub-disease-full-v1"
    cancer_id = "matrix-hub-cancer-full-v1"
    if not hub_longform_ready(data_root, disease_id):
        pytest.skip("matrix-hub-disease-full-v1 not present")
    if not hub_longform_ready(data_root, cancer_id):
        pytest.skip("matrix-hub-cancer-full-v1 not present")

    for mid, min_labels in ((disease_id, 5), (cancer_id, 5)):
        idx = data_root / "canonical" / "matrices" / mid / "sample_index.parquet"
        side = data_root / "canonical" / "matrices" / mid / "sample_phenotypes.parquet"
        sample_ids = pd.read_parquet(idx, columns=["sample_id"])["sample_id"].astype(str).tolist()
        assert sample_ids
        maps = load_longform_multilabel(side, sample_ids=sample_ids[:64], min_count=5)
        assert len(maps.label_names) >= min_labels
        # At least one of the first 64 samples should have an observed label.
        assert any(bool(m.any()) for m in maps.masks.values())
        # Unknown sample stays all-masked.
        unknown = load_longform_multilabel(
            side, sample_ids=["__no_such_gsm__"], min_count=5, label_names=list(maps.label_names)
        )
        assert not unknown.masks["__no_such_gsm__"].any()
        assert float(unknown.targets["__no_such_gsm__"].sum()) == 0.0
