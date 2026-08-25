"""Milestone 7E′ Hub multitask + hygiene unit tests (synthetic only)."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest
import yaml

from mbs.matrix.store import (
    create_betas_zarr,
    write_locus_index,
    write_sample_index,
)
from mbs.matrix.virtual_hub_store import (
    PACK_PRIORITY,
    VIRTUAL_MATRIX_ID,
    build_virtual_hub_store,
    open_betas_for_matrix,
)
from mbs.training.controls import evaluate_metadata_only_ceiling, fit_metadata_only
from mbs.training.encoder_config import resolve_encoder
from mbs.training.phenotype_table import (
    _build_hub_union_rows,
    build_hub_union_phenotype_table,
    build_tissue_ontology,
    default_sex_ontology,
)
from mbs.training.phenotypes import SamplePhenotype


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    root = scratch_base / f"7ep-{uuid4().hex}"
    root.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(root / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(root / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(root / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(root / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(root / "docker"))
    return root


def _write_pack(
    data_root: Path,
    *,
    family: str,
    sample_ids: list[str],
    locus_ids: list[int],
    fill: float,
) -> Path:
    mid = f"matrix-hub-{family}-full-v1"
    root = data_root / "canonical" / "matrices" / mid
    root.mkdir(parents=True, exist_ok=True)
    write_sample_index(
        root / "sample_index.parquet",
        sample_ids=sample_ids,
        source_sample_ids=sample_ids,
    )
    n_loci = len(locus_ids)
    write_locus_index(
        root / "locus_index.parquet",
        locus_ids=np.asarray(locus_ids, dtype=np.uint64),
        canonical_keys=np.asarray([f"k{i}" for i in locus_ids], dtype=object),
        probe_ids=np.asarray([f"cg{i}" for i in locus_ids], dtype=object),
    )
    arr = create_betas_zarr(root / "betas.zarr", n_samples=len(sample_ids), n_loci=n_loci)
    for i in range(len(sample_ids)):
        arr[i, :] = np.full(n_loci, fill + 0.01 * i, dtype=np.float32)
    # Minimal manifest (skip schema validate — discover only needs betas + sample_index).
    (root / "matrix_manifest.json").write_text(
        json.dumps(
            {
                "matrix_id": mid,
                "shape": [len(sample_ids), n_loci],
                "platform_id": "HM450",
            }
        ),
        encoding="utf-8",
    )
    pheno = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "study_id": [f"GSE_{family}_{i}" for i in range(len(sample_ids))],
            "platform": ["450K"] * len(sample_ids),
            "phenotype_family": [family] * len(sample_ids),
        }
    )
    if family == "age":
        pheno["age"] = [30.0 + i for i in range(len(sample_ids))]
        pheno["phenotype_value_numeric"] = pheno["age"]
    elif family == "tissue":
        pheno["tissue"] = ["blood" if i % 2 == 0 else "brain" for i in range(len(sample_ids))]
    elif family == "sex":
        pheno["sex"] = ["Female" if i % 2 == 0 else "Male" for i in range(len(sample_ids))]
    elif family == "disease":
        pheno["sample_type"] = ["case"] * len(sample_ids)
        pheno["phenotype_value"] = ["asthma"] * len(sample_ids)
    pheno.to_parquet(root / "sample_phenotypes.parquet", index=False)
    return root


def test_virtual_routing_overlap_and_locus_remap(workspace: Path) -> None:
    data_root = workspace / "data"
    # Age has loci 1,2,3; tissue has 2,3,4 → intersection 2,3 ordered as age.
    _write_pack(
        data_root,
        family="age",
        sample_ids=["GSM_A", "GSM_BOTH"],
        locus_ids=[1, 2, 3],
        fill=0.1,
    )
    _write_pack(
        data_root,
        family="tissue",
        sample_ids=["GSM_T", "GSM_BOTH"],
        locus_ids=[2, 3, 4],
        fill=0.7,
    )
    result = build_virtual_hub_store(data_root=data_root)
    assert result.matrix_id == VIRTUAL_MATRIX_ID
    assert result.n_samples == 3
    assert result.n_loci == 2  # intersection
    route = pd.read_parquet(result.route_path)
    by_sid = {str(r["sample_id"]): r for r in route.to_dict(orient="records")}
    assert by_sid["GSM_BOTH"]["family"] == "age"  # pack priority
    assert PACK_PRIORITY[0] == "age"
    assert by_sid["GSM_T"]["family"] == "tissue"
    betas = open_betas_for_matrix(result.output_dir)
    assert betas.shape == (3, 2)
    # GSM_BOTH from age pack fill 0.1 + 0.01*1 (second sample)
    both_row = int(by_sid["GSM_BOTH"]["row_index"])
    row = np.asarray(betas[both_row], dtype=np.float32)
    assert np.allclose(row, 0.11, atol=1e-5)
    # Missing locus → NaN would appear if intersection failed; both cols present.
    assert not np.isnan(row).any()
    # np.ix_ path used by Level-1
    block = betas[np.ix_([both_row], [0, 1])]
    assert block.shape == (1, 2)


def test_hub_union_unknown_ne_control(workspace: Path) -> None:
    data_root = workspace / "data"
    matrices = data_root / "canonical" / "matrices"
    # Age-labeled sample also in disease pack; disease-only sample unlabeled for age.
    _write_pack(data_root, family="age", sample_ids=["GSM_AGE"], locus_ids=[1, 2], fill=0.2)
    _write_pack(
        data_root,
        family="tissue",
        sample_ids=["GSM_AGE"],
        locus_ids=[1, 2],
        fill=0.25,
    )
    _write_pack(
        data_root,
        family="disease",
        sample_ids=["GSM_AGE", "GSM_DIS"],
        locus_ids=[1, 2],
        fill=0.3,
    )
    sample_index = pd.DataFrame(
        {
            "row_index": [0, 1],
            "sample_id": ["GSM_AGE", "GSM_DIS"],
            "source_sample_id": ["GSM_AGE", "GSM_DIS"],
        }
    )
    assert (matrices / "matrix-hub-disease-full-v1" / "sample_index.parquet").is_file()
    result = build_hub_union_phenotype_table(
        data_root=data_root,
        sample_index=sample_index,
        matrix_id="matrix-hub-nine-pack-virtual-v1",
        phenotype_table_path=data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_hub_nine_pack_v1.parquet",
        tissue_ontology_path=data_root
        / "canonical"
        / "phenotypes"
        / "tissue_ontology_hub_nine_pack_v1.yaml",
        sex_ontology_path=data_root
        / "canonical"
        / "phenotypes"
        / "sex_ontology_hub_nine_pack_v1.yaml",
        min_tissue_n=1,
    )
    frame = pd.read_parquet(result.phenotype_table_path)
    by = {str(r["sample_id"]): r for r in frame.to_dict(orient="records")}
    assert by["GSM_AGE"]["age_mask"] is True or by["GSM_AGE"]["age_mask"] == True  # noqa: E712
    assert by["GSM_DIS"]["age_mask"] is False or by["GSM_DIS"]["age_mask"] == False  # noqa: E712
    assert by["GSM_DIS"]["disease_mask"] is True or by["GSM_DIS"]["disease_mask"] == True  # noqa: E712
    # Missing cancer pack membership → cancer_mask False (unknown ≠ control).
    assert by["GSM_AGE"]["cancer_mask"] is False or by["GSM_AGE"]["cancer_mask"] == False  # noqa: E712


def test_build_hub_union_rows_disease_only_no_age() -> None:
    ontology = build_tissue_ontology(["blood", "brain"], min_n=1)
    sex_ontology = default_sex_ontology()
    sample_index = pd.DataFrame(
        {
            "row_index": [0],
            "sample_id": ["GSM_X"],
            "source_sample_id": ["GSM_X"],
        }
    )
    rows = _build_hub_union_rows(
        sample_index=sample_index,
        matrix_id="virt",
        age_pheno=pd.DataFrame(columns=["sample_id"]),
        tissue_pheno=pd.DataFrame(columns=["sample_id"]),
        sex_pheno=pd.DataFrame(columns=["sample_id"]),
        ontology=ontology,
        sex_ontology=sex_ontology,
        study_by={"GSM_X": "GSE1"},
        platform_by={"GSM_X": "HM450"},
        disease_members={"GSM_X"},
        cancer_members=set(),
    )
    assert len(rows) == 1
    assert rows[0]["age_mask"] is False
    assert rows[0]["disease_mask"] is True
    assert rows[0]["cancer_mask"] is False


def test_metadata_only_holdout_not_insample() -> None:
    # Distinct studies so one-hot cannot perfectly transfer → holdout ≠ train score.
    train = [
        SamplePhenotype(
            sample_id=f"tr{i}",
            cell_type="blood",
            donor_id=None,
            title="",
            class_index=0,
            study_id="GSE_TRAIN",
            age=40.0 + i,
            platform="HM450",
            age_mask=True,
            tissue_mask=True,
            sex_mask=True,
            sex_class_index=i % 2,
        )
        for i in range(8)
    ]
    holdout = [
        SamplePhenotype(
            sample_id=f"ho{i}",
            cell_type="blood",
            donor_id=None,
            title="",
            class_index=0,
            study_id="GSE_HOLD",
            age=55.0 + i,
            platform="HM450",
            age_mask=True,
            tissue_mask=True,
            sex_mask=True,
            sex_class_index=i % 2,
        )
        for i in range(4)
    ]
    ceiling = evaluate_metadata_only_ceiling(
        train=train,
        eval_sets={"external_test": holdout},
    )
    assert ceiling["protocol"] == "fit_train_score_holdout"
    assert "age" in ceiling["external_test"]
    hold_mae = float(ceiling["external_test"]["age"]["mae"])
    # In-sample train score (legacy path) should be much tighter.
    in_sample = fit_metadata_only(
        study_ids=["GSE_TRAIN"] * 8,
        platforms=["HM450"] * 8,
        tissues=["blood"] * 8,
        y=np.array([40.0 + i for i in range(8)]),
        task="regression",
    )
    assert hold_mae > float(in_sample["mae"])


def test_hub_yamls_share_matched_encoder() -> None:
    root = _repo_root()
    flat = yaml.safe_load(
        (root / "configs/experiment/stage0_flat_hub_multitask_v1.yaml").read_text(encoding="utf-8")
    )
    hier = yaml.safe_load(
        (root / "configs/experiment/stage0_hier_hub_multitask_v1.yaml").read_text(encoding="utf-8")
    )
    enc_f = resolve_encoder(flat["model"])
    enc_h = resolve_encoder(hier["model"])
    assert enc_f == enc_h
    assert enc_f["activation"] == "gelu"
    assert enc_f["dropout"] == 0.1
    assert enc_f["layer_norm"] is True
    assert enc_f["cpg_hidden_dim"] == 64
    assert flat["heads"]["disease"]["enabled"] is True
    assert flat["heads"]["cancer"]["enabled"] is True
    assert float(flat["loss"]["lambda_blood"]) == 0.0
    assert flat["heads"]["blood"]["enabled"] is False
    assert flat["pilot"].get("empty_as_control") is False
    assert hier["pilot"]["split_id"] == "hub-nine-pack-full-auto-v1"
    assert hier["pilot"].get("reuse_flat_split") is True


def test_pack_priority_order() -> None:
    assert PACK_PRIORITY[:5] == ("age", "tissue", "sex", "disease", "cancer")
