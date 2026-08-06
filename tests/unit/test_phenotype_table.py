"""Unit tests for sample phenotype table + tissue ontology (Milestone 5c)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mbs.evaluation.splits import assert_no_study_leakage, build_study_grouped_split
from mbs.training.phenotype_table import (
    build_sample_phenotype_rows,
    build_tissue_ontology,
    normalize_platform,
    write_sample_phenotype_table,
    write_tissue_ontology,
)
from mbs.training.phenotypes import load_multitask_phenotypes


def test_normalize_platform_450k() -> None:
    assert normalize_platform("450K") == "HM450"
    assert normalize_platform("HM450") == "HM450"
    assert normalize_platform(None) is None


def test_build_table_gsm_dedupe_and_masks(tmp_path: Path) -> None:
    age_pheno = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE_A",
                "platform": "450K",
                "phenotype_value_numeric": 42.0,
                "tissue": "blood",
            },
            {
                "sample_id": "GSM2",
                "study_id": "GSE_A",
                "platform": "450K",
                "phenotype_value_numeric": 55.0,
                "tissue": "brain",
            },
        ]
    )
    tissue_pheno = pd.DataFrame(
        [
            {
                "sample_id": "GSM2",
                "study_id": "GSE_A",
                "platform": "450K",
                "phenotype_value": "brain",
            },
            {
                "sample_id": "GSM3",
                "study_id": "GSE_T",
                "platform": "450K",
                "phenotype_value": "lung",
            },
        ]
    )
    sample_index = pd.DataFrame(
        [
            {"row_index": 0, "sample_id": "GSM1", "source_sample_id": "GSM1"},
            {"row_index": 1, "sample_id": "GSM2", "source_sample_id": "GSM2"},
            {"row_index": 2, "sample_id": "GSM3", "source_sample_id": "GSM3"},
        ]
    )
    ontology = build_tissue_ontology(
        ["brain", "lung", "brain", "lung"],
        min_n=1,
    )
    table = build_sample_phenotype_rows(
        age_pheno=age_pheno,
        tissue_pheno=tissue_pheno,
        sample_index=sample_index,
        matrix_id="matrix-test",
        ontology=ontology,
    )
    assert len(table) == 3
    by_id = table.set_index("sample_id")
    assert bool(by_id.loc["GSM1", "age_mask"]) is True
    assert bool(by_id.loc["GSM1", "tissue_mask"]) is False  # blood filtered out
    assert by_id.loc["GSM1", "platform_id"] == "HM450"
    assert bool(by_id.loc["GSM2", "age_mask"]) is True
    assert bool(by_id.loc["GSM2", "tissue_mask"]) is True
    assert by_id.loc["GSM2", "phenotype_family"] == "multi"
    assert bool(by_id.loc["GSM3", "age_mask"]) is False
    assert bool(by_id.loc["GSM3", "tissue_mask"]) is True

    out = tmp_path / "sample_phenotype_table.parquet"
    write_sample_phenotype_table(out, table)
    ont_path = tmp_path / "tissue_ontology.yaml"
    write_tissue_ontology(ont_path, ontology)

    phenotypes, class_names = load_multitask_phenotypes(out, class_names=ontology.class_names)
    assert class_names == ontology.class_names
    assert len(phenotypes) == 3
    p2 = next(p for p in phenotypes if p.sample_id == "GSM2")
    assert p2.age_mask and p2.tissue_mask
    assert p2.age == 55.0

    sample_rows = [
        {"sample_id": p.sample_id, "study_id": p.study_id or p.donor_id} for p in phenotypes
    ]
    split = build_study_grouped_split(
        sample_rows,
        train_studies=["GSE_A"],
        validation_studies=["GSE_T"],
        external_test_studies=[],
        split_id="test-split",
    )
    assert_no_study_leakage(split)
