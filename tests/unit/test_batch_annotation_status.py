"""Unit tests for MethylationBatch annotation-status masks."""

from __future__ import annotations

import torch

from mbs.batch import (
    ANNOTATION_STATUS_AMBIGUOUS,
    ANNOTATION_STATUS_MAPPED,
    ANNOTATION_STATUS_MULTI_MAPPED,
    ANNOTATION_STATUS_UNMAPPED,
    MethylationBatch,
    annotation_status_masks,
)


def test_annotation_status_masks() -> None:
    status = torch.tensor(
        [
            ANNOTATION_STATUS_MAPPED,
            ANNOTATION_STATUS_UNMAPPED,
            ANNOTATION_STATUS_AMBIGUOUS,
            ANNOTATION_STATUS_MULTI_MAPPED,
        ],
        dtype=torch.long,
    )
    masks = annotation_status_masks(status)
    assert masks["mapped"].tolist() == [True, False, False, False]
    assert masks["unmapped"].tolist() == [False, True, False, False]
    assert masks["ambiguous"].tolist() == [False, False, True, False]
    assert masks["multi_mapped"].tolist() == [False, False, False, True]


def test_methylation_batch_validate_with_residual() -> None:
    batch = MethylationBatch(
        sample_ids=["s0"],
        cpg_features=torch.randn(2, 3),
        locus_row=torch.tensor([0, 1]),
        cpg_sample_index=torch.tensor([0, 0]),
        cpg_to_region=torch.tensor([0, 0]),
        region_type=torch.tensor([1]),
        region_to_gene=torch.tensor([0]),
        gene_sample_index=torch.tensor([0]),
        gene_panel_index=torch.tensor([0]),
        targets={"age": torch.tensor([40.0])},
        target_masks={"age": torch.tensor([True])},
        covariates={},
        annotation_status=torch.tensor(
            [ANNOTATION_STATUS_MAPPED, ANNOTATION_STATUS_MAPPED],
            dtype=torch.long,
        ),
        residual_features=torch.randn(1, 3),
        residual_sample_index=torch.tensor([0]),
    )
    batch.validate()
    masks = batch.annotation_masks()
    assert bool(masks["mapped"].all())
