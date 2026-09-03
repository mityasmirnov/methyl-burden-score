"""Unit tests for FlatDeepSetRegion Stage A evaluation helpers."""

from __future__ import annotations

import numpy as np
import torch

from mbs.models import FlatDeepSetRegion
from mbs.training.flat_region_features import flat_region_input_dim
from mbs.training.flat_stage_a_eval import (
    evaluate_flat_mbs_e2e,
    evaluate_flat_mbs_enet,
    evaluate_flat_mbs_linear_probe,
)
from mbs.training.multitask import MultitaskHeads
from mbs.training.phenotypes import SamplePhenotype


def _ph(i: int, *, tissue: int = 0) -> SamplePhenotype:
    return SamplePhenotype(
        sample_id=f"s{i}",
        cell_type="t",
        donor_id=None,
        title=f"s{i}",
        class_index=tissue,
        study_id=f"st{i // 3}",
        platform="p",
        age=40.0 + i,
        age_mask=True,
        tissue_mask=True,
        sex_mask=True,
        sex_class_index=i % 2,
    )


def test_flat_stage_a_eval_shapes() -> None:
    n_genes = 8
    n_train, n_test = 12, 6
    rng = np.random.default_rng(0)
    mbs = rng.normal(size=(n_train + n_test, n_genes)).astype(np.float32)
    present = np.ones_like(mbs, dtype=bool)
    train = [_ph(i, tissue=i % 3) for i in range(n_train)]
    test = [_ph(n_train + i, tissue=i % 3) for i in range(n_test)]
    arrays = {
        "age": np.asarray([float(p.age or 0.0) for p in train + test], dtype=np.float64),
        "age_mask": np.ones(n_train + n_test, dtype=bool),
        "tissue": np.asarray([p.class_index for p in train + test], dtype=np.int64),
        "tissue_mask": np.ones(n_train + n_test, dtype=bool),
        "sex": np.asarray([int(p.sex_class_index or 0) for p in train + test], dtype=np.int64),
        "sex_mask": np.ones(n_train + n_test, dtype=bool),
        "study_ids": np.asarray([str(p.study_id) for p in train + test], dtype=object),
    }
    train_idx = np.arange(n_train, dtype=np.int64)
    test_idx = np.arange(n_train, n_train + n_test, dtype=np.int64)
    class_names = ["a", "b", "c"]
    probe = evaluate_flat_mbs_linear_probe(
        mbs_all=mbs,
        train_idx=train_idx,
        test_idx=test_idx,
        arrays=arrays,
        class_names=class_names,
    )
    enet = evaluate_flat_mbs_enet(
        mbs_all=mbs,
        train_idx=train_idx,
        test_idx=test_idx,
        arrays=arrays,
        class_names=class_names,
    )
    assert probe["eval_split"] == "test"
    assert enet["eval_split"] == "test"
    assert "tissue" in probe["metrics"]
    assert "tissue" in enet["metrics"]

    model = FlatDeepSetRegion(flat_region_input_dim(), phi_hidden_dim=16, rho_hidden_dim=8)
    heads = MultitaskHeads(n_genes, 3, sex_enabled=True)
    # Untrained heads still return structured metrics on noise MBS.
    e2e = evaluate_flat_mbs_e2e(
        heads=heads,
        mbs_test=mbs[test_idx],
        present_test=present[test_idx],
        phenotypes_test=test,
        phenotypes_train=train,
        class_names=class_names,
        device=torch.device("cpu"),
        age_mean=40.0,
        age_std=10.0,
    )
    assert e2e["eval_split"] == "test"
    assert e2e["evaluation"] == "mbs_e2e"
    assert "macro_f1" in e2e["metrics"]["tissue"]
    del model
