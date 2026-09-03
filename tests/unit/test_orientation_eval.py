"""Orientation contract v2 and flat_region feature path tests."""

from __future__ import annotations

import numpy as np
import torch

from mbs.scoring.orientation import (
    accumulate_signed_gene_mean_m,
    orient_mbs_array,
    signed_gene_mean_m,
)
from mbs.training.feature_schema import FLAT_REGION, m_column_index, observed_column_index
from mbs.training.flat_region_features import (
    apply_flat_region_feature_mode,
    flat_region_input_dim,
    gather_flat_region_features,
)
from mbs.training.multitask import MultitaskHeads


def test_flat_region_m_column_is_zero() -> None:
    assert m_column_index(feature_schema=FLAT_REGION, include_m_value=True) == 0
    assert observed_column_index(feature_schema=FLAT_REGION) == -1


def test_observed_mask_excludes_unobserved_from_signed_m() -> None:
    m = np.array([1.0, 99.0, 3.0], dtype=np.float64)
    genes = np.array([0, 0, 1], dtype=np.int64)
    obs = np.array([True, False, True], dtype=bool)
    out = signed_gene_mean_m(m, genes, n_genes=2, observed_mask=obs)
    np.testing.assert_allclose(out, [1.0, 3.0])


def test_accumulate_signed_gene_mean_m_observed_batches() -> None:
    out = accumulate_signed_gene_mean_m(
        n_genes=2,
        cpg_m_batches=[np.array([1.0, 5.0]), np.array([3.0])],
        cpg_to_gene_batches=[np.array([0, 0]), np.array([1])],
        observed_batches=[np.array([True, False]), np.array([True])],
    )
    np.testing.assert_allclose(out, [1.0, 3.0])


def test_orient_mbs_array_flips() -> None:
    mbs = np.array([[0.2, 0.8]], dtype=np.float32)
    flipped = orient_mbs_array(mbs, score_polarity="flipped")
    np.testing.assert_allclose(flipped, [[0.8, 0.2]])


def test_multitask_head_logit_invariance_under_orientation() -> None:
    """W·(x-0.5) equals (-W)·((1-x)-0.5) with unchanged bias (ADR 0008 pairing)."""
    n_genes = 4
    heads = MultitaskHeads(n_genes, n_tissue_classes=3, sex_enabled=False)
    with torch.no_grad():
        heads.tissue_head.gene_weight.normal_(0.0, 0.1)
        heads.tissue_head.bias.zero_()
    mbs = torch.rand(2, n_genes)
    present = torch.ones(2, n_genes, dtype=torch.bool)
    logits_raw = heads.forward_tissue(mbs, present)
    with torch.no_grad():
        w = heads.tissue_head.gene_weight.clone()
        heads.tissue_head.gene_weight.copy_(-w)
    logits_neg_on_flip = heads.forward_tissue(1.0 - mbs, present)
    assert torch.allclose(logits_raw, logits_neg_on_flip, atol=1e-5)


def test_gather_flat_region_filters_unobserved_edges() -> None:
    from mbs.training.cascade_assign import build_cascade_assignment
    from mbs.training.cascade_loop import make_synthetic_cascade_tables
    from mbs.training.flat_region_features import build_flat_region_gene_index

    tables = make_synthetic_cascade_tables(seed=7)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    beta = np.asarray(tables["betas"][0], dtype=np.float32)
    beta[0] = np.nan
    feats, genes = gather_flat_region_features(beta_row=beta, index=index)
    assert feats.shape[0] == genes.shape[0]
    assert feats.shape[0] < index.n_edges
    assert np.all(feats[:, -1] > 0.5)


def test_flat_region_feature_mode_ablation_zeros_role_block() -> None:
    dim = flat_region_input_dim()
    feats = np.zeros((2, dim), dtype=np.float32)
    feats[:, 1] = 1.0
    out = apply_flat_region_feature_mode(feats, "m_only")
    assert float(out[:, 1].sum()) == 0.0


# ── evaluate_flat_mbs_e2e orientation contract ────────────────────────────────


def _make_eval_inputs(
    n_samples: int = 6, n_genes: int = 4, n_tissue: int = 3
) -> tuple:
    """Synthetic inputs for evaluate_flat_mbs_e2e tests."""
    from mbs.training.multitask import MultitaskHeads
    from mbs.training.phenotypes import SamplePhenotype

    rng = np.random.default_rng(0)
    mbs = rng.random((n_samples, n_genes)).astype(np.float32)
    present = rng.integers(0, 2, (n_samples, n_genes)).astype(bool)

    def _pheno(i: int) -> SamplePhenotype:
        return SamplePhenotype(
            sample_id=f"s{i}",
            cell_type="tissue",
            donor_id=None,
            title=f"s{i}",
            class_index=i % n_tissue,
            study_id="st0",
            age=float(30 + i * 5),
            age_mask=True,
            tissue_mask=True,
            sex_mask=True,
            sex_class_index=i % 2,
        )

    phenotypes = [_pheno(i) for i in range(n_samples)]
    heads = MultitaskHeads(n_genes, n_tissue_classes=n_tissue, sex_enabled=True)
    return mbs, present, phenotypes, heads


def test_v2_contract_heads_see_raw_mbs_when_flipped() -> None:
    """Contract v2: phenotype heads receive raw MBS even when score_polarity='flipped'."""
    from mbs.training.flat_stage_a_eval import evaluate_flat_mbs_e2e
    import torch

    mbs, present, phenotypes, heads = _make_eval_inputs()

    class _CapturingHeads:
        """Thin wrapper that records what score tensor heads receive."""
        def __init__(self, inner: MultitaskHeads) -> None:
            self._inner = inner
            self.received: list[np.ndarray] = []

        def eval(self) -> None:
            self._inner.eval()

        def forward_age(self, mbs_t: torch.Tensor, present_t: torch.Tensor) -> torch.Tensor:
            self.received.append(mbs_t.detach().cpu().numpy())
            return self._inner.forward_age(mbs_t, present_t)

        def forward_tissue(self, mbs_t: torch.Tensor, present_t: torch.Tensor) -> torch.Tensor:
            return self._inner.forward_tissue(mbs_t, present_t)

        def forward_sex(self, mbs_t: torch.Tensor, present_t: torch.Tensor):
            return self._inner.forward_sex(mbs_t, present_t)

        @property
        def sex_head(self):
            return self._inner.sex_head

    capturing = _CapturingHeads(heads)
    evaluate_flat_mbs_e2e(
        mbs_test=mbs,
        present_test=present,
        phenotypes_train=phenotypes,
        phenotypes_test=phenotypes,
        heads=capturing,  # type: ignore[arg-type]
        class_names=["t0", "t1", "t2"],
        device=torch.device("cpu"),
        age_mean=0.0,
        age_std=1.0,
        score_polarity="flipped",
        legacy_negated_heads=False,
    )
    assert capturing.received, "forward_age was never called"
    received = capturing.received[0]
    # Heads must receive raw MBS, not 1-MBS
    np.testing.assert_allclose(
        received, mbs, atol=1e-6,
        err_msg="v2 contract: heads should receive raw MBS, not 1-MBS"
    )


def test_legacy_path_heads_see_one_minus_mbs_when_flipped() -> None:
    """Legacy repair: heads receive 1-MBS when legacy_negated_heads=True and flipped."""
    from mbs.training.flat_stage_a_eval import evaluate_flat_mbs_e2e
    import torch

    mbs, present, phenotypes, heads = _make_eval_inputs()

    received_scores: list[np.ndarray] = []
    _orig_age = heads.forward_age

    def _capture_age(mbs_t: torch.Tensor, present_t: torch.Tensor) -> torch.Tensor:
        received_scores.append(mbs_t.detach().cpu().numpy())
        return _orig_age(mbs_t, present_t)

    heads.forward_age = _capture_age  # type: ignore[method-assign]
    evaluate_flat_mbs_e2e(
        mbs_test=mbs,
        present_test=present,
        phenotypes_train=phenotypes,
        phenotypes_test=phenotypes,
        heads=heads,
        class_names=["t0", "t1", "t2"],
        device=torch.device("cpu"),
        age_mean=0.0,
        age_std=1.0,
        score_polarity="flipped",
        legacy_negated_heads=True,
    )
    assert received_scores
    np.testing.assert_allclose(
        received_scores[0], 1.0 - mbs, atol=1e-6,
        err_msg="legacy path: heads should receive 1-MBS when legacy_negated_heads=True"
    )


def test_hyper_aligned_heads_always_see_raw() -> None:
    """hyper_aligned: heads receive raw MBS regardless of legacy flag."""
    from mbs.training.flat_stage_a_eval import evaluate_flat_mbs_e2e
    import torch

    mbs, present, phenotypes, heads = _make_eval_inputs()
    received_scores: list[np.ndarray] = []
    _orig_age = heads.forward_age

    def _capture(mbs_t: torch.Tensor, present_t: torch.Tensor) -> torch.Tensor:
        received_scores.append(mbs_t.detach().cpu().numpy())
        return _orig_age(mbs_t, present_t)

    heads.forward_age = _capture  # type: ignore[method-assign]
    evaluate_flat_mbs_e2e(
        mbs_test=mbs,
        present_test=present,
        phenotypes_train=phenotypes,
        phenotypes_test=phenotypes,
        heads=heads,
        class_names=["t0", "t1", "t2"],
        device=torch.device("cpu"),
        age_mean=0.0,
        age_std=1.0,
        score_polarity="hyper_aligned",
        legacy_negated_heads=False,
    )
    assert received_scores
    np.testing.assert_allclose(received_scores[0], mbs, atol=1e-6)
