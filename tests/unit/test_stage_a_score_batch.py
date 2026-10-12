"""Stage A scoring must keep edge-rows per batch at what training calibrated, not batch size."""

from __future__ import annotations

from mbs.training.flat_stage_a_eval import score_batch_size_for_graph


def test_k8_full_graph_scales_batch_down() -> None:
    # 2026-10-12 k8-rescore2 OOM: batch 111 calibrated on 147,438 edges, scored on 440,903.
    b = score_batch_size_for_graph(111, calibrated_edges=147_438, score_edges=440_903)
    assert b == 37
    assert b * 440_903 <= 111 * 147_438


def test_same_or_narrower_graph_keeps_batch() -> None:
    assert score_batch_size_for_graph(111, calibrated_edges=440_903, score_edges=440_903) == 111
    assert score_batch_size_for_graph(111, calibrated_edges=440_903, score_edges=100_000) == 111


def test_unknown_edges_and_floor() -> None:
    assert score_batch_size_for_graph(8, calibrated_edges=None, score_edges=1_000) == 8
    assert score_batch_size_for_graph(1, calibrated_edges=10, score_edges=1_000_000) == 1
