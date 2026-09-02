"""Comparable panel × eval-mode ranking rows."""

from __future__ import annotations

from mbs.inspection.comparable_metrics import COMPARISON_SPECS, load_comparable_rows


def test_comparison_specs_cover_gene_and_full_panel() -> None:
    panels = {s["panel"] for s in COMPARISON_SPECS}
    modes = {s["eval_mode"] for s in COMPARISON_SPECS}
    assert "gene-linked" in panels
    assert "65k prefix" in panels
    assert "mbs_e2e" in modes
    assert "fusion_full" in modes


def test_load_comparable_rows_shape() -> None:
    rows = load_comparable_rows(
        __import__("pathlib").Path("/nonexistent"),
        classical_baselines_path=__import__("pathlib").Path("/nonexistent"),
    )
    assert len(rows) == len(COMPARISON_SPECS)
    assert all("arm" in r and "eval_mode" in r for r in rows)
