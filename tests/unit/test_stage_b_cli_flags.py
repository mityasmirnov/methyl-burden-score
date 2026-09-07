"""CLI helpers for Stage B runner flags."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _stage_b_module():
    root = Path(__file__).resolve().parents[2]
    path = root / "scripts" / "run_7g_prime_stage_b.py"
    spec = importlib.util.spec_from_file_location("run_7g_prime_stage_b", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parse_folds_arg_default_all() -> None:
    mod = _stage_b_module()
    assert mod.parse_folds_arg(None, 3) == [0, 1, 2]
    assert mod.parse_folds_arg("", 3) == [0, 1, 2]


def test_parse_folds_arg_csv() -> None:
    mod = _stage_b_module()
    assert mod.parse_folds_arg("0", 3) == [0]
    assert mod.parse_folds_arg("0,2", 3) == [0, 2]


def test_parse_folds_arg_out_of_range() -> None:
    mod = _stage_b_module()
    with pytest.raises(ValueError, match="out of range"):
        mod.parse_folds_arg("3", 3)
