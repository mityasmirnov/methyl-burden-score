"""Unit tests for the 7G′ weekend supervisor helpers."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def _load_sup():
    import sys

    path = ROOT / "scripts" / "supervise_7g_prime_weekend.py"
    spec = importlib.util.spec_from_file_location("supervise_7g_prime_weekend", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def sup():
    return _load_sup()


def _write_metrics(path: Path, *, best_epoch: int, max_epochs: int | None, tissue: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = {
        "best_epoch": best_epoch,
        "checkpoint_selection": {"best_epoch": best_epoch, "max_epochs": max_epochs},
        "evaluations": {
            "mbs_e2e": {
                "eval_split": "test",
                "metrics": {
                    "tissue": {"macro_f1": tissue},
                    "age": {"mae": 18.0},
                    "sex": {"auroc": 0.7},
                },
            }
        },
    }
    (path.parent / "scores").mkdir(parents=True, exist_ok=True)
    (path.parent / "scores" / "mbs.npy").write_bytes(b"x")
    path.write_text(json.dumps(blob), encoding="utf-8")


def test_validate_fold_ok(tmp_path: Path, monkeypatch, sup) -> None:
    monkeypatch.setattr(sup, "ART", tmp_path)
    root = tmp_path / "runA" / "fold_0"
    _write_metrics(root / "metrics.json", best_epoch=12, max_epochs=16, tissue=0.35)
    st = sup.validate_fold(
        arm_id="A",
        fold=0,
        kind="cascade",
        run_id="runA",
        min_ceiling=15,
    )
    assert st.status == "ok"
    assert st.matched16 is True


def test_validate_fold_rejects_5ep(tmp_path: Path, monkeypatch, sup) -> None:
    monkeypatch.setattr(sup, "ART", tmp_path)
    root = tmp_path / "runB" / "fold_0"
    _write_metrics(root / "metrics.json", best_epoch=4, max_epochs=5, tissue=0.35)
    st = sup.validate_fold(
        arm_id="B",
        fold=0,
        kind="cascade",
        run_id="runB",
        min_ceiling=15,
    )
    assert st.status == "incomplete"
    assert "not matched-16" in (st.reason or "")


def test_validate_fold_rejects_wrong_split(tmp_path: Path, monkeypatch, sup) -> None:
    monkeypatch.setattr(sup, "ART", tmp_path)
    root = tmp_path / "runC" / "fold_0"
    _write_metrics(root / "metrics.json", best_epoch=16, max_epochs=16, tissue=0.4)
    blob = json.loads((root / "metrics.json").read_text())
    blob["evaluations"]["mbs_e2e"]["eval_split"] = "val"
    (root / "metrics.json").write_text(json.dumps(blob), encoding="utf-8")
    st = sup.validate_fold(
        arm_id="C",
        fold=0,
        kind="cascade",
        run_id="runC",
        min_ceiling=15,
    )
    assert st.status == "incomplete"
    assert "eval_split" in (st.reason or "")
