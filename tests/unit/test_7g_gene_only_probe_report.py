"""Report writer smoke for 7G′ Stage A."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _report_module():
    root = Path(__file__).resolve().parents[2]
    path = root / "scripts" / "write_7g_gene_only_probe_report.py"
    spec = importlib.util.spec_from_file_location("write_7g_gene_only_probe_report", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_write_analysis_empty(tmp_path: Path) -> None:
    mod = _report_module()
    report_dir = tmp_path / "report"
    (report_dir / "per_arm").mkdir(parents=True)
    lock = mod.build_lock_recommendation([], [])
    mod.write_analysis(report_dir, lock=lock)
    assert (report_dir / "analysis.md").is_file()
    assert (report_dir / "lock_recommendation.json").is_file()
    text = (report_dir / "analysis.md").read_text(encoding="utf-8")
    assert "mbs_e2e" in text
    assert "No architecture lock" in text


def test_lock_blocked_without_valid_mbs_e2e() -> None:
    mod = _report_module()
    contaminated_fold = {
        "evaluations": {
            "mbs_e2e": {
                "metrics": {"tissue": {"macro_f1": 0.69}},
            }
        }
    }
    valid_classical = {
        "folds": [
            {
                "arms": {
                    "C-mvalue-enet-G": {"tissue": {"macro_f1": 0.39}},
                }
            }
        ]
    }
    lock = mod.build_lock_recommendation(
        [{"arm_id": "P4-G", "mbs_e2e_f1": 0.69}],
        [{"arm_id": "C-mvalue-enet-G", "tissue_f1": 0.39}],
        cascade_folds_by_arm={"P4-G": [contaminated_fold]},
        classical_payload=valid_classical,
    )
    assert lock["locked_cascade_arm"] is None
    assert lock["mbs_e2e_valid"] is False
    assert "test-only" in str(lock["lock_blocked_reason"])


def test_lock_allowed_with_valid_mbs_e2e() -> None:
    mod = _report_module()
    valid_fold = {
        "evaluations": {
            "mbs_e2e": {
                "eval_split": "test",
                "metrics": {"tissue": {"macro_f1": 0.40}},
            }
        }
    }
    valid_classical = {
        "folds": [
            {
                "arms": {
                    "C-mvalue-enet-G": {"tissue": {"macro_f1": 0.35}},
                }
            }
        ]
    }
    lock = mod.build_lock_recommendation(
        [{"arm_id": "P2-G", "mbs_e2e_f1": 0.40}],
        [{"arm_id": "C-mvalue-enet-G", "tissue_f1": 0.35}],
        cascade_folds_by_arm={"P2-G": [valid_fold, valid_fold, valid_fold]},
        classical_payload=valid_classical,
    )
    assert lock["locked_cascade_arm"] == "P2-G"
    assert lock["mbs_e2e_valid"] is True
