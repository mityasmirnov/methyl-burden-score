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
    assert "no architecture lock" in text.lower()
    assert "trait/seed-gene" in text.lower() or "Next real gate" in text


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
    assert lock["architecture_locked"] is False
    assert lock["mbs_e2e_valid"] is False
    assert lock["best_landed_cascade_arm"] is None
    assert "seed" in str(lock["lock_blocked_reason"]).lower() or "lock" in str(
        lock["lock_blocked_reason"]
    ).lower()


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
    # P2-G is a reference, not a retained architecture lock.
    assert lock["locked_cascade_arm"] is None
    assert lock["architecture_locked"] is False
    assert lock["best_landed_cascade_arm"] == "P2-G"
    assert lock["mbs_e2e_valid"] is True
    assert lock.get("next_gate") == "trait_seed_gene_stage_a_repeat"


def test_cascade_mode_row_metric_specific_fold_count() -> None:
    mod = _report_module()
    folds = [
        {
            "evaluations": {
                "mbs_enet": {"metrics": {"tissue": {"macro_f1": 0.28}, "age": {"mae": 20.0}}},
                "mbs_e2e": {
                    "eval_split": "test",
                    "metrics": {"tissue": {"macro_f1": 0.30}, "age": {"mae": 18.0}},
                },
            }
        },
        {
            "evaluations": {
                "mbs_e2e": {
                    "eval_split": "test",
                    "metrics": {"tissue": {"macro_f1": 0.31}, "age": {"mae": 19.0}},
                },
            }
        },
        {
            "evaluations": {
                "mbs_e2e": {
                    "eval_split": "test",
                    "metrics": {"tissue": {"macro_f1": 0.32}, "age": {"mae": 17.0}},
                },
            }
        },
    ]
    enet_row = mod._cascade_mode_row(folds, "N-light-gene-mean", "mbs_enet")
    assert enet_row is not None
    assert enet_row["n_folds"] == 1
    assert enet_row["n_folds_total"] == 3
    e2e_row = mod._cascade_mode_row(folds, "N-light-gene-mean", "mbs_e2e")
    assert e2e_row is not None
    assert e2e_row["n_folds"] == 3


def test_fold_epoch_helpers_and_training_epochs_section(tmp_path: Path) -> None:
    mod = _report_module()
    fold = {
        "best_epoch": 16,
        "checkpoint_selection": {
            "best_epoch": 16,
            "selection": "validation_tissue_macro_f1_then_age_mae",
            "epochs_trained": 21,
        },
        "history": [{"epoch": i} for i in range(1, 22)],
        "evaluations": {
            "mbs_e2e": {
                "eval_split": "test",
                "metrics": {"tissue": {"macro_f1": 0.30}},
            }
        },
        "score_manifest": {"orientation_contract_version": "2"},
    }
    assert mod._fold_best_epoch(fold) == 16
    assert mod._fold_epochs_trained(fold) == 21

    report_dir = tmp_path / "report"
    per_arm = report_dir / "per_arm"
    per_arm.mkdir(parents=True)
    (per_arm / "N-light-gene-max.json").write_text(
        __import__("json").dumps({"arm_id": "N-light-gene-max", "folds": [fold]}),
        encoding="utf-8",
    )
    lock = {"locked_cascade_arm": None, "lock_blocked_reason": "test"}
    mod.write_analysis(report_dir, lock=lock)
    text = (report_dir / "analysis.md").read_text(encoding="utf-8")
    assert "## Training epochs (ceiling / ran / best)" in text
    assert "N-light-gene-max" in text
    assert "best ep" in text
    assert "| 16 |" in text or "16" in text
    assert "21" in text
