#!/usr/bin/env python3
"""Apply 16-epoch promotion decision rules to gene-only probe per_arm metrics.

Writes ``reports/inspection/stage0_7g_gene_only_probe/promotion_decision.json``
and appends a short section to analysis.md when enough 16-ep folds exist.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports/inspection/stage0_7g_gene_only_probe"
P2_REF_TISSUE = 0.373
CLASSICAL_TISSUE = 0.388


def _e2e_tissue(folds: list[dict[str, Any]]) -> list[float]:
    out: list[float] = []
    for f in folds:
        ev = ((f.get("evaluations") or {}).get("mbs_e2e") or {}).get("metrics") or {}
        v = (ev.get("tissue") or {}).get("macro_f1")
        if v is not None:
            out.append(float(v))
    return out


def _e2e_age(folds: list[dict[str, Any]]) -> list[float]:
    out: list[float] = []
    for f in folds:
        ev = ((f.get("evaluations") or {}).get("mbs_e2e") or {}).get("metrics") or {}
        v = (ev.get("age") or {}).get("mae")
        if v is not None:
            out.append(float(v))
    return out


def _e2e_sex(folds: list[dict[str, Any]]) -> list[float]:
    out: list[float] = []
    for f in folds:
        ev = ((f.get("evaluations") or {}).get("mbs_e2e") or {}).get("metrics") or {}
        v = (ev.get("sex") or {}).get("auroc")
        if v is not None:
            out.append(float(v))
    return out


def _ceiling(folds: list[dict[str, Any]]) -> int | None:
    vals: list[int] = []
    for f in folds:
        ckpt = f.get("checkpoint_selection") or {}
        if isinstance(ckpt, dict) and ckpt.get("max_epochs") is not None:
            vals.append(int(ckpt["max_epochs"]))
        elif f.get("max_epochs") is not None:
            vals.append(int(f["max_epochs"]))
    return max(vals) if vals else None


def _load_arm(report_dir: Path, arm_id: str) -> dict[str, Any]:
    path = report_dir / "per_arm" / f"{arm_id}.json"
    if not path.is_file():
        return {"arm_id": arm_id, "folds": [], "n_folds": 0}
    return json.loads(path.read_text(encoding="utf-8"))


def decide(report_dir: Path) -> dict[str, Any]:
    arms = {
        "P2-G": _load_arm(report_dir, "P2-G"),
        "N-light-gene-max": _load_arm(report_dir, "N-light-gene-max"),
        "N-light-gene-mean": _load_arm(report_dir, "N-light-gene-mean"),
        "N-cascade-scalar-mean-max": _load_arm(report_dir, "N-cascade-scalar-mean-max"),
        "N-cascade-scalar-max-mean": _load_arm(report_dir, "N-cascade-scalar-max-mean"),
        "N-cascade-vector-mean-max": _load_arm(report_dir, "N-cascade-vector-mean-max"),
    }
    summary: dict[str, Any] = {"arms": {}, "rules": {}, "next_gate": None}
    for arm_id, payload in arms.items():
        folds = payload.get("folds") or []
        tissue = _e2e_tissue(folds)
        age = _e2e_age(folds)
        sex = _e2e_sex(folds)
        summary["arms"][arm_id] = {
            "n_folds": len(folds),
            "ceiling": _ceiling(folds),
            "tissue_f1_mean": float(np.mean(tissue)) if tissue else None,
            "age_mae_mean": float(np.mean(age)) if age else None,
            "sex_auroc_mean": float(np.mean(sex)) if sex else None,
            "matched_16ep": (_ceiling(folds) or 0) >= 15 and len(tissue) >= 3,
        }

    light = summary["arms"]["N-light-gene-max"]
    vec = summary["arms"]["N-cascade-vector-mean-max"]
    mean_max = summary["arms"]["N-cascade-scalar-mean-max"]
    max_mean = summary["arms"]["N-cascade-scalar-max-mean"]
    p2 = summary["arms"]["P2-G"]
    p2_f1 = p2.get("tissue_f1_mean") or P2_REF_TISSUE

    # Rule 1: one-hop max ≈ P2-G
    r1 = False
    if light.get("matched_16ep") and light.get("tissue_f1_mean") is not None:
        r1 = abs(float(light["tissue_f1_mean"]) - float(p2_f1)) <= 0.03
    summary["rules"]["one_hop_max_near_p2"] = {
        "fired": r1,
        "detail": "Prefer smaller DeepRVAT-like one-hop max if within 0.03 tissue F1 of P2-G.",
        "light_f1": light.get("tissue_f1_mean"),
        "p2_f1": p2_f1,
    }

    # Rule 2: vector improves age/sex e2e but loses after gene pool → typed-RBS
    r2 = False
    if vec.get("matched_16ep") and vec.get("age_mae_mean") is not None:
        # Compare to 5-ep baseline age ~22.75; improvement if age MAE drops toward P2 (~15.6)
        r2 = float(vec["age_mae_mean"]) < 18.0 and float(vec.get("tissue_f1_mean") or 0) < p2_f1 - 0.01
    summary["rules"]["vector_age_improves_but_gene_pool_loses"] = {
        "fired": r2,
        "detail": "Proceed with typed-RBS aggregation, not scalar MBS.",
        "vector_age_mae": vec.get("age_mae_mean"),
        "vector_tissue_f1": vec.get("tissue_f1_mean"),
    }

    # Rule 3: scalar mixed closes gap to P2 → retain full 2×2
    r3 = False
    mixed_f1s = [
        x
        for x in (mean_max.get("tissue_f1_mean"), max_mean.get("tissue_f1_mean"))
        if x is not None and (mean_max.get("matched_16ep") or max_mean.get("matched_16ep"))
    ]
    if mixed_f1s:
        r3 = any(abs(float(f) - float(p2_f1)) <= 0.02 for f in mixed_f1s)
    summary["rules"]["scalar_mixed_closes_gap"] = {
        "fired": r3,
        "detail": "Retain full 2×2 pooling; max/max is not locked.",
        "mean_max_f1": mean_max.get("tissue_f1_mean"),
        "max_mean_f1": max_mean.get("tissue_f1_mean"),
    }

    # Rule 4: nothing beats P2 or classical → stop sweeps, start seed-mask
    promoted = [
        a
        for a in (
            light.get("tissue_f1_mean"),
            mean_max.get("tissue_f1_mean"),
            max_mean.get("tissue_f1_mean"),
            vec.get("tissue_f1_mean"),
        )
        if a is not None
    ]
    ready = all(
        summary["arms"][k].get("matched_16ep")
        for k in (
            "N-light-gene-max",
            "N-light-gene-mean",
            "N-cascade-scalar-mean-max",
            "N-cascade-scalar-max-mean",
            "N-cascade-vector-mean-max",
        )
    )
    r4 = False
    if ready and promoted:
        r4 = max(promoted) < max(p2_f1, CLASSICAL_TISSUE) - 0.01 and not (r1 or r2 or r3)
    summary["rules"]["nothing_beats_p2_or_classical"] = {
        "fired": r4,
        "detail": "Stop architecture sweeps; start age-primary seed-gene experiment.",
        "best_promoted_f1": max(promoted) if promoted else None,
    }

    summary["screen_complete"] = ready
    if not ready:
        summary["next_gate"] = "matched_16ep_promotion_screen"
        summary["recommendation"] = "Wait for all 16-ep promotion folds before deciding."
    elif r1:
        summary["next_gate"] = "prefer_one_hop_max"
        summary["recommendation"] = "One-hop max is preferred smaller architecture."
    elif r2:
        summary["next_gate"] = "typed_rbs_aggregation"
        summary["recommendation"] = "Pursue typed-RBS aggregation (not scalar MBS)."
    elif r3:
        summary["next_gate"] = "retain_pooling_2x2"
        summary["recommendation"] = "Retain full 2×2 pooling result; no pooling lock."
    elif r4:
        summary["next_gate"] = "age_primary_seed_mask"
        summary["recommendation"] = "No promoted architecture beats P2-G/classical; start seed-mask."
    else:
        summary["next_gate"] = "manual_review"
        summary["recommendation"] = "Screen complete but no rule fired cleanly; review Pareto."

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, default=REPORT)
    args = parser.parse_args()
    report_dir = args.report_dir if args.report_dir.is_absolute() else ROOT / args.report_dir
    summary = decide(report_dir)
    out = report_dir / "promotion_decision.json"
    out.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
