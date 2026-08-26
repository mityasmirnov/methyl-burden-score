#!/usr/bin/env python3
"""Merge parallel 7E bake-off report slices into stage0_7e_dev_cv/."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mbs.annotation.manifest import write_json  # noqa: E402
from mbs.training.dev_cv import _select_winner  # noqa: E402


def main() -> int:
    report_dirs = [
        ROOT / "reports/inspection/stage0_7e_dev_cv_flat",
        ROOT / "reports/inspection/stage0_7e_dev_cv_hier",
        ROOT / "reports/inspection/stage0_7e_dev_cv_mp",
    ]
    out = ROOT / "reports/inspection/stage0_7e_dev_cv"
    rows: list[dict] = []
    metadata: list[dict] = []
    fold_pack: dict = {}
    missing = []
    for d in report_dirs:
        summary_path = d / "summary.json"
        if not summary_path.is_file():
            missing.append(str(summary_path))
            continue
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        rows.extend(payload.get("results") or [])
        metadata.extend(payload.get("metadata_controls") or [])
        fold_pack.setdefault("split_id", payload.get("split_id"))
        fold_pack.setdefault("n_folds", payload.get("n_folds"))
        fold_pack.setdefault("sha256", payload.get("fold_sha256"))
    if missing:
        print("waiting on slices:\n  " + "\n  ".join(missing), file=sys.stderr)
        if not rows:
            return 1
    seen: set[tuple] = set()
    dedup: list[dict] = []
    for row in rows:
        key = (row.get("arm"), row.get("fold"), row.get("restart"))
        if key in seen:
            continue
        seen.add(key)
        dedup.append(row)
    seen_folds: set[int] = set()
    meta_dedup: list[dict] = []
    for mc in metadata:
        f = int(mc.get("fold", -1))
        if f in seen_folds:
            continue
        seen_folds.add(f)
        meta_dedup.append(mc)
    winner = _select_winner(dedup)
    summary = {
        "milestone": "7E",
        "split_id": fold_pack.get("split_id"),
        "n_folds": fold_pack.get("n_folds"),
        "fold_sha256": fold_pack.get("sha256"),
        "n_results": len(dedup),
        "winner": winner,
        "results": dedup,
        "metadata_controls": meta_dedup,
        "selection_rule": (
            "Among neural architecture arms (flat/hier/gene-direct/multipath): "
            "highest mean tissue macro-F1 on held-out studies, ties broken by "
            "lowest age RMSE. Transparent and metadata-only are ceilings."
        ),
        "cv_budget": {"max_loci": 8192, "max_epochs": 2},
        "merged_from": [str(d) for d in report_dirs if (d / "summary.json").is_file()],
    }
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "summary.json", summary)
    lines = [
        "# Milestone 7E development CV",
        "",
        f"Split: `{summary['split_id']}` ({summary['n_folds']} outer folds).",
        f"CV budget: max_loci={summary['cv_budget']['max_loci']}, "
        f"max_epochs={summary['cv_budget']['max_epochs']}.",
        "",
        f"**Winner for Milestone 7:** `{winner.get('arm')}` "
        f"(tissue macro-F1={winner.get('tissue_macro_f1')}, "
        f"age RMSE={winner.get('age_rmse')}).",
        "",
        f"Selection: {summary['selection_rule']}",
        "",
        "## Arms",
        "",
    ]
    lines.extend(
        [
            f"- `{row.get('arm')}` fold={row.get('fold')} restart={row.get('restart')} "
            f"L1={row.get('level1')} CpGPT={row.get('cpgpt')}: "
            f"tissue_f1={row.get('tissue_macro_f1')} age_rmse={row.get('age_rmse')}"
            for row in dedup
        ]
    )
    lines.extend(["", "## Metadata-only controls", ""])
    lines.extend(
        f"- fold={mc.get('fold')}: `{json.dumps(mc.get('metrics'), default=str)}`"
        for mc in meta_dedup
    )
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(winner, indent=2))
    print(f"wrote {out / 'summary.md'} ({len(dedup)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
