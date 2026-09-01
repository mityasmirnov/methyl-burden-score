#!/usr/bin/env python3
"""Write 7G′ Stage B matched-panel report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mbs.annotation.manifest import write_json
from mbs.paths import DataPaths

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_prime_matched_probe"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    report_dir = args.report_dir
    if not report_dir.is_absolute():
        report_dir = paths.project_root / report_dir
    summary_path = report_dir / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    lines = [
        "# 7G′ Stage B — fold-selected panel + full model",
        "",
        "Arms: `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), "
        "`N-full`, `N-mbs-direct-only`.",
        "",
        f"Lock from Stage A: `{json.dumps(summary.get('lock_from_stage_a', {}), sort_keys=True)}`",
        "",
        "## Per-fold panels",
        "",
    ]
    for fold in summary.get("folds") or []:
        panel = fold.get("panel") or {}
        lines.append(
            f"- Fold {fold.get('fold')}: {panel.get('n_seed')} seeds → {panel.get('n_panel')} panel CpGs"
        )
    lines.extend(
        [
            "",
            "## Recommendation",
            "",
            "Compare `C-mvalue-enetS` vs `N-cascade-S` vs `N-light-type` tissue macro-F1 on "
            "identical fold-selected panels. Use `N-full` vs `N-mbs-direct-only` for orphan RBS.",
            "Product export includes `direct_cpg.zarr` from cascade score writes when direct "
            "loci are present.",
            "",
        ]
    )
    (report_dir / "analysis.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(report_dir / "analysis_manifest.json", {"status": "stage_b_report"})
    print(f"wrote {report_dir / 'analysis.md'}", flush=True)


if __name__ == "__main__":
    main()
