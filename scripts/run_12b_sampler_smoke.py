#!/usr/bin/env python3
"""Milestone 12b §3.5 smoke: full-width + within-gene CpG sampler (max_cpgs=16).

Validates that capping edges per gene raises feasible batch vs dense full-width
(batch 55 at 440k edges). Not product OOF acceptance.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from run_7g_gene_only_probe import train_flat_region_arm  # noqa: E402

from mbs.paths import DataPaths  # noqa: E402

CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_sampler_smoke.yaml"
RUN_PREFIX = "stage0-12b-sampler-smoke"
REPORT = ROOT / "reports" / "inspection" / "stage0_12b_sampler_smoke"


def main() -> None:
    paths = DataPaths.from_environment()
    REPORT.mkdir(parents=True, exist_ok=True)
    rows = train_flat_region_arm(
        paths=paths,
        config_path=CONFIG,
        run_prefix=RUN_PREFIX,
        device="cuda",
        report_dir=REPORT,
        fold_filter=0,
    )
    print(f"[12b-sampler] done rows={rows!r}", flush=True)


if __name__ == "__main__":
    main()
