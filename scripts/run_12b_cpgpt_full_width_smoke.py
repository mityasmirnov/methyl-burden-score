#!/usr/bin/env python3
"""Milestone 12b probe: CpGPT2M static embeddings at full gene-linked width
(no 65,536-locus prefix) -- does the G2 signal seen at 65k hold, or does it
collapse the way the no-CpGPT full-width smoke did? See
stage0_12b_cpgpt_full_width_smoke.yaml for the regularization rationale.

Single fold (fold 0 of hub-nine-pack-3fold-v1), 16 epochs. Informal probe,
not 12b recipe acceptance (no within-gene CpG sampling cap yet).
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

CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_cpgpt_full_width_smoke.yaml"
RUN_PREFIX = "stage0-12b-cpgpt-full-width-smoke"
REPORT = ROOT / "reports" / "inspection" / "stage0_12b_gene_expansion_smoke"


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
    print(f"[12b-cpgpt-full-width-smoke] done rows={rows}", flush=True)


if __name__ == "__main__":
    main()
