#!/usr/bin/env python3
"""Milestone 12b §3.5 smoke: full-width + sparse gather + within-gene sampler.

Default config is K=8 (`stage0_12b_sampler_smoke_k8.yaml`). Override with
MBS_SAMPLER_SMOKE_CONFIG or argv[1]. Not product OOF acceptance.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from run_7g_gene_only_probe import train_flat_region_arm  # noqa: E402

from mbs.paths import DataPaths  # noqa: E402

DEFAULT_CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_sampler_smoke_k8.yaml"
RUN_PREFIX = "stage0-12b-sampler-smoke-k8"
REPORT = ROOT / "reports" / "inspection" / "stage0_12b_sampler_smoke_k8"


def main() -> None:
    config = Path(
        sys.argv[1]
        if len(sys.argv) > 1
        else os.environ.get("MBS_SAMPLER_SMOKE_CONFIG", str(DEFAULT_CONFIG))
    )
    paths = DataPaths.from_environment()
    REPORT.mkdir(parents=True, exist_ok=True)
    rows = train_flat_region_arm(
        paths=paths,
        config_path=config,
        run_prefix=RUN_PREFIX,
        device="cuda",
        report_dir=REPORT,
        fold_filter=0,
    )
    print(f"[12b-sampler-k8] done rows={rows!r} config={config}", flush=True)


if __name__ == "__main__":
    main()
