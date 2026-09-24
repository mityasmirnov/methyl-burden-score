#!/usr/bin/env python3
"""Milestone 12b go/no-go smoke: N-light on the full gene-linked column
universe (no 65,536-prefix cap) vs. the Milestone-10 65k-prefix baseline.

Single fold (fold 0 of hub-nine-pack-3fold-v1), 16 epochs -- same budget as
stage0_7h_nine_pack_m_only_wide.yaml. Not the full 12b recipe (no within-gene
CpG sampling cap, no role embeddings); just answers whether wider gene
coverage moves mbs_e2e at all before investing in that build.
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

CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_gene_expansion_nlight_smoke.yaml"
RUN_PREFIX = "stage0-12b-gene-expansion-nlight-smoke"
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
    print(f"[12b-nlight-smoke] done rows={rows}", flush=True)


if __name__ == "__main__":
    main()
