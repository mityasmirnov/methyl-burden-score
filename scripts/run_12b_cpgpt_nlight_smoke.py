#!/usr/bin/env python3
"""Milestone 12b go/no-go smoke: does concatenating precomputed CpGPT2M
sequence-adapter embeddings onto the current best N-light architecture
(flat_region, m_only, 65,536-prefix) move mbs_e2e at all?

Single fold (fold 0 of hub-nine-pack-3fold-v1), 16 epochs -- same budget and
scale as stage0_7h_nine_pack_m_only_wide.yaml, isolating CpGPT as the only
variable. See stage0_12b_cpgpt_nlight_smoke.yaml for the full rationale and
the new flat_region_features.py plumbing this relies on.
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

CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_cpgpt_nlight_smoke.yaml"
RUN_PREFIX = "stage0-12b-cpgpt-nlight-smoke"
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
    print(f"[12b-cpgpt-smoke] done rows={rows}", flush=True)


if __name__ == "__main__":
    main()
