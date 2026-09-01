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
