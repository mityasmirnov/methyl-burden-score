from __future__ import annotations

from pathlib import Path
from typing import Any

import pyarrow as pa

from mbs.inspect_cpgcorpus import inspect_cpgcorpus_gpl, write_cpgcorpus_report


def _write_arrow(path: Path, table: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        writer = pa.ipc.new_file(handle, table.schema)
        writer.write_table(table)
        writer.close()


def test_inspect_cpgcorpus_gpl_alignment_and_beta_qc(tmp_path: Path) -> None:
    root = tmp_path / "GSE1" / "GPL1"
    betas = pa.table(
        {
            "GSM_ID": ["GSM1", "GSM2"],
            "cg00000029": [0.1, 0.2],
            "cg00000109": [0.3, float("nan")],
        }
    )
    metadata = pa.table(
        {
            "GSM_ID": ["GSM1", "GSM2"],
            "platform_id": ["GPL1", "GPL1"],
            "organism_ch1": ["Homo sapiens", "Homo sapiens"],
            "Sex:ch1": ["F", "M"],
            "age (years):ch1": ["10", "20"],
        }
    )
    _write_arrow(root / "betas" / "QCDPB.arrow", betas)
    _write_arrow(root / "metadata" / "metadata.arrow", metadata)

    report = inspect_cpgcorpus_gpl(root, gse="GSE1", gpl="GPL1")
    assert report["sample_alignment"]["perfect_alignment"] is True
    assert report["value_qc"]["n_samples"] == 2
    assert report["value_qc"]["n_probes"] == 2
    assert report["value_qc"]["missing_fraction"] == 0.25
    assert report["metadata_counts"]["fields"]["platform_id"]["GPL1"] == 2
    assert report["metadata_counts"]["age_numeric"]["n_numeric"] == 2

    out = write_cpgcorpus_report(report, tmp_path / "report")
    assert (out / "summary.json").exists()
    assert (out / "summary.md").exists()
    assert (out / "warnings.json").exists()


def test_inspect_cpgcorpus_gpl_detects_misalignment(tmp_path: Path) -> None:
    root = tmp_path / "GSE1" / "GPL1"
    betas = pa.table({"GSM_ID": ["GSM1"], "cg00000029": [0.5]})
    metadata = pa.table(
        {
            "GSM_ID": ["GSM2"],
            "platform_id": ["GPL1"],
        }
    )
    _write_arrow(root / "betas" / "QCDPB.arrow", betas)
    _write_arrow(root / "metadata" / "metadata.arrow", metadata)

    report = inspect_cpgcorpus_gpl(root, gse="GSE1", gpl="GPL1")
    assert report["sample_alignment"]["perfect_alignment"] is False
    assert report["sample_alignment"]["beta_only_count"] == 1
    assert report["sample_alignment"]["metadata_only_count"] == 1
    assert any("aligned" in warning for warning in report["warnings"])
