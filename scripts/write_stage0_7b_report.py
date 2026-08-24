#!/usr/bin/env python3
"""Write Milestone 7B Hub full-pack matrix inspection report."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from mbs.matrix.hub_pack import SUPPORTED_PACK_FAMILIES
from mbs.paths import DataPaths

SEVEN_B_FAMILIES = ("disease", "cancer", "blood", "brain", "bmi", "ancestry")
FROZEN_FAMILIES = ("age", "tissue", "sex")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _matrix_row(data_root: Path, family: str) -> dict:
    matrix_id = f"matrix-hub-{family}-full-v1"
    root = data_root / "canonical" / "matrices" / matrix_id
    man_path = root / "matrix_manifest.json"
    stats_path = root / "conversion_stats.json"
    if not man_path.is_file():
        return {"family": family, "matrix_id": matrix_id, "present": False}
    man = _load_json(man_path)
    stats = _load_json(stats_path) if stats_path.is_file() else {}
    source = (man.get("source_files") or [{}])[0]
    sha = str(source.get("sha256") or "")
    name_size = f"{Path(str(source.get('path') or '')).name}:{source.get('byte_size')}"
    return {
        "family": family,
        "matrix_id": matrix_id,
        "present": True,
        "shape": man.get("shape"),
        "n_samples": (man.get("shape") or [None])[0],
        "n_loci": (man.get("shape") or [None, None])[1],
        "n_phenotype_rows": stats.get("n_phenotype_rows"),
        "n_collapsed_probes": stats.get("n_collapsed_probes"),
        "platform_id": man.get("platform_id"),
        "compression": man.get("compression"),
        "source_sha256": sha,
        "source_byte_size": source.get("byte_size"),
        "checksum_is_content": bool(sha) and name_size not in sha,
        "created_at": man.get("created_at"),
    }


def main() -> None:
    paths = DataPaths.from_environment()
    data_root = paths.data_root
    report_dir = (
        Path(__file__).resolve().parents[1] / "reports" / "inspection" / "stage0_7b_hub_matrices"
    )
    report_dir.mkdir(parents=True, exist_ok=True)

    seven_b = [_matrix_row(data_root, f) for f in SEVEN_B_FAMILIES]
    frozen = [_matrix_row(data_root, f) for f in FROZEN_FAMILIES]
    all_full = [_matrix_row(data_root, f) for f in SUPPORTED_PACK_FAMILIES]

    index_path = data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
    overlap_path = report_dir / "overlap_concordance.json"
    index_meta: dict = {"present": index_path.is_file()}
    if index_path.is_file():
        index = pd.read_parquet(index_path)
        index_meta.update(
            {
                "n_rows": len(index),
                "n_families": int(index["family"].nunique()) if not index.empty else 0,
                "n_unique_gsm": int(index["sample_id"].nunique()) if not index.empty else 0,
            }
        )
    overlap = _load_json(overlap_path) if overlap_path.is_file() else None

    payload = {
        "milestone": "7B",
        "seven_b_matrices": seven_b,
        "frozen_age_tissue_sex": frozen,
        "all_full_matrices": all_full,
        "virtual_index": index_meta,
        "overlap": overlap,
        "all_seven_b_present": all(row["present"] for row in seven_b),
    }
    (report_dir / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Milestone 7B: complete Hub pack matrices",
        "",
        "Per-pack Zarr stores (no dense nine-pack union). Age/tissue/sex full",
        "matrices were not reconverted.",
        "",
        "## 7B matrices",
        "",
        "| Family | Matrix ID | Samples | Loci | Phenotype rows | Platform | Content sha256 |",
        "|--------|-----------|--------:|-----:|---------------:|----------|----------------|",
    ]
    for row in seven_b:
        if not row["present"]:
            lines.append(f"| `{row['family']}` | `{row['matrix_id']}` | — | — | — | — | missing |")
            continue
        sha = str(row["source_sha256"])[:12]
        lines.append(
            f"| `{row['family']}` | `{row['matrix_id']}` | {row['n_samples']} | "
            f"{row['n_loci']} | {row['n_phenotype_rows']} | `{row['platform_id']}` | "
            f"`{sha}…` |"
        )
    lines.extend(
        [
            "",
            "Disease/cancer: `n_samples` is unique GSM (matrix rows);",
            "`n_phenotype_rows` is long-form sample-info (may exceed unique GSM).",
            "",
            "## Virtual index",
            "",
            "- Path: `canonical/matrices/hub_pack_matrix_index.parquet`",
            f"- Present: {index_meta.get('present')}",
            f"- Rows / unique GSM / families: "
            f"{index_meta.get('n_rows')} / {index_meta.get('n_unique_gsm')} / "
            f"{index_meta.get('n_families')}",
            "",
            "## Overlap concordance",
            "",
        ]
    )
    if overlap:
        lines.append(f"- Status: `{overlap.get('status')}`")
        lines.append(f"- Shared GSM: {overlap.get('n_shared_gsm')}")
        n_pairs = overlap.get("n_pairs_checked")
        n_disc = overlap.get("n_discordant")
        lines.append(f"- Pairs checked / discordant: {n_pairs} / {n_disc}")
        lines.append(f"- Max abs diff: {overlap.get('max_abs_diff')}")
        lines.append(f"- Merge allowed: {overlap.get('merge_allowed')}")
    else:
        lines.append("- Overlap report not written yet (`overlap_concordance.json`).")
    lines.extend(
        [
            "",
            "## Frozen 5d matrices (not overwritten)",
            "",
        ]
    )
    lines.extend(
        f"- `{row['matrix_id']}` shape `{row['shape']}` platform `{row['platform_id']}`"
        for row in frozen
        if row["present"]
    )
    (report_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    sys.stdout.write(
        json.dumps(
            {"report_dir": str(report_dir), "all_seven_b_present": payload["all_seven_b_present"]}
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
