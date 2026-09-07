#!/usr/bin/env python3
"""Write disease/cancer label_status into the nine-pack phenotype table.

Joins Hub ``*_sample_info.parquet`` ``sample_type`` through
``SAMPLE_TYPE_CASE_CONTROL``. Does **not** treat pack-membership masks as
controls. Blood/brain are left deferred (no case/control recoverable).

Default writes beside the live table as
``sample_phenotype_table_hub_nine_pack_v1.label_status.parquet`` and a compact
census under ``reports/inspection/stage0_7h_nine_pack_smoke/``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from mbs.datahub_census import SAMPLE_TYPE_CASE_CONTROL
from mbs.paths import DataPaths

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLE = (
    "canonical/phenotypes/sample_phenotype_table_hub_nine_pack_v1.parquet"
)
REPORT = ROOT / "reports" / "inspection" / "stage0_7h_nine_pack_smoke"


def _map_status(sample_type: object) -> str | None:
    if sample_type is None or (isinstance(sample_type, float) and pd.isna(sample_type)):
        return None
    return SAMPLE_TYPE_CASE_CONTROL.get(str(sample_type).strip().lower())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the live nine-pack phenotype parquet (makes a .bak first).",
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    table_path = paths.data_root / DEFAULT_TABLE
    disease_info = paths.data_root / "canonical/phenotypes/disease_sample_info.parquet"
    cancer_info = paths.data_root / "canonical/phenotypes/cancer_sample_info.parquet"
    df = pd.read_parquet(table_path)
    dis = pd.read_parquet(disease_info, columns=["sample_id", "sample_type"])
    can = pd.read_parquet(cancer_info, columns=["sample_id", "sample_type"])
    # disease_sample_info.parquet / cancer_sample_info.parquet have duplicate
    # sample_id rows (multiple source annotation passes). sample_type is
    # consistent across duplicates for every id (verified), but merging
    # without deduplicating first fans out rows -- confirmed this produced
    # 37,257 rows instead of 34,234 (2,011 duplicated sample_ids) before this
    # fix. Dedup first so the merge can't multiply the phenotype table.
    for name, frame in (("disease", dis), ("cancer", can)):
        n_types = frame.groupby("sample_id")["sample_type"].nunique()
        conflicting = n_types[n_types > 1]
        if not conflicting.empty:
            raise ValueError(
                f"{name}: {len(conflicting)} sample_id(s) have conflicting sample_type, "
                f"e.g. {conflicting.index[:5].tolist()} -- refusing to silently pick one"
            )
    dis = dis.drop_duplicates(subset="sample_id")
    can = can.drop_duplicates(subset="sample_id")
    dis = dis.rename(columns={"sample_type": "disease_sample_type"})
    can = can.rename(columns={"sample_type": "cancer_sample_type"})
    out = df.merge(dis, on="sample_id", how="left").merge(can, on="sample_id", how="left")
    if len(out) != len(df):
        raise RuntimeError(f"merge changed row count {len(df)} -> {len(out)} -- dedup failed")
    out["disease_label_status"] = out["disease_sample_type"].map(_map_status)
    out["cancer_label_status"] = out["cancer_sample_type"].map(_map_status)
    # Trainable binary masks: only confirmed case/control (exclude adjacent_normal).
    out["disease_case_control_mask"] = out["disease_label_status"].isin(["case", "control"])
    out["cancer_case_control_mask"] = out["cancer_label_status"].isin(["case", "control"])
    out["disease_is_case"] = out["disease_label_status"].eq("case")
    out["cancer_is_case"] = out["cancer_label_status"].eq("case")

    census = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_samples": int(len(out)),
        "disease": {
            "label_status": out["disease_label_status"].value_counts(dropna=False).to_dict(),
            "raw_sample_type": out["disease_sample_type"].value_counts(dropna=False).to_dict(),
            "n_case_control_mask": int(out["disease_case_control_mask"].sum()),
        },
        "cancer": {
            "label_status": out["cancer_label_status"].value_counts(dropna=False).to_dict(),
            "raw_sample_type": out["cancer_sample_type"].value_counts(dropna=False).to_dict(),
            "n_case_control_mask": int(out["cancer_case_control_mask"].sum()),
        },
        "note": (
            "Pack-membership disease_mask/cancer_mask are unchanged and must not be "
            "used as control labels. Blood/brain heads remain deferred."
        ),
    }
    REPORT.mkdir(parents=True, exist_ok=True)
    (REPORT / "label_status_census.json").write_text(json.dumps(census, indent=2) + "\n")
    lines = [
        "# Nine-pack disease/cancer label_status census",
        "",
        f"Generated: `{census['generated_at']}`",
        "",
        census["note"],
        "",
        "## Disease",
        f"- label_status: `{census['disease']['label_status']}`",
        f"- case/control trainable: **{census['disease']['n_case_control_mask']}**",
        "",
        "## Cancer",
        f"- label_status: `{census['cancer']['label_status']}`",
        f"- case/control trainable: **{census['cancer']['n_case_control_mask']}**",
        "",
    ]
    (REPORT / "label_status_census.md").write_text("\n".join(lines) + "\n")

    if args.in_place:
        bak = table_path.with_suffix(table_path.suffix + ".pre_label_status.bak")
        if not bak.is_file():
            table_path.replace(bak)
            print(f"backed up live table -> {bak}")
            out.to_parquet(table_path, index=False)
            dest = table_path
        else:
            out.to_parquet(table_path, index=False)
            dest = table_path
    else:
        dest = table_path.with_name(
            "sample_phenotype_table_hub_nine_pack_v1.label_status.parquet"
        )
        out.to_parquet(dest, index=False)
    print(f"wrote {dest}")
    print(f"wrote {REPORT / 'label_status_census.md'}")


if __name__ == "__main__":
    main()
