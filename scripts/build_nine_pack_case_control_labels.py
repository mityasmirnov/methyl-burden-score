#!/usr/bin/env python3
"""Add honest disease/cancer case-control labels to the nine-pack phenotype table.

Milestone 10c policy (see ``reports/inspection/stage0_7h_nine_pack_smoke/
trait_hygiene.md``): the nine-pack table's own ``disease_mask``/``cancer_mask``
columns mean "sample routed into that pack," not "confirmed control" when
false — training a case/control head on those directly would silently
conflate "unknown" with "control." The real label lives in each pack's own
``sample_type`` column (``disease_sample_info.parquet`` /
``cancer_sample_info.parquet``), mapped through the already-established
``mbs.datahub_census.SAMPLE_TYPE_CASE_CONTROL``.

This script does not touch GPU training or rewrite ``..._v1.parquet`` in
place -- it writes a new ``..._v2.parquet`` superset (all v1 columns plus
the new label/mask columns) so nothing that reads v1 is affected.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.datahub_census import SAMPLE_TYPE_CASE_CONTROL
from mbs.paths import DataPaths

FAMILIES = ("disease", "cancer")


def _case_control_labels(sample_info_path: Path, family: str) -> pd.DataFrame:
    """One row per sample_id: {family}_label_status in {case, control, adjacent_normal, None}."""
    info = pd.read_parquet(sample_info_path, columns=["sample_id", "sample_type"])
    # sample_id repeats within a pack (multi-row source annotations); verified
    # zero conflicting sample_type values across duplicates for disease/cancer
    # (checked interactively before writing this script) -- safe to take the
    # first non-null value per sample_id rather than silently picking blind.
    conflicts = info.groupby("sample_id")["sample_type"].nunique()
    bad = conflicts[conflicts > 1]
    if not bad.empty:
        raise ValueError(
            f"{family}: {len(bad)} sample_id(s) have conflicting sample_type values, "
            f"e.g. {bad.index[:5].tolist()} -- refusing to silently pick one"
        )
    dedup = info.drop_duplicates(subset="sample_id").set_index("sample_id")
    status = dedup["sample_type"].map(
        lambda v: SAMPLE_TYPE_CASE_CONTROL.get(str(v).strip().lower()) if pd.notna(v) else None
    )
    return status.rename(f"{family}_label_status").reset_index()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phenotype-table",
        default="canonical/phenotypes/sample_phenotype_table_hub_nine_pack_v1.parquet",
    )
    parser.add_argument(
        "--out-table",
        default="canonical/phenotypes/sample_phenotype_table_hub_nine_pack_v2.parquet",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    in_path = paths.data_root / args.phenotype_table
    out_path = paths.data_root / args.out_table
    if out_path.exists():
        raise SystemExit(f"{out_path} already exists -- refusing to overwrite")

    table = pd.read_parquet(in_path)
    print(f"[labels] loaded {len(table)} rows from {in_path}", flush=True)

    report: dict[str, object] = {"generated_at": utc_now_iso(), "families": {}}
    for family in FAMILIES:
        info_path = paths.data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
        labels = _case_control_labels(info_path, family)
        before = len(table)
        table = table.merge(labels, on="sample_id", how="left")
        if len(table) != before:
            raise RuntimeError(f"{family} merge changed row count {before} -> {len(table)}")

        status_col = f"{family}_label_status"
        mask_col = f"{family}_case_control_mask"
        table[mask_col] = table[status_col].isin(["case", "control"])

        pack_mask_col = f"{family}_mask"
        # Every row the old pack_mask marked True must have gotten a real label --
        # else the join silently dropped samples the pack-membership mask says exist.
        in_pack = table[pack_mask_col].fillna(False).astype(bool)
        missing_label = in_pack & table[status_col].isna()
        if missing_label.any():
            raise RuntimeError(
                f"{family}: {int(missing_label.sum())} pack-member samples got no "
                f"sample_type label after merge -- join is incomplete"
            )

        counts = table.loc[in_pack, status_col].value_counts(dropna=False).to_dict()
        n_usable = int(table[mask_col].sum())
        report["families"][family] = {  # type: ignore[index]
            "n_in_pack": int(in_pack.sum()),
            "status_counts": {str(k): int(v) for k, v in counts.items()},
            "n_case_or_control": n_usable,
        }
        print(f"[labels] {family}: n_in_pack={int(in_pack.sum())} usable(case+control)={n_usable}", flush=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out_path, index=False)
    report["n_rows"] = len(table)
    report["n_cols"] = len(table.columns)
    report["source_table"] = str(in_path)
    report["source_table_sha256"] = sha256_file(in_path)
    report["out_table"] = str(out_path)
    report["out_table_sha256"] = sha256_file(out_path)
    report["git_commit"] = git_commit(paths.project_root)
    manifest_path = out_path.with_suffix(".manifest.json")
    write_json(manifest_path, report)
    print(f"[labels] wrote {out_path} ({len(table)} rows, {len(table.columns)} cols)", flush=True)
    print(f"[labels] wrote {manifest_path}", flush=True)


if __name__ == "__main__":
    main()
