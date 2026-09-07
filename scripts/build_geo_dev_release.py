#!/usr/bin/env python3
"""Build immutable deepmat-data-geo-dev-v1 phenotype release (not ATS).

Arms: hub_only / geo_only / hub_geo / metadata_only.

Only samples with methylation (Hub membership and/or EWAS_db assay on disk)
and ≥1 acceptable observed label (age/sex/tissue). Disease/cancer columns may
appear for audit; ``disease_training_ok`` / ``cancer_training_ok`` stay false.

Does **not** mutate Hub packs or matrix-hub-age-tissue-sex-full-v1.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json
from mbs.geo_metadata import GEO_SOURCE_FAMILY
from mbs.paths import DataPaths
from mbs.release import HUB_FAMILIES

RELEASE_ID = "deepmat-data-geo-dev-v1"
ACCEPT = frozenset({"age", "sex", "tissue"})


def _blank(value: object) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return True
    return str(value).strip() == ""


def _hub_wide(data_root: Path) -> pd.DataFrame:
    """Union Hub pack sample_info tables into one wide phenotype frame."""
    rows: dict[str, dict[str, Any]] = {}
    for family in HUB_FAMILIES:
        path = data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
        if not path.is_file():
            continue
        frame = pd.read_parquet(path)
        for rec in frame.to_dict(orient="records"):
            sid = str(rec["sample_id"]).strip()
            if not sid:
                continue
            slot = rows.setdefault(
                sid,
                {
                    "sample_id": sid,
                    "study_id": str(rec.get("study_id") or ""),
                    "platform": rec.get("platform"),
                    "age": None,
                    "sex": None,
                    "tissue": None,
                    "disease": None,
                    "cancer": None,
                    "hub_families": set(),
                    "has_age": False,
                    "has_sex": False,
                    "has_tissue": False,
                },
            )
            slot["hub_families"].add(family)
            if _blank(slot["study_id"]) and not _blank(rec.get("study_id")):
                slot["study_id"] = str(rec["study_id"])
            if _blank(slot["platform"]) and not _blank(rec.get("platform")):
                slot["platform"] = rec.get("platform")

            # Primary pack value
            if family == "age" and not _blank(rec.get("phenotype_value_numeric")):
                slot["age"] = float(rec["phenotype_value_numeric"])
                slot["has_age"] = True
            elif family == "sex" and not _blank(rec.get("phenotype_value")):
                slot["sex"] = str(rec["phenotype_value"]).strip()
                slot["has_sex"] = True
            elif family == "tissue" and not _blank(rec.get("phenotype_value")):
                slot["tissue"] = str(rec["phenotype_value"]).strip()
                slot["has_tissue"] = True
            elif family == "disease" and not _blank(rec.get("phenotype_value")):
                slot["disease"] = str(rec["phenotype_value"]).strip()
            elif family == "cancer" and not _blank(rec.get("phenotype_value")):
                slot["cancer"] = str(rec["phenotype_value"]).strip()

            # Cross-pack aux columns on the same row
            if not slot["has_age"] and not _blank(rec.get("age")):
                try:
                    slot["age"] = float(rec["age"])
                    slot["has_age"] = True
                except (TypeError, ValueError):
                    pass
            if not slot["has_sex"] and not _blank(rec.get("sex")):
                slot["sex"] = str(rec["sex"]).strip()
                slot["has_sex"] = True
            if not slot["has_tissue"] and not _blank(rec.get("tissue")):
                slot["tissue"] = str(rec["tissue"]).strip()
                slot["has_tissue"] = True

    out = []
    for slot in rows.values():
        fams = sorted(slot.pop("hub_families"))
        slot["source_arm"] = "hub_only"
        slot["hub_families"] = json.dumps(fams)
        slot["has_acceptable_label"] = bool(
            slot["has_age"] or slot["has_sex"] or slot["has_tissue"]
        )
        slot["disease_training_ok"] = False
        slot["cancer_training_ok"] = False
        out.append(slot)
    return pd.DataFrame(out)


def _geo_wide(con: duckdb.DuckDBPyConnection, samples: pd.DataFrame) -> pd.DataFrame:
    pheno = con.execute(
        """
        SELECT sample_id, phenotype_id, label_status,
               categorical_value, numeric_value, ontology_id
        FROM sample_phenotype
        WHERE is_observed AND source_family = ?
        """,
        [GEO_SOURCE_FAMILY],
    ).fetchdf()
    if pheno.empty:
        return pd.DataFrame()

    meta = samples.set_index("sample_id", drop=False)
    rows: dict[str, dict[str, Any]] = {}
    for rec in pheno.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        slot = rows.setdefault(
            sid,
            {
                "sample_id": sid,
                "study_id": None,
                "platform": None,
                "age": None,
                "sex": None,
                "tissue": None,
                "tissue_ontology_id": None,
                "disease": None,
                "disease_label_status": None,
                "cancer": None,
                "cancer_label_status": None,
                "has_age": False,
                "has_sex": False,
                "has_tissue": False,
            },
        )
        pid = str(rec["phenotype_id"])
        if pid == "age" and rec.get("numeric_value") is not None:
            slot["age"] = float(rec["numeric_value"])
            slot["has_age"] = True
        elif pid == "sex" and not _blank(rec.get("categorical_value")):
            slot["sex"] = str(rec["categorical_value"])
            slot["has_sex"] = True
        elif pid == "tissue" and not _blank(rec.get("categorical_value")):
            slot["tissue"] = str(rec["categorical_value"])
            slot["tissue_ontology_id"] = rec.get("ontology_id")
            slot["has_tissue"] = True
        elif pid == "disease" and not _blank(rec.get("categorical_value")):
            slot["disease"] = str(rec["categorical_value"])
            slot["disease_label_status"] = str(rec.get("label_status") or "")
        elif pid == "cancer" and not _blank(rec.get("categorical_value")):
            slot["cancer"] = str(rec["categorical_value"])
            slot["cancer_label_status"] = str(rec.get("label_status") or "")

    out = []
    for sid, slot in rows.items():
        if sid in meta.index:
            m = meta.loc[sid]
            if isinstance(m, pd.DataFrame):
                m = m.iloc[0]
            slot["study_id"] = m.get("study_id")
            slot["platform"] = m.get("platform_id")
        slot["source_arm"] = "geo_only"
        slot["has_acceptable_label"] = bool(
            slot["has_age"] or slot["has_sex"] or slot["has_tissue"]
        )
        slot["disease_training_ok"] = False
        slot["cancer_training_ok"] = False
        out.append(slot)
    return pd.DataFrame(out)


def _merge_hub_geo(hub: pd.DataFrame, geo: pd.DataFrame) -> pd.DataFrame:
    """Union sample_ids; Hub labels win; GEO fills blanks."""
    hub_idx = hub.set_index("sample_id", drop=False) if not hub.empty else None
    geo_idx = geo.set_index("sample_id", drop=False) if not geo.empty else None
    ids = set()
    if hub_idx is not None:
        ids |= set(hub_idx.index.astype(str))
    if geo_idx is not None:
        ids |= set(geo_idx.index.astype(str))
    fill_keys = (
        "age",
        "sex",
        "tissue",
        "tissue_ontology_id",
        "disease",
        "disease_label_status",
        "cancer",
        "cancer_label_status",
        "study_id",
        "platform",
    )
    out: list[dict[str, Any]] = []
    for sid in sorted(ids):
        if hub_idx is not None and sid in hub_idx.index:
            base = hub_idx.loc[sid]
            if isinstance(base, pd.DataFrame):
                base = base.iloc[0]
            row = base.to_dict()
            row["source_arm"] = "hub_geo"
        else:
            base = geo_idx.loc[sid]  # type: ignore[union-attr]
            if isinstance(base, pd.DataFrame):
                base = base.iloc[0]
            row = base.to_dict()
            row["source_arm"] = "hub_geo"
        if geo_idx is not None and sid in geo_idx.index:
            g = geo_idx.loc[sid]
            if isinstance(g, pd.DataFrame):
                g = g.iloc[0]
            for key in fill_keys:
                if _blank(row.get(key)) and not _blank(g.get(key)):
                    row[key] = g.get(key)
            if not row.get("has_age") and not _blank(row.get("age")):
                row["has_age"] = True
            if not row.get("has_sex") and not _blank(row.get("sex")):
                row["has_sex"] = True
            if not row.get("has_tissue") and not _blank(row.get("tissue")):
                row["has_tissue"] = True
        row["has_acceptable_label"] = bool(
            row.get("has_age") or row.get("has_sex") or row.get("has_tissue")
        )
        row["disease_training_ok"] = False
        row["cancer_training_ok"] = False
        out.append(row)
    frame = pd.DataFrame(out)
    if frame.empty:
        return frame
    return frame.loc[frame["has_acceptable_label"]].reset_index(drop=True)


def main() -> None:
    paths = DataPaths.from_environment()
    catalog = (
        paths.data_root
        / "canonical"
        / "releases"
        / "deepmat-data-v1"
        / "catalog"
        / "catalog.duckdb"
    )
    if not catalog.is_file():
        raise SystemExit(f"missing catalog: {catalog}")

    con = duckdb.connect(str(catalog), read_only=True)
    try:
        assay = con.execute(
            """
            SELECT upper(regexp_extract(path, '(GSM[0-9]+)\\.txt$', 1)) AS sample_id
            FROM assay_file
            WHERE path LIKE '%/GSM%.txt'
              AND regexp_extract(path, '(GSM[0-9]+)\\.txt$', 1) != ''
            """
        ).fetchdf()
        hub_mem = con.execute(
            "SELECT DISTINCT sample_id FROM sample_source_membership"
        ).fetchdf()
        samples = con.execute(
            """
            SELECT s.sample_id, s.study_id, st.platform_id
            FROM sample s
            LEFT JOIN study st USING (study_id)
            WHERE s.sample_id LIKE 'GSM%'
            """
        ).fetchdf()
        geo = _geo_wide(con, samples)
    finally:
        con.close()

    assay_ids = set(assay["sample_id"].astype(str)) if not assay.empty else set()
    hub_ids = set(hub_mem["sample_id"].astype(str)) if not hub_mem.empty else set()

    hub = _hub_wide(paths.data_root)
    hub = hub.loc[hub["has_acceptable_label"]].reset_index(drop=True)
    # Hub packs imply methylation; still require membership id present.
    hub = hub.loc[hub["sample_id"].astype(str).isin(hub_ids)].reset_index(drop=True)

    geo_only = geo.copy() if not geo.empty else pd.DataFrame()
    if not geo_only.empty:
        geo_only = geo_only.loc[geo_only["has_acceptable_label"]]
        geo_only = geo_only.loc[geo_only["sample_id"].astype(str).isin(assay_ids)]
        geo_only = geo_only.loc[~geo_only["sample_id"].astype(str).isin(hub_ids)]
        geo_only = geo_only.reset_index(drop=True)

    # For hub_geo fill, allow GEO rows with assay OR hub (hub already has meth).
    geo_for_union = geo.copy() if not geo.empty else pd.DataFrame()
    if not geo_for_union.empty:
        geo_for_union = geo_for_union.loc[geo_for_union["has_acceptable_label"]]
        geo_for_union = geo_for_union.loc[
            geo_for_union["sample_id"].astype(str).isin(assay_ids | hub_ids)
        ].reset_index(drop=True)

    hub_geo = _merge_hub_geo(hub, geo_for_union)

    metadata_only = pd.DataFrame(
        [
            {
                "sample_id": r["sample_id"],
                "study_id": r.get("study_id"),
                "platform": r.get("platform"),
                "phenotype_value": r.get("study_id"),
                "phenotype_family": "metadata_only_confounding",
                "arm": "metadata_only",
            }
            for r in (hub_geo.to_dict(orient="records") if not hub_geo.empty else [])
        ]
    )

    root = paths.data_root / "canonical" / "releases" / RELEASE_ID
    pheno_dir = root / "phenotypes"
    pheno_dir.mkdir(parents=True, exist_ok=True)

    arms = {
        "hub_only": hub,
        "geo_only": geo_only,
        "hub_geo": hub_geo,
        "metadata_only": metadata_only,
    }
    checksums: dict[str, str] = {}
    counts: dict[str, int] = {}
    for name, frame in arms.items():
        out = pheno_dir / f"sample_info_{name}.parquet"
        frame.to_parquet(out, index=False)
        checksums[name] = sha256_file(out)
        counts[name] = len(frame)

    cfg_dir = paths.project_root / "configs" / "experiment" / "geo_dev"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    for arm in arms:
        (cfg_dir / f"{arm}.yaml").write_text(
            "\n".join(
                [
                    f"# GEO-dev comparison arm: {arm}",
                    f"# Release: {RELEASE_ID}",
                    "# Phenotype tables only — no training launch until gates allow.",
                    f"release_id: {RELEASE_ID}",
                    f"arm: {arm}",
                    "phenotype_table: "
                    f"canonical/releases/{RELEASE_ID}/phenotypes/"
                    f"sample_info_{arm}.parquet",
                    "mutate_ats: false",
                    "disease_training_ok: false",
                    "cancer_training_ok: false",
                    "",
                ]
            ),
            encoding="utf-8",
        )

    manifest = {
        "artifact_version": "geo-dev-release-v1",
        "release_id": RELEASE_ID,
        "created_at": utc_now_iso(),
        "parent_catalog_release": "deepmat-data-v1",
        "mutates_ats": False,
        "mutates_hub_packs": False,
        "n_assay_gsm": len(assay_ids),
        "n_hub_membership_gsm": len(hub_ids),
        "arm_counts": counts,
        "parquet_sha256": checksums,
        "acceptable_label_traits": sorted(ACCEPT),
        "notes": [
            "Phenotype-only development release (no matrix convert).",
            "GEO-only requires EWAS_db assay .txt on disk + acceptable age/sex/tissue.",
            "Disease/cancer training_ok flags are false in this builder.",
            "Do not enlarge matrix-hub-age-tissue-sex-full-v1.",
            "Larger GEO crawl and training launch remain gated.",
        ],
    }
    write_json(root / "manifest.json", manifest)

    report_dir = paths.project_root / "reports" / "inspection" / "deepmat_data_geo_dev_v1"
    report_dir.mkdir(parents=True, exist_ok=True)
    write_json(report_dir / "summary.json", manifest)
    lines = [
        f"# {RELEASE_ID}",
        "",
        f"- Generated: `{manifest['created_at']}`",
        "- Mutates ATS / Hub packs: **false** / **false**",
        "",
        "## Arm counts",
        "",
        "| arm | n_samples |",
        "| --- | ---: |",
    ]
    for arm, n in counts.items():
        lines.append(f"| `{arm}` | {n} |")
    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            f"- `{root / 'manifest.json'}`",
            f"- `{pheno_dir}/sample_info_*.parquet`",
            "- `configs/experiment/geo_dev/*.yaml` (stubs; no training launch)",
            "",
            "## Notes",
            "",
        ]
    )
    lines.extend(f"- {n}" for n in manifest["notes"])
    lines.append("")
    (report_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    summary_line = json.dumps(
        {"release_id": RELEASE_ID, "arm_counts": counts, "root": str(root)},
        indent=2,
    )
    print(summary_line)


if __name__ == "__main__":
    main()
