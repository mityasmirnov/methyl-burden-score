#!/usr/bin/env python3
"""Build full-catalog sample overview (Hub ∪ EWAS_db) parquet + CSV + manifest.

One row per catalog ``sample_id`` (~173k). Does not mutate ATS freeze artifacts.
See ``docs/plans/label-prep-bmi-ancestry-sample-overview.md``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.paths import DataPaths
from mbs.sample_overview_enrich import fill_species_for_rows, infer_platforms_for_rows

HUB_FAMILIES = (
    "age",
    "tissue",
    "sex",
    "disease",
    "cancer",
    "blood",
    "brain",
    "bmi",
    "ancestry",
)

OVERVIEW_STEM = "sample_overview_hub_geo_v1"


def _first_non_null(*values: object) -> object | None:
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and pd.isna(value):
            continue
        if isinstance(value, str) and value.strip() == "":
            continue
        return value
    return None


def build_sample_overview(
    *,
    data_root: Path,
    catalog_release: str = "deepmat-data-v1",
) -> pd.DataFrame:
    """Join catalog + lane flags + Hub packs + GEO into one wide overview."""
    data_root = data_root.resolve()
    tables = data_root / "canonical" / "releases" / catalog_release / "catalog" / "tables"
    phenotypes = data_root / "canonical" / "phenotypes"

    sample = pd.read_parquet(tables / "sample.parquet")
    study = pd.read_parquet(tables / "study.parquet")
    flags = pd.read_parquet(tables / "sample_lane_flags.parquet")
    membership = pd.read_parquet(tables / "sample_source_membership.parquet")

    study_slim = study[
        ["study_id", "platform_id", "processing_level", "genome_build", "gse_id"]
    ].rename(
        columns={
            "platform_id": "catalog_platform_id",
            "processing_level": "processing_level",
        }
    )

    out = sample.merge(flags, on=["sample_id", "study_id"], how="left", suffixes=("", "_flag"))
    out = out.merge(study_slim, on="study_id", how="left")

    # Per-pack membership booleans + primary matrix_id (first seen).
    memb = membership.copy()
    memb["phenotype_family"] = memb["phenotype_family"].astype(str)
    for family in HUB_FAMILIES:
        members = set(memb.loc[memb["phenotype_family"] == family, "sample_id"].astype(str))
        out[f"in_hub_{family}"] = out["sample_id"].astype(str).isin(members)

    matrix_by = (
        memb.dropna(subset=["matrix_id"])
        .drop_duplicates("sample_id", keep="first")
        .set_index("sample_id")["matrix_id"]
        .astype(str)
        .to_dict()
    )
    out["hub_matrix_id"] = out["sample_id"].map(matrix_by)

    # Nine-pack Hub phenotype labels (Hub-wins when present).
    hub_wide_path = phenotypes / "sample_phenotype_table_hub_nine_pack_v1.parquet"
    hub_cols = [
        "sample_id",
        "age_years",
        "age_mask",
        "tissue_label",
        "tissue_mask",
        "sex_label",
        "sex_mask",
        "bmi",
        "bmi_mask",
        "ancestry_label",
        "ancestry_mask",
        "brain_label",
        "brain_mask",
        "disease_mask",
        "cancer_mask",
        "platform_id",
    ]
    if hub_wide_path.is_file():
        hub_wide = pd.read_parquet(hub_wide_path, columns=hub_cols)
        hub_wide = hub_wide.rename(
            columns={
                "age_years": "hub_age_years",
                "age_mask": "hub_age_mask",
                "tissue_label": "hub_tissue_label",
                "tissue_mask": "hub_tissue_mask",
                "sex_label": "hub_sex_label",
                "sex_mask": "hub_sex_mask",
                "bmi": "hub_bmi",
                "bmi_mask": "hub_bmi_mask",
                "ancestry_label": "hub_ancestry_label",
                "ancestry_mask": "hub_ancestry_mask",
                "brain_label": "hub_brain_label",
                "brain_mask": "hub_brain_mask",
                "disease_mask": "hub_disease_mask",
                "cancer_mask": "hub_cancer_mask",
                "platform_id": "hub_platform_id",
            }
        )
        out = out.merge(hub_wide, on="sample_id", how="left")
    else:
        for col in (
            "hub_age_years",
            "hub_age_mask",
            "hub_tissue_label",
            "hub_tissue_mask",
            "hub_sex_label",
            "hub_sex_mask",
            "hub_bmi",
            "hub_bmi_mask",
            "hub_ancestry_label",
            "hub_ancestry_mask",
            "hub_brain_label",
            "hub_brain_mask",
            "hub_disease_mask",
            "hub_cancer_mask",
            "hub_platform_id",
        ):
            out[col] = None if "mask" not in col else False

    # GEO sample + series enrichment.
    geo_path = phenotypes / "geo_sample_metadata.parquet"
    if geo_path.is_file():
        geo = pd.read_parquet(
            geo_path,
            columns=[
                "sample_id",
                "platform_id",
                "catalog_platform_id",
                "pubmed_ids",
                "age",
                "sex",
                "tissue",
                "tissue_raw",
                "species_status",
                "taxon_id",
                "organism",
                "disease_label_status",
                "cancer_label_status",
            ],
        ).rename(
            columns={
                "platform_id": "geo_gpl_id",
                "catalog_platform_id": "geo_catalog_platform_id",
                "pubmed_ids": "geo_sample_pubmed_ids",
                "age": "geo_age",
                "sex": "geo_sex",
                "tissue": "geo_tissue",
                "tissue_raw": "geo_tissue_raw",
                "species_status": "species_status",
                "taxon_id": "taxon_id",
                "organism": "organism",
                "disease_label_status": "geo_disease_label_status",
                "cancer_label_status": "geo_cancer_label_status",
            }
        )
        out = out.merge(geo, on="sample_id", how="left")
    else:
        for col in (
            "geo_gpl_id",
            "geo_catalog_platform_id",
            "geo_sample_pubmed_ids",
            "geo_age",
            "geo_sex",
            "geo_tissue",
            "geo_tissue_raw",
            "species_status",
            "taxon_id",
            "organism",
            "geo_disease_label_status",
            "geo_cancer_label_status",
        ):
            out[col] = None

    series_path = phenotypes / "geo_series_metadata.parquet"
    if series_path.is_file():
        series = pd.read_parquet(
            series_path,
            columns=["study_id", "title", "summary", "pubmed_ids", "platform_title", "taxon"],
        ).rename(
            columns={
                "title": "series_title",
                "summary": "series_summary",
                "pubmed_ids": "geo_series_pubmed_ids",
                "platform_title": "geo_platform_title",
                "taxon": "series_taxon",
            }
        )
        out = out.merge(series, on="study_id", how="left")
    else:
        out["series_title"] = None
        out["series_summary"] = None
        out["geo_series_pubmed_ids"] = None
        out["geo_platform_title"] = None
        out["series_taxon"] = None

    # Blood cell composition (sparse Hub/DataHub census strings).
    census_path = phenotypes / "ewas_datahub_sample_census.parquet"
    if census_path.is_file():
        census = pd.read_parquet(census_path, columns=["sample_id", "cell component"])
        census = census.rename(columns={"cell component": "cell_component"})
        census = census.drop_duplicates("sample_id", keep="first")
        out = out.merge(census, on="sample_id", how="left")
    else:
        out["cell_component"] = None

    # Harmonized display phenotypes: Hub wins, GEO fills blanks.
    out["age_years"] = [
        _first_non_null(h, g, c)
        for h, g, c in zip(out["hub_age_years"], out["geo_age"], out["age"], strict=True)
    ]
    out["sex"] = [
        _first_non_null(h, g, c)
        for h, g, c in zip(out["hub_sex_label"], out["geo_sex"], out["sex"], strict=True)
    ]
    out["tissue"] = [
        _first_non_null(h, g, c)
        for h, g, c in zip(
            out["hub_tissue_label"], out["geo_tissue"], out["tissue_raw"], strict=True
        )
    ]
    out["bmi"] = out["hub_bmi"]
    out["ancestry_label"] = out["hub_ancestry_label"]
    out["brain_label"] = out["hub_brain_label"]
    out["platform_id"] = [
        _first_non_null(h, g, c)
        for h, g, c in zip(
            out["hub_platform_id"],
            out["geo_catalog_platform_id"],
            out["catalog_platform_id"],
            strict=True,
        )
    ]
    out["pubmed_ids"] = [
        _first_non_null(s, r)
        for s, r in zip(out["geo_sample_pubmed_ids"], out["geo_series_pubmed_ids"], strict=True)
    ]
    out["cpgpt_sample_embedding_status"] = "unavailable"
    out["in_geo_soft"] = out["geo_gpl_id"].notna() | out["species_status"].notna()

    # Fill platform / species gaps (DataHub JSON, GPL, n_probes, series taxon, Hub priors).
    out = infer_platforms_for_rows(out, data_root=data_root)
    out = fill_species_for_rows(out)

    # Truncate long free text for CSV friendliness (parquet keeps full via copy).
    if "series_summary" in out.columns:
        out["series_summary_short"] = out["series_summary"].map(
            lambda x: None
            if x is None or (isinstance(x, float) and pd.isna(x))
            else (str(x)[:500] + ("…" if len(str(x)) > 500 else ""))
        )

    # Stable column order (identity first).
    preferred = [
        "sample_id",
        "study_id",
        "gse_id",
        "source_sample_id",
        "donor_id",
        "replicate_group",
        "in_hub_baseline",
        "in_ewas_db",
        "in_geo_soft",
        "hub_families",
        *[f"in_hub_{f}" for f in HUB_FAMILIES],
        "hub_matrix_id",
        "platform_id",
        "platform_source",
        "n_probes",
        "catalog_platform_id",
        "geo_gpl_id",
        "geo_catalog_platform_id",
        "geo_platform_title",
        "processing_level",
        "genome_build",
        "age_years",
        "sex",
        "tissue",
        "bmi",
        "ancestry_label",
        "brain_label",
        "cell_component",
        "case_control",
        "hub_age_mask",
        "hub_tissue_mask",
        "hub_sex_mask",
        "hub_bmi_mask",
        "hub_ancestry_mask",
        "hub_brain_mask",
        "hub_disease_mask",
        "hub_cancer_mask",
        "species_status",
        "species_source",
        "taxon_id",
        "organism",
        "series_taxon",
        "pubmed_ids",
        "geo_sample_pubmed_ids",
        "geo_series_pubmed_ids",
        "series_title",
        "series_summary_short",
        "cpgpt_sample_embedding_status",
    ]
    rest = [c for c in out.columns if c not in preferred]
    # Drop bulky raw duplicates from export spine (keep short summary).
    drop = {"series_summary", "metadata_json", "characteristics_raw"}
    cols = [c for c in preferred + rest if c in out.columns and c not in drop]
    return out.loc[:, cols]


def write_overview_artifacts(
    frame: pd.DataFrame,
    *,
    out_dir: Path,
    stem: str = OVERVIEW_STEM,
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = out_dir / f"{stem}.parquet"
    csv_path = out_dir / f"{stem}.csv.gz"
    manifest_path = out_dir / f"{stem}.manifest.json"
    frame.to_parquet(parquet_path, index=False)
    frame.to_csv(csv_path, index=False, compression="gzip")
    payload = {
        "artifact": stem,
        "created_at": utc_now_iso(),
        "git_commit": git_commit(Path(__file__).resolve().parents[1]),
        "n_rows": len(frame),
        "n_columns": int(frame.shape[1]),
        "parquet": str(parquet_path),
        "csv_gz": str(csv_path),
        "parquet_sha256": sha256_file(parquet_path),
        "csv_gz_sha256": sha256_file(csv_path),
        "cpgpt_sample_embedding_status": "unavailable",
        "notes": (
            "GSM=sample_id, GSE=study_id, GPL=geo_gpl_id. "
            "Hub phenotype columns win; GEO fills blanks. "
            "CpGPT sample CLS embeddings are not exported (locus-only today)."
        ),
    }
    write_json(manifest_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog-release",
        default="deepmat-data-v1",
        help="Catalog release id under canonical/releases/",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Default: $MBS_DATA_ROOT/canonical/phenotypes",
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    out_dir = args.out_dir or (paths.data_root / "canonical" / "phenotypes")
    frame = build_sample_overview(
        data_root=paths.data_root,
        catalog_release=args.catalog_release,
    )
    if len(frame) == 0:
        raise SystemExit("overview frame is empty")
    payload = write_overview_artifacts(frame, out_dir=out_dir)
    print(json.dumps(payload, indent=2, default=str))  # noqa: T201


if __name__ == "__main__":
    main()
