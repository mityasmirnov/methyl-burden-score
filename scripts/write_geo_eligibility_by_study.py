#!/usr/bin/env python3
"""GEO eligibility-by-study deep audit (batch-50 / geo_metadata_backfill).

Reports trait × study coverage, Hub overlap, tissue map status, and which GSM
have EWAS_db assay files on disk. Catalog-global eligibility is not enough to train.
"""

from __future__ import annotations

from pathlib import Path

import duckdb

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.geo_metadata import GEO_SOURCE_FAMILY, load_geo_frame
from mbs.paths import DataPaths

REPORT_DIR = Path("reports/inspection/deepmat_data_v1/geo_backfill_batch")


def main() -> None:
    paths = DataPaths.from_environment()
    report_dir = paths.project_root / REPORT_DIR
    report_dir.mkdir(parents=True, exist_ok=True)
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

    geo = load_geo_frame(paths.data_root)
    tissue_by_study: dict[str, dict[str, int]] = {}
    if not geo.empty and "tissue_map_status" in geo.columns:
        for study_id, grp in geo.groupby(geo["study_id"].astype(str), sort=True):
            counts = grp["tissue_map_status"].fillna("empty").value_counts().to_dict()
            tissue_by_study[str(study_id)] = {str(k): int(v) for k, v in counts.items()}

    con = duckdb.connect(str(catalog), read_only=True)
    try:
        # Assay GSM on disk (EWAS_db beta txt).
        assay = con.execute(
            """
            SELECT upper(regexp_extract(path, '(GSM[0-9]+)\\.txt$', 1)) AS sample_id,
                   study_id,
                   path
            FROM assay_file
            WHERE path LIKE '%/GSM%.txt'
              AND regexp_extract(path, '(GSM[0-9]+)\\.txt$', 1) IS NOT NULL
              AND regexp_extract(path, '(GSM[0-9]+)\\.txt$', 1) != ''
            """
        ).fetchdf()
        hub = con.execute(
            """
            SELECT DISTINCT sample_id
            FROM sample_source_membership
            """
        ).fetchdf()
        pheno = con.execute(
            """
            SELECT sample_id, phenotype_id, label_status, is_observed,
                   categorical_value, numeric_value, ontology_id
            FROM sample_phenotype
            WHERE source_family = ?
              AND is_observed
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchdf()
        sample_study = con.execute(
            """
            SELECT sample_id, study_id
            FROM sample
            WHERE sample_id LIKE 'GSM%'
            """
        ).fetchdf()
        global_elig = con.execute(
            """
            SELECT phenotype_id, n_samples, n_cases, n_controls,
                   eligible_core_task, eligible_auxiliary_task, exclusion_reason
            FROM trait_eligibility
            WHERE phenotype_family = ?
            ORDER BY phenotype_id
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchdf()
    finally:
        con.close()

    hub_ids = set(hub["sample_id"].astype(str)) if not hub.empty else set()
    assay_ids = set(assay["sample_id"].astype(str)) if not assay.empty else set()
    study_of = dict(
        zip(
            sample_study["sample_id"].astype(str),
            sample_study["study_id"].astype(str),
            strict=False,
        )
    )

    # Attach study_id to phenotype rows (prefer catalog sample.study_id).
    pheno = pheno.copy()
    pheno["sample_id"] = pheno["sample_id"].astype(str)
    pheno["study_id"] = pheno["sample_id"].map(study_of)
    pheno["has_hub_membership"] = pheno["sample_id"].isin(hub_ids)
    pheno["has_ewas_assay"] = pheno["sample_id"].isin(assay_ids)

    per_study: list[dict[str, object]] = []
    for study_id, grp in pheno.groupby(pheno["study_id"].fillna("UNKNOWN").astype(str), sort=True):
        gsm = set(grp["sample_id"])
        row: dict[str, object] = {
            "study_id": study_id,
            "n_geo_phenotype_rows": int(len(grp)),
            "n_geo_gsm": int(len(gsm)),
            "n_with_ewas_assay": int(sum(1 for s in gsm if s in assay_ids)),
            "n_with_hub_membership": int(sum(1 for s in gsm if s in hub_ids)),
            "n_methylation_and_label": int(
                sum(1 for s in gsm if s in assay_ids or s in hub_ids)
            ),
            "traits": {},
            "tissue_map_parquet": tissue_by_study.get(str(study_id), {}),
        }
        traits: dict[str, dict[str, object]] = {}
        for pid, pgrp in grp.groupby("phenotype_id"):
            status = pgrp["label_status"].fillna("observed").astype(str)
            traits[str(pid)] = {
                "n_rows": int(len(pgrp)),
                "n_gsm": int(pgrp["sample_id"].nunique()),
                "n_cases": int((status == "case").sum()),
                "n_controls": int((status == "control").sum()),
                "n_with_assay": int(pgrp["has_ewas_assay"].sum()),
                "n_hub_overlap_rows": int(pgrp["has_hub_membership"].sum()),
            }
        row["traits"] = traits
        per_study.append(row)

    # Training-candidate GSM: GEO-observed label + methylation on disk (assay or Hub).
    candidates = pheno.loc[pheno["has_ewas_assay"] | pheno["has_hub_membership"]]
    by_trait_candidates = {
        str(pid): {
            "n_gsm": int(g["sample_id"].nunique()),
            "n_geo_only_gsm": int(
                g.loc[~g["has_hub_membership"], "sample_id"].nunique()
            ),
            "n_hub_overlap_gsm": int(
                g.loc[g["has_hub_membership"], "sample_id"].nunique()
            ),
        }
        for pid, g in candidates.groupby("phenotype_id")
    }

    summary = {
        "generated_at": utc_now_iso(),
        "geo_source_family": GEO_SOURCE_FAMILY,
        "n_geo_observed_rows": int(len(pheno)),
        "n_geo_observed_gsm": int(pheno["sample_id"].nunique()),
        "n_assay_gsm_on_disk": int(len(assay_ids)),
        "n_hub_membership_gsm": int(len(hub_ids)),
        "n_geo_gsm_with_assay": int(pheno.loc[pheno["has_ewas_assay"], "sample_id"].nunique()),
        "n_geo_gsm_hub_overlap": int(
            pheno.loc[pheno["has_hub_membership"], "sample_id"].nunique()
        ),
        "global_trait_eligibility": global_elig.to_dict(orient="records")
        if not global_elig.empty
        else [],
        "candidate_gsm_by_trait": by_trait_candidates,
        "per_study": per_study,
        "notes": [
            "Global trait_eligibility can look green while many studies lack assay files or case/control balance.",
            "Training requires a separate deepmat-data-geo-dev-v1 release; do not mutate ATS.",
            "Disease/cancer: inspect per-study case/control before any head; diagnosis-free ≠ control.",
            "Larger GEO crawl and wiring GEO into age/tissue training remain gated.",
        ],
    }
    write_json(report_dir / "eligibility_by_study.json", summary)

    lines = [
        "# GEO eligibility by study (deep audit)",
        "",
        f"- Generated: `{summary['generated_at']}`",
        f"- GEO observed rows / GSM: **{summary['n_geo_observed_rows']}** / "
        f"**{summary['n_geo_observed_gsm']}**",
        f"- EWAS_db assay GSM on disk: **{summary['n_assay_gsm_on_disk']}**",
        f"- GEO GSM with assay: **{summary['n_geo_gsm_with_assay']}**",
        f"- GEO GSM with Hub membership (should be omitted from GEO-only): "
        f"**{summary['n_geo_gsm_hub_overlap']}**",
        "",
        "## Global catalog eligibility (`geo_metadata_backfill`)",
        "",
        "| phenotype_id | n | cases | controls | core | aux |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    for rec in summary["global_trait_eligibility"]:
        lines.append(
            f"| `{rec['phenotype_id']}` | {rec['n_samples']} | {rec.get('n_cases')} | "
            f"{rec.get('n_controls')} | {rec['eligible_core_task']} | "
            f"{rec['eligible_auxiliary_task']} |"
        )
        if rec.get("exclusion_reason") and not rec.get("eligible_core_task"):
            lines.append(f"  - {rec['exclusion_reason']}")
    lines.extend(
        [
            "",
            "## Candidate GSM by trait (GEO label + assay or Hub)",
            "",
            "| phenotype_id | n_gsm | geo_only | hub_overlap |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for pid, rec in sorted(by_trait_candidates.items()):
        lines.append(
            f"| `{pid}` | {rec['n_gsm']} | {rec['n_geo_only_gsm']} | {rec['n_hub_overlap_gsm']} |"
        )
    lines.extend(
        [
            "",
            "## Per study (top by GEO GSM)",
            "",
            "| study_id | geo_gsm | with_assay | hub | age | sex | tissue | disease c/ctrl | cancer c/ctrl |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
        ]
    )
    ranked = sorted(per_study, key=lambda r: -int(r["n_geo_gsm"]))  # type: ignore[arg-type]
    for rec in ranked[:60]:
        traits = rec["traits"]  # type: ignore[assignment]
        assert isinstance(traits, dict)
        def _n(pid: str) -> int:
            t = traits.get(pid) or {}
            return int(t.get("n_gsm") or 0)

        def _cc(pid: str) -> str:
            t = traits.get(pid) or {}
            return f"{t.get('n_cases', 0)}/{t.get('n_controls', 0)}"

        lines.append(
            f"| `{rec['study_id']}` | {rec['n_geo_gsm']} | {rec['n_with_ewas_assay']} | "
            f"{rec['n_with_hub_membership']} | {_n('age')} | {_n('sex')} | {_n('tissue')} | "
            f"{_cc('disease')} | {_cc('cancer')} |"
        )
    lines.extend(["", "## Notes", ""])
    lines.extend(f"- {n}" for n in summary["notes"])
    lines.append("")
    (report_dir / "eligibility_by_study.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {report_dir / 'eligibility_by_study.json'}")
    print(f"wrote {report_dir / 'eligibility_by_study.md'}")
    print(
        f"geo_gsm={summary['n_geo_observed_gsm']} "
        f"with_assay={summary['n_geo_gsm_with_assay']} "
        f"hub_overlap={summary['n_geo_gsm_hub_overlap']}"
    )


if __name__ == "__main__":
    main()
