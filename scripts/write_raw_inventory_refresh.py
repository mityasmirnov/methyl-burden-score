#!/usr/bin/env python3
"""Refresh shallow raw_inventory sizes for Hub packs, Atlas, manifests, CpGCorpus."""

from __future__ import annotations

import json
import os
import subprocess
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import pyarrow.parquet as pq

HUB_PACKS = (
    ("age", "age_methylation_v1.zip", 11.73),
    ("tissue", "tissue_methylation_v1.zip", 7.7),
    ("sex", "sex_methylation_v1.zip", 4.33),
    ("blood", "blood_methylation_v1.zip", 4.86),
    ("brain", "brain_methylation_v1.zip", 2.77),
    ("bmi", "bmi_methylation_v1.zip", 3.06),
    ("ancestry", "ancestry_category_methylation_v1.zip", 1.96),
    ("cancer", "cancer_methylation_v1.zip", 16.07),
    ("disease", "disease_methylation_v1.zip", 20.11),
)

STAGE0_CPGCORPUS = (
    ("GSE116992", "GPL13534"),
    ("GSE116992", "GPL21145"),
    ("GSE125367", "GPL21145"),
    ("GSE35069", "GPL13534"),
)

SAMPLE_INFO_ZIPS = (
    ("age", "sample_age_methylation_v1.zip"),
    ("tissue", "sample_tissue_methylation_v1.zip"),
    ("sex", "sample_sex_methylation_v1.zip"),
    ("blood", "sample_blood_methylation_v1.zip"),
    ("brain", "sample_brain_methylation_v1.zip"),
    ("bmi", "sample_bmi_methylation_v1.zip"),
    ("ancestry", "sample_ancestry_category_methylation_v1.zip"),
    ("cancer", "sample_cancer_methylation_v1.zip"),
    ("disease", "sample_disease_methylation_v1.zip"),
)


def _sample_zip_sizes(hub_dl: Path) -> list[dict]:
    rows: list[dict] = []
    for family, fname in SAMPLE_INFO_ZIPS:
        path = hub_dl / fname
        rows.append(
            {
                "family": family,
                "filename": fname,
                "bytes": int(path.stat().st_size) if path.is_file() else None,
                "exists": path.is_file(),
            }
        )
    return rows


def _sample_info_stats(data_root: Path) -> list[dict]:
    """Unique GSM / study counts from canonical sample-info Parquet."""
    ph_dir = data_root / "canonical" / "phenotypes"
    rows: list[dict] = []
    for family, _fname in SAMPLE_INFO_ZIPS:
        path = ph_dir / f"{family}_sample_info.parquet"
        if not path.is_file():
            rows.append({"family": family, "exists": False})
            continue
        df = pq.read_table(path).to_pandas()
        sid = "sample_id" if "sample_id" in df.columns else None
        study = "study_id" if "study_id" in df.columns else None
        rows.append(
            {
                "family": family,
                "exists": True,
                "n_rows": len(df),
                "n_unique_sample_id": int(df[sid].nunique()) if sid else None,
                "n_unique_study_id": int(df[study].nunique()) if study else None,
                "path": str(path),
            }
        )
    return rows


def _bytes(path: Path) -> int | None:
    """File size, or directory total via `du -sb` (avoids slow Python walks)."""
    if not path.is_file() and not path.is_dir():
        return None
    if path.is_file():
        return int(path.stat().st_size)
    try:
        out = subprocess.check_output(  # noqa: S603
            ["/usr/bin/du", "-sb", str(path)],
            text=True,
            stderr=subprocess.DEVNULL,
        )
        return int(out.split()[0])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError):
        total = 0
        for root, _dirs, files in os.walk(path):
            for name in files:
                total += (Path(root) / name).stat().st_size
        return total

def _zip_status(path: Path) -> dict:
    if not path.is_file():
        parent = path.parent
        stem = path.name
        matches = sorted(parent.glob(stem + "*")) if parent.is_dir() else []
        corrupt = [p for p in matches if "corrupt" in p.name.lower()]
        if corrupt:
            c = corrupt[-1]
            return {
                "path": str(c),
                "bytes": int(c.stat().st_size),
                "exists": True,
                "zip_readable": False,
                "status": "quarantined_corrupt",
            }
        return {
            "path": str(path),
            "bytes": None,
            "exists": False,
            "zip_readable": False,
            "status": "missing",
        }
    size = int(path.stat().st_size)
    try:
        with zipfile.ZipFile(path) as zf:
            n = len(zf.namelist())
        return {
            "path": str(path),
            "bytes": size,
            "exists": True,
            "zip_readable": True,
            "n_members": n,
            "status": "ok",
        }
    except zipfile.BadZipFile:
        return {
            "path": str(path),
            "bytes": size,
            "exists": True,
            "zip_readable": False,
            "status": "bad_zip",
        }


def _dir_size_shallow(path: Path) -> dict:
    if not path.is_dir():
        return {"path": str(path), "exists": False, "bytes": None}
    return {"path": str(path), "exists": True, "bytes": _bytes(path)}


def build_inventory(data_root: Path) -> dict:
    hub_dl = data_root / "raw" / "ewas_datahub" / "download"
    packs = []
    for family, fname, advertised_gb in HUB_PACKS:
        info = _zip_status(hub_dl / fname)
        info["family"] = family
        info["filename"] = fname
        info["advertised_gb"] = advertised_gb
        packs.append(info)

    atlas = data_root / "raw" / "ewas_atlas"
    if atlas.is_dir():
        atlas_files = [
            {"name": p.name, "bytes": int(p.stat().st_size)}
            for p in sorted(atlas.glob("*"))
            if p.is_file()
        ]
    else:
        atlas_files = []

    manifests = data_root / "raw" / "manifests" / "epicv2"
    manifest_files = (
        [
            {"name": p.name, "bytes": int(p.stat().st_size)}
            for p in sorted(manifests.glob("*"))
            if p.is_file()
        ]
        if manifests.is_dir()
        else []
    )

    cpg = []
    for gse, gpl in STAGE0_CPGCORPUS:
        root = data_root / "raw" / "cpgcorpus" / gse / gpl
        cpg.append(
            {
                "gse": gse,
                "gpl": gpl,
                "bytes": _bytes(root) if root.exists() else None,
                "exists": root.exists(),
            }
        )

    trees = {
        "ewas_datahub": _dir_size_shallow(data_root / "raw" / "ewas_datahub"),
        "ewas_atlas": _dir_size_shallow(data_root / "raw" / "ewas_atlas"),
        "cpgcorpus": _dir_size_shallow(data_root / "raw" / "cpgcorpus"),
        "manifests": _dir_size_shallow(data_root / "raw" / "manifests"),
        "raw_total": _dir_size_shallow(data_root / "raw"),
    }

    ewas_db = data_root / "raw" / "ewas_datahub" / "EWAS_db"
    n_studies = 0
    if ewas_db.is_dir():
        n_studies = sum(1 for p in ewas_db.iterdir() if p.is_dir())

    sample_info = _sample_info_stats(data_root)
    sample_zips = _sample_zip_sizes(data_root / "raw" / "ewas_datahub" / "download")

    return {
        "generated_at": datetime.now(UTC).isoformat(),
        "hub_profile_packs": packs,
        "hub_sample_info_zips": sample_zips,
        "hub_sample_info_parquet": sample_info,
        "atlas_files": atlas_files,
        "manifest_files": manifest_files,
        "cpgcorpus_stage0": cpg,
        "trees": trees,
        "ewas_db_n_study_dirs": n_studies,
        "ewas_db_remote_study_count": 1989,
        "gmqn": _zip_status(hub_dl / "GMQN.zip"),
        "download_notes": {
            "disease_profile_zip": "complete_2026-08-11 (size+EOCD match remote)",
            "ewas_db": "in_progress via scripts/download_ewas_datahub.sh EWAS_db",
            "host_disk_free_note": "check df -h /data; Hub packs alone ~73 GiB",
        },
    }


def _fmt_size(n: int | None) -> str:
    if n is None:
        return "—"
    if n < 1024**2:
        return f"{n / 1024:.0f} KiB"
    if n < 1024**3:
        return f"{n / (1024**2):.1f} MiB"
    return f"{n / (1024**3):.2f} GiB"


def _gib(n: int | None) -> str:
    return _fmt_size(n)


def markdown(payload: dict) -> str:
    lines = [
        "# Raw data inventory (refreshed)",
        "",
        f"Inspected: **{payload['generated_at']}**. Machine-readable: "
        "[`summary.json`](summary.json).",
        "",
        "Shallow sizes for Hub profile packs, Atlas, EPICv2 manifests, and Stage 0 "
        "CpGCorpus GSEs. Does not recurse into every `EWAS_db` GSM file.",
        "",
        "## Totals under `data/raw/`",
        "",
        "| Tree | Bytes (approx.) |",
        "|------|----------------:|",
    ]
    for key in ("ewas_datahub", "cpgcorpus", "manifests", "ewas_atlas", "raw_total"):
        t = payload["trees"][key]
        lines.append(f"| `{key}/` | {_gib(t.get('bytes'))} |")
    lines.extend(
        [
            "",
            f"EWAS_db study directories: **{payload['ewas_db_n_study_dirs']}**",
            "",
            "## Figures",
            "",
            "![Raw tree sizes](figures/raw_tree_sizes.png)",
            "",
            "![Hub pack sizes](figures/hub_pack_sizes.png)",
            "",
            "![Hub sample counts](figures/hub_sample_counts.png)",
            "",
            "## Hub profile packs (`raw/ewas_datahub/download/`)",
            "",
            "| Family | File | Advertised GB | On-disk | Zip OK | Status |",
            "|--------|------|--------------:|--------:|:------:|--------|",
        ]
    )
    lines.extend(
        [
            (
                f"| {p['family']} | `{p['filename']}` | {p['advertised_gb']} | "
                f"{_gib(p.get('bytes'))} | {p.get('zip_readable')} | `{p['status']}` |"
            )
            for p in payload["hub_profile_packs"]
        ]
    )
    g = payload["gmqn"]
    lines.extend(
        [
            "",
            f"GMQN.zip: {_gib(g.get('bytes'))} status=`{g.get('status')}`",
            "",
            "## Hub sample-info zips (phenotypes, not betas)",
            "",
            "| Family | File | On-disk |",
            "|--------|------|--------:|",
        ]
    )
    lines.extend(
        f"| {z['family']} | `{z['filename']}` | {_gib(z.get('bytes'))} |"
        for z in payload.get("hub_sample_info_zips") or []
    )
    lines.extend(
        [
            "",
            "## Hub sample-info Parquet (`canonical/phenotypes/*_sample_info.parquet`)",
            "",
            "Row counts can exceed unique GSM (duplicate rows in Hub R tables). "
            "Use **unique `sample_id`** as training N.",
            "",
            "| Family | Rows | Unique GSM | Unique studies |",
            "|--------|-----:|-----------:|---------------:|",
        ]
    )
    for row in payload.get("hub_sample_info_parquet") or []:
        if not row.get("exists"):
            lines.append(f"| {row['family']} | — | — | — |")
            continue
        lines.append(
            f"| {row['family']} | {row['n_rows']:,} | "
            f"{row['n_unique_sample_id']:,} | {row['n_unique_study_id']:,} |"
        )
    n_loc = payload["ewas_db_n_study_dirs"]
    n_rem = payload.get("ewas_db_remote_study_count")
    lines.extend(
        [
            "",
            "## EWAS_db All-Data tree (in progress)",
            "",
            f"Local study directories: **{n_loc}** / advertised remote **{n_rem}** "
            f"({100.0 * n_loc / n_rem:.1f}% of study folders if the remote count is stable).",
            "Per-GSM text files under `raw/ewas_datahub/EWAS_db/{GSE}/`. "
            "Resume: `bash scripts/download_ewas_datahub.sh EWAS_db`.",
            "",
            "## EWAS Atlas files",
            "",
            "| File | Bytes |",
            "|------|------:|",
        ]
    )
    lines.extend(f"| `{f['name']}` | {f['bytes']:,} |" for f in payload["atlas_files"])
    lines.extend(
        [
            "",
            "## EPICv2 manifests",
            "",
            "| File | Bytes |",
            "|------|------:|",
        ]
    )
    lines.extend(f"| `{f['name']}` | {f['bytes']:,} |" for f in payload["manifest_files"])
    lines.extend(
        [
            "",
            "## CpGCorpus Stage 0 GSEs",
            "",
            "| GSE | GPL | On-disk |",
            "|-----|-----|--------:|",
        ]
    )
    lines.extend(
        f"| {c['gse']} | {c['gpl']} | {_gib(c.get('bytes'))} |" for c in payload["cpgcorpus_stage0"]
    )
    lines.extend(
        [
            "",
            "## Regenerate",
            "",
            "```bash",
            "uv run python scripts/write_raw_inventory_refresh.py",
            "uv sync --extra analysis  # once, for matplotlib",
            "uv run python scripts/write_pipeline_doc_figures.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    data_root = Path(os.environ.get("MBS_DATA_ROOT", str(project_root / "data")))
    payload = build_inventory(data_root)
    out = project_root / "reports" / "inspection" / "raw_inventory"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out / "summary.md").write_text(markdown(payload), encoding="utf-8")
    print(  # noqa: T201
        json.dumps(
            {"report_dir": str(out), "generated_at": payload["generated_at"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
