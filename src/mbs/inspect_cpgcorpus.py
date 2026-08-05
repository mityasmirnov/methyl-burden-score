"""CpGCorpus GSE/GPL scientific inspection (sanitized aggregates only)."""

from __future__ import annotations

import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TypedDict

import numpy as np
import pyarrow as pa


class CpgCorpusGplReport(TypedDict):
    source_id: str
    gse: str
    gpl: str
    generated_at: str
    layout: dict[str, Any]
    sample_alignment: dict[str, Any]
    value_qc: dict[str, Any]
    metadata_counts: dict[str, Any]
    warnings: list[str]


_SAFE_METADATA_FIELDS = (
    "platform_id",
    "organism_ch1",
    "molecule_ch1",
    "type",
    "Sex:ch1",
    "sex:ch1",
    "gender:ch1",
    "age (years):ch1",
    "age:ch1",
)


def _open_arrow(path: Path) -> pa.RecordBatchFileReader:
    source = pa.memory_map(str(path), "r")
    return pa.ipc.open_file(source)


def _string_ids(column: Any) -> list[str]:
    values = column.to_pylist()
    return ["" if value is None else str(value) for value in values]


def _value_counts(column: Any, *, max_levels: int = 50) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in column.to_pylist():
        key = "__NULL__" if value is None else str(value)
        counts[key] = counts.get(key, 0) + 1
    if len(counts) > max_levels:
        ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        kept = dict(ranked[:max_levels])
        kept["__truncated_levels__"] = len(counts) - max_levels
        return kept
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def _beta_value_qc(betas_path: Path, *, probe_batch_size: int = 2_048) -> dict[str, Any]:
    reader = _open_arrow(betas_path)
    schema = reader.schema
    if not schema.names or schema.names[0] != "GSM_ID":
        first = schema.names[0] if schema.names else None
        raise ValueError(f"expected first beta column GSM_ID, found {first!r}")

    table = reader.read_all()
    n_samples = table.num_rows
    probe_names = schema.names[1:]
    n_probes = len(probe_names)

    n_values = 0
    n_missing = 0
    n_below_0 = 0
    n_above_1 = 0
    n_finite = 0
    total_sum = 0.0
    total_sq = 0.0
    global_min = math.inf
    global_max = -math.inf
    sample_missing = np.zeros(n_samples, dtype=np.int64)
    sample_sum = np.zeros(n_samples, dtype=np.float64)
    sample_count = np.zeros(n_samples, dtype=np.int64)
    zero_variance_loci = 0
    duplicate_probe_names = n_probes - len(set(probe_names))

    for start in range(0, n_probes, probe_batch_size):
        stop = min(start + probe_batch_size, n_probes)
        columns = probe_names[start:stop]
        batch = table.select(columns)
        array = np.column_stack(
            [batch.column(name).to_numpy(zero_copy_only=False) for name in columns]
        ).astype(np.float64, copy=False)

        missing = ~np.isfinite(array)
        n_missing += int(missing.sum())
        n_values += int(array.size)
        sample_missing += missing.sum(axis=1)

        finite = array[~missing]
        if finite.size:
            n_finite += int(finite.size)
            total_sum += float(finite.sum())
            total_sq += float(np.square(finite).sum())
            global_min = min(global_min, float(finite.min()))
            global_max = max(global_max, float(finite.max()))
            n_below_0 += int((finite < 0.0).sum())
            n_above_1 += int((finite > 1.0).sum())

        valid = np.where(missing, np.nan, array)
        with np.errstate(all="ignore"):
            sample_sum += np.nansum(valid, axis=1)
            sample_count += np.sum(~missing, axis=1)
            col_std = np.nanstd(valid, axis=0)
        zero_variance_loci += int(np.sum(col_std == 0.0))

    mean = total_sum / n_finite if n_finite else None
    variance = (total_sq / n_finite - mean * mean) if n_finite and mean is not None else None
    std = math.sqrt(variance) if variance is not None and variance >= 0 else None
    sample_means = np.divide(
        sample_sum,
        sample_count,
        out=np.full(n_samples, np.nan),
        where=sample_count > 0,
    )
    sample_missing_frac = sample_missing / max(n_probes, 1)

    return {
        "n_samples": n_samples,
        "n_probes": n_probes,
        "n_values": n_values,
        "missing_fraction": (n_missing / n_values) if n_values else None,
        "fraction_below_0": (n_below_0 / n_finite) if n_finite else None,
        "fraction_above_1": (n_above_1 / n_finite) if n_finite else None,
        "min": None if global_min is math.inf else global_min,
        "max": None if global_max is -math.inf else global_max,
        "mean": mean,
        "std": std,
        "zero_variance_loci": zero_variance_loci,
        "duplicate_probe_names": duplicate_probe_names,
        "sample_missing_fraction": {
            "min": float(np.min(sample_missing_frac)) if n_samples else None,
            "median": float(np.median(sample_missing_frac)) if n_samples else None,
            "max": float(np.max(sample_missing_frac)) if n_samples else None,
        },
        "sample_mean_beta": {
            "min": float(np.nanmin(sample_means)) if n_samples else None,
            "median": float(np.nanmedian(sample_means)) if n_samples else None,
            "max": float(np.nanmax(sample_means)) if n_samples else None,
        },
    }


def inspect_cpgcorpus_gpl(
    raw_gpl_root: Path,
    *,
    gse: str,
    gpl: str,
) -> CpgCorpusGplReport:
    """Inspect one CpGCorpus ``GSE/GPL`` directory and return sanitized aggregates."""
    raw_gpl_root = raw_gpl_root.absolute()
    generated_at = datetime.now(UTC).isoformat()
    warnings: list[str] = []
    source_id = f"{gse}_{gpl}"

    betas_path = raw_gpl_root / "betas" / "QCDPB.arrow"
    metadata_path = raw_gpl_root / "metadata" / "metadata.arrow"
    layout = {
        "raw_gpl_root": str(raw_gpl_root),
        "exists": raw_gpl_root.is_dir(),
        "betas_path": str(betas_path),
        "betas_exists": betas_path.is_file(),
        "betas_bytes": betas_path.stat().st_size if betas_path.is_file() else None,
        "metadata_path": str(metadata_path),
        "metadata_exists": metadata_path.is_file(),
        "metadata_bytes": metadata_path.stat().st_size if metadata_path.is_file() else None,
        "orientation": "samples_as_rows_probes_as_columns",
        "beta_id_column": "GSM_ID",
        "processing_level": "QCDPB",
    }

    if not betas_path.is_file() or not metadata_path.is_file():
        warnings.append("missing betas and/or metadata Arrow files")
        return {
            "source_id": source_id,
            "gse": gse,
            "gpl": gpl,
            "generated_at": generated_at,
            "layout": layout,
            "sample_alignment": {},
            "value_qc": {},
            "metadata_counts": {},
            "warnings": warnings,
        }

    beta_reader = _open_arrow(betas_path)
    beta_ids = _string_ids(beta_reader.read_all().column("GSM_ID"))
    meta_table = _open_arrow(metadata_path).read_all()

    if "GSM_ID" not in meta_table.column_names:
        raise ValueError("metadata.arrow is missing GSM_ID")
    meta_ids = _string_ids(meta_table.column("GSM_ID"))

    beta_set = set(beta_ids)
    meta_set = set(meta_ids)
    alignment = {
        "n_beta_samples": len(beta_ids),
        "n_metadata_samples": len(meta_ids),
        "n_unique_beta_ids": len(beta_set),
        "n_unique_metadata_ids": len(meta_set),
        "n_duplicate_beta_ids": len(beta_ids) - len(beta_set),
        "n_duplicate_metadata_ids": len(meta_ids) - len(meta_set),
        "intersection_size": len(beta_set & meta_set),
        "beta_only_count": len(beta_set - meta_set),
        "metadata_only_count": len(meta_set - beta_set),
        "perfect_alignment": beta_set == meta_set and len(beta_ids) == len(beta_set),
    }
    if not alignment["perfect_alignment"]:
        warnings.append("beta and metadata sample IDs are not perfectly aligned")

    metadata_counts: dict[str, Any] = {
        "n_columns": meta_table.num_columns,
        "column_names": meta_table.column_names,
        "fields": {},
    }
    for field in _SAFE_METADATA_FIELDS:
        if field in meta_table.column_names:
            metadata_counts["fields"][field] = _value_counts(meta_table.column(field))

    # Age numeric conversion summary without storing raw ages per sample.
    age_field = next(
        (name for name in ("age (years):ch1", "age:ch1") if name in meta_table.column_names),
        None,
    )
    if age_field is not None:
        raw_ages = meta_table.column(age_field).to_pylist()
        parsed: list[float] = []
        n_missing_age = 0
        n_unparsed = 0
        for value in raw_ages:
            if value is None or str(value).strip() == "":
                n_missing_age += 1
                continue
            try:
                parsed.append(float(str(value).strip()))
            except ValueError:
                n_unparsed += 1
        metadata_counts["age_numeric"] = {
            "field": age_field,
            "n_total": len(raw_ages),
            "n_missing": n_missing_age,
            "n_unparsed": n_unparsed,
            "n_numeric": len(parsed),
            "min": min(parsed) if parsed else None,
            "median": float(np.median(parsed)) if parsed else None,
            "max": max(parsed) if parsed else None,
        }

    if "platform_id" in meta_table.column_names:
        platforms = {str(v) for v in meta_table.column("platform_id").to_pylist() if v is not None}
        if platforms and platforms != {gpl}:
            warnings.append(
                f"metadata platform_id values {sorted(platforms)} differ from path GPL {gpl}"
            )

    value_qc = _beta_value_qc(betas_path)
    if value_qc["fraction_below_0"] and value_qc["fraction_below_0"] > 0:
        warnings.append("observed beta values < 0")
    if value_qc["fraction_above_1"] and value_qc["fraction_above_1"] > 0:
        warnings.append("observed beta values > 1")

    return {
        "source_id": source_id,
        "gse": gse,
        "gpl": gpl,
        "generated_at": generated_at,
        "layout": layout,
        "sample_alignment": alignment,
        "value_qc": value_qc,
        "metadata_counts": metadata_counts,
        "warnings": warnings,
    }


def _markdown(report: CpgCorpusGplReport) -> str:
    lines = [
        f"# Source inspection: `{report['source_id']}`",
        "",
        f"- Generated at: `{report['generated_at']}`",
        f"- GSE: `{report['gse']}`",
        f"- GPL: `{report['gpl']}`",
        "",
        "## Layout",
        "",
        "```json",
        json.dumps(report["layout"], indent=2, sort_keys=True),
        "```",
        "",
        "## Sample alignment",
        "",
        "```json",
        json.dumps(report["sample_alignment"], indent=2, sort_keys=True),
        "```",
        "",
        "## Beta value QC",
        "",
        "```json",
        json.dumps(report["value_qc"], indent=2, sort_keys=True),
        "```",
        "",
        "## Metadata counts (sanitized fields)",
        "",
        "```json",
        json.dumps(report["metadata_counts"], indent=2, sort_keys=True),
        "```",
        "",
        "## Warnings",
        "",
    ]
    if report["warnings"]:
        lines.extend(f"- {warning}" for warning in report["warnings"])
    else:
        lines.append("_None._")
    lines.append("")
    return "\n".join(lines)


def write_cpgcorpus_report(report: CpgCorpusGplReport, report_dir: Path) -> Path:
    """Write sanitized JSON/Markdown inspection outputs."""
    report_dir = report_dir.absolute()
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (report_dir / "summary.md").write_text(_markdown(report), encoding="utf-8")
    (report_dir / "warnings.json").write_text(
        json.dumps(report["warnings"], indent=2) + "\n",
        encoding="utf-8",
    )
    (report_dir / "schema.json").write_text(
        json.dumps(
            {
                "orientation": report["layout"].get("orientation"),
                "beta_id_column": report["layout"].get("beta_id_column"),
                "metadata_columns": report["metadata_counts"].get("column_names", []),
                "n_probes": report["value_qc"].get("n_probes"),
                "n_samples": report["value_qc"].get("n_samples"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return report_dir
