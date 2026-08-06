#!/usr/bin/env python3
"""Compare GSE35069 betas: EWAS Data Hub EWAS_db vs CpGCorpus Arrow.

Writes sanitized aggregates under reports/inspection/GSE35069_hub_vs_cpgcorpus/.
"""

from __future__ import annotations

import json
import math
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pyarrow as pa

from mbs.annotation.manifest import write_json
from mbs.matrix.ewas_db import list_ewas_db_sample_files, read_ewas_db_sample
from mbs.paths import DataPaths


def _open_arrow(path: Path) -> pa.RecordBatchFileReader:
    return pa.ipc.open_file(pa.memory_map(str(path), "r"))


def load_cpgcorpus_betas(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Return (sample_ids, probe_ids, betas[n_samples, n_probes] float64)."""
    table = _open_arrow(path).read_all()
    names = table.schema.names
    if not names or names[0] != "GSM_ID":
        raise ValueError(f"expected GSM_ID first column in {path}")
    sample_ids = ["" if v is None else str(v) for v in table.column("GSM_ID").to_pylist()]
    probe_ids = list(names[1:])
    cols = [
        table.column(name).to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        for name in probe_ids
    ]
    betas = np.column_stack(cols)
    return sample_ids, probe_ids, betas


def main() -> None:
    paths = DataPaths.from_environment()
    hub_dir = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db" / "GSE35069"
    cpg_betas = (
        paths.data_root / "raw" / "cpgcorpus" / "GSE35069" / "GPL13534" / "betas" / "gse_betas.arrow"
    )
    report_dir = paths.project_root / "reports" / "inspection" / "GSE35069_hub_vs_cpgcorpus"
    report_dir.mkdir(parents=True, exist_ok=True)

    if not hub_dir.is_dir():
        raise FileNotFoundError(hub_dir)
    if not cpg_betas.is_file():
        raise FileNotFoundError(cpg_betas)

    print("Loading CpGCorpus arrow…")
    cpg_samples, cpg_probes, cpg_mat = load_cpgcorpus_betas(cpg_betas)
    cpg_probe_index = {p: i for i, p in enumerate(cpg_probes)}
    cpg_sample_index = {s: i for i, s in enumerate(cpg_samples)}

    hub_files = list_ewas_db_sample_files(hub_dir)
    hub_samples = [f.sample_id for f in hub_files]

    shared_samples = sorted(set(hub_samples) & set(cpg_samples))
    hub_only = sorted(set(hub_samples) - set(cpg_samples))
    cpg_only = sorted(set(cpg_samples) - set(hub_samples))

    print(f"Hub samples={len(hub_samples)} CpGCorpus={len(cpg_samples)} shared={len(shared_samples)}")

    # Probe vocabulary from first Hub file + union check via first file size
    first = read_ewas_db_sample(hub_files[0].path)
    hub_probe_set = {str(p) for p in first.probe_ids}
    # Confirm later samples don't add many new probes (sample a few)
    for sample in hub_files[1:5]:
        table = read_ewas_db_sample(sample.path)
        hub_probe_set.update(str(p) for p in table.probe_ids)

    shared_probes = sorted(hub_probe_set & set(cpg_probes))
    hub_only_probes = sorted(hub_probe_set - set(cpg_probes))
    cpg_only_probes = sorted(set(cpg_probes) - hub_probe_set)

    # Pairwise compare on shared sample × shared probe grid (stream Hub samples)
    n_cells = 0
    n_both_finite = 0
    n_hub_missing = 0
    n_cpg_missing = 0
    n_nan_asym = 0
    n_exact_f32 = 0
    n_within_1e4 = 0
    n_within_1e3 = 0
    n_within_1e2 = 0
    abs_sum = 0.0
    abs_sq = 0.0
    max_abs = 0.0
    max_abs_sample = ""
    max_abs_probe = ""
    corr_x_sum = 0.0
    corr_y_sum = 0.0
    corr_xx = 0.0
    corr_yy = 0.0
    corr_xy = 0.0

    shared_probe_arr = np.asarray(shared_probes, dtype=object)
    cpg_cols = np.asarray([cpg_probe_index[p] for p in shared_probes], dtype=np.int64)

    per_sample: list[dict[str, float | int | str]] = []

    for hub_file in hub_files:
        sid = hub_file.sample_id
        if sid not in cpg_sample_index:
            continue
        table = read_ewas_db_sample(hub_file.path, sample_id=sid)
        hub_map = {
            str(pid): float(val) if np.isfinite(val) else float("nan")
            for pid, val in zip(table.probe_ids, table.betas, strict=True)
        }
        hub_vals = np.asarray(
            [hub_map.get(str(p), float("nan")) for p in shared_probe_arr],
            dtype=np.float64,
        )
        cpg_vals = cpg_mat[cpg_sample_index[sid], cpg_cols]

        both = np.isfinite(hub_vals) & np.isfinite(cpg_vals)
        hub_miss = ~np.isfinite(hub_vals) & np.isfinite(cpg_vals)
        cpg_miss = np.isfinite(hub_vals) & ~np.isfinite(cpg_vals)
        both_miss = ~np.isfinite(hub_vals) & ~np.isfinite(cpg_vals)

        n = int(hub_vals.size)
        n_cells += n
        n_both_finite += int(both.sum())
        n_hub_missing += int((~np.isfinite(hub_vals)).sum())
        n_cpg_missing += int((~np.isfinite(cpg_vals)).sum())
        n_nan_asym += int(hub_miss.sum() + cpg_miss.sum())

        local_max = 0.0
        if both.any():
            diff = np.abs(hub_vals[both] - cpg_vals[both])
            hub_f = hub_vals[both]
            cpg_f = cpg_vals[both]
            n_exact_f32 += int((diff == 0.0).sum())
            n_within_1e4 += int((diff <= 1e-4).sum())
            n_within_1e3 += int((diff <= 1e-3).sum())
            n_within_1e2 += int((diff <= 1e-2).sum())
            abs_sum += float(diff.sum())
            abs_sq += float(np.square(diff).sum())
            local_max = float(diff.max())
            if local_max > max_abs:
                max_abs = local_max
                max_abs_sample = sid
                idx = int(np.argmax(diff))
                max_abs_probe = str(shared_probe_arr[both][idx])
            corr_x_sum += float(hub_f.sum())
            corr_y_sum += float(cpg_f.sum())
            corr_xx += float(np.square(hub_f).sum())
            corr_yy += float(np.square(cpg_f).sum())
            corr_xy += float((hub_f * cpg_f).sum())

        per_sample.append(
            {
                "sample_id": sid,
                "n_shared_probes": n,
                "n_both_finite": int(both.sum()),
                "n_hub_only_missing": int(hub_miss.sum()),
                "n_cpg_only_missing": int(cpg_miss.sum()),
                "n_both_missing": int(both_miss.sum()),
                "mae": float(diff.mean()) if both.any() else float("nan"),
                "max_abs_diff": local_max if both.any() else float("nan"),
            }
        )

    n_bf = max(n_both_finite, 1)
    mae = abs_sum / n_bf
    rmse = math.sqrt(abs_sq / n_bf)
    # Pearson on pooled finite pairs
    if n_both_finite > 1:
        mx = corr_x_sum / n_both_finite
        my = corr_y_sum / n_both_finite
        cov = corr_xy / n_both_finite - mx * my
        vx = corr_xx / n_both_finite - mx * mx
        vy = corr_yy / n_both_finite - my * my
        pearson = float(cov / math.sqrt(vx * vy)) if vx > 0 and vy > 0 else float("nan")
    else:
        pearson = float("nan")

    summary = {
        "generated_at": datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "study_id": "GSE35069",
        "hub_path": str(hub_dir),
        "cpgcorpus_betas_path": str(cpg_betas),
        "cpgcorpus_platform": "GPL13534",
        "samples": {
            "hub_n": len(hub_samples),
            "cpgcorpus_n": len(cpg_samples),
            "shared_n": len(shared_samples),
            "hub_only": hub_only,
            "cpgcorpus_only": cpg_only,
            "shared_ids_head": shared_samples[:10],
        },
        "probes": {
            "hub_n_observed_from_sample_scan": len(hub_probe_set),
            "cpgcorpus_n": len(cpg_probes),
            "shared_n": len(shared_probes),
            "hub_only_n": len(hub_only_probes),
            "cpgcorpus_only_n": len(cpg_only_probes),
            "hub_only_head": hub_only_probes[:20],
            "cpgcorpus_only_head": cpg_only_probes[:20],
        },
        "beta_comparison": {
            "n_shared_sample_probe_cells": n_cells,
            "n_both_finite": n_both_finite,
            "n_hub_missing_cells": n_hub_missing,
            "n_cpgcorpus_missing_cells": n_cpg_missing,
            "n_missingness_asymmetric": n_nan_asym,
            "mae": mae,
            "rmse": rmse,
            "pearson_r": pearson,
            "frac_exact": n_exact_f32 / n_bf,
            "frac_abs_diff_le_1e-4": n_within_1e4 / n_bf,
            "frac_abs_diff_le_1e-3": n_within_1e3 / n_bf,
            "frac_abs_diff_le_1e-2": n_within_1e2 / n_bf,
            "max_abs_diff": max_abs,
            "max_abs_diff_sample": max_abs_sample,
            "max_abs_diff_probe": max_abs_probe,
        },
        "notes": [
            "Hub All Data files are GMQN-normalized probe\\tbeta text (ADR 0002).",
            "CpGCorpus Arrow is GEO-derived betas (gse_betas.arrow) for GPL13534.",
            "Exact equality is not expected when processing pipelines differ.",
            "Probe scan used all Hub samples for sample IDs; probe union from first 5 Hub files plus first-file base set.",
        ],
        "per_sample_mae_head": sorted(
            per_sample, key=lambda r: float("-inf") if math.isnan(float(r["mae"])) else float(r["mae"])
        )[:5],
        "per_sample_mae_tail": sorted(
            per_sample,
            key=lambda r: float("inf") if math.isnan(float(r["mae"])) else float(r["mae"]),
            reverse=True,
        )[:5],
    }

    # Full probe union from all Hub samples for accurate probe counts
    print("Scanning all Hub samples for full probe union…")
    full_hub_probes: set[str] = set()
    for hub_file in hub_files:
        table = read_ewas_db_sample(hub_file.path)
        full_hub_probes.update(str(p) for p in table.probe_ids)
    summary["probes"]["hub_n"] = len(full_hub_probes)
    summary["probes"]["shared_n"] = len(full_hub_probes & set(cpg_probes))
    summary["probes"]["hub_only_n"] = len(full_hub_probes - set(cpg_probes))
    summary["probes"]["cpgcorpus_only_n"] = len(set(cpg_probes) - full_hub_probes)
    summary["probes"]["hub_only_head"] = sorted(full_hub_probes - set(cpg_probes))[:20]
    summary["probes"]["cpgcorpus_only_head"] = sorted(set(cpg_probes) - full_hub_probes)[:20]
    del summary["probes"]["hub_n_observed_from_sample_scan"]

    write_json(report_dir / "summary.json", summary)

    bc = summary["beta_comparison"]
    pr = summary["probes"]
    sm = summary["samples"]
    md = [
        "# GSE35069: EWAS Data Hub vs CpGCorpus",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        "## Sources",
        "",
        f"- Hub: `{hub_dir}` ({sm['hub_n']} GSM*.txt)",
        f"- CpGCorpus: `{cpg_betas}` (GPL13534, {sm['cpgcorpus_n']} samples × {pr['cpgcorpus_n']} probes)",
        "",
        "## Sample overlap",
        "",
        f"- Shared: **{sm['shared_n']}** / Hub {sm['hub_n']} / CpGCorpus {sm['cpgcorpus_n']}",
        f"- Hub-only: {sm['hub_only'] or 'none'}",
        f"- CpGCorpus-only: {sm['cpgcorpus_only'] or 'none'}",
        "",
        "## Probe overlap",
        "",
        f"- Shared: **{pr['shared_n']}**",
        f"- Hub-only: {pr['hub_only_n']} (e.g. {pr['hub_only_head'][:5]})",
        f"- CpGCorpus-only: {pr['cpgcorpus_only_n']} (e.g. {pr['cpgcorpus_only_head'][:5]})",
        "",
        "## Beta values (shared sample × shared probe, both finite)",
        "",
        f"- Cells compared (both finite): **{bc['n_both_finite']:,}**",
        f"- Pearson r: **{bc['pearson_r']:.6f}**",
        f"- MAE: **{bc['mae']:.6g}**",
        f"- RMSE: **{bc['rmse']:.6g}**",
        f"- Exact match fraction: {bc['frac_exact']:.4%}",
        f"- |Δ| ≤ 1e-4: {bc['frac_abs_diff_le_1e-4']:.4%}",
        f"- |Δ| ≤ 1e-3: {bc['frac_abs_diff_le_1e-3']:.4%}",
        f"- |Δ| ≤ 1e-2: {bc['frac_abs_diff_le_1e-2']:.4%}",
        f"- Max |Δ|: {bc['max_abs_diff']:.6g} at sample `{bc['max_abs_diff_sample']}` probe `{bc['max_abs_diff_probe']}`",
        f"- Missingness asymmetric cells: {bc['n_missingness_asymmetric']:,}",
        "",
        "## Interpretation",
        "",
        "Hub All Data is GMQN-normalized; CpGCorpus stores GEO-derived betas.",
        "High correlation with non-zero MAE means the same study/samples/probes",
        "are aligned but values are not bit-identical — expected across pipelines.",
        "Do not treat CpGCorpus as the Stage 0 pilot default (ADR 0002).",
        "",
    ]
    (report_dir / "summary.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps({"report_dir": str(report_dir), "pearson_r": bc["pearson_r"], "mae": bc["mae"]}, indent=2))


if __name__ == "__main__":
    main()
