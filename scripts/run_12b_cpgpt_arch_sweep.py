#!/usr/bin/env python3
"""CpGPT architecture sweep (autoresearch-style: fixed budget, one metric, ledger).

Motivation: with CpGPT2M enabled the flat_region per-CpG input width jumps
from **24 -> 152** dims (24 annotation + 128 CpGPT), but phi/rho stayed at
**64**. That makes the encoder a bottleneck that compresses 152 -> 64, so the
65k G2 probe may be *understating* CpGPT by starving it of capacity. This
sweep widens/deepens phi/rho (and cpg_hidden) at fixed budget to find out.

Protocol (every variant identical, so results are comparable):
  * fold 0 of hub-nine-pack-3fold-v1, 65,536-locus prefix, 16 epochs
  * CpGPT2M static embeddings ON in every arm (this sweep is about capacity,
    not about re-litigating whether CpGPT helps -- that was the G2 probe)
  * ranked on in-loop ``mbs_e2e`` **age MAE** (user-designated primary metric).
    ``mbs_enet_nested`` is NOT computed inline (too slow for a sweep); run
    ``eval_mbs_enet_from_scores.py --nested`` on the winner(s) afterwards
    before promoting anything. Smoke-level ranking only -- not a gate.

Ledger: reports/inspection/stage0_12b_cpgpt_arch_sweep/ledger.json (one row
per variant, appended as each finishes, so a crash keeps prior results).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from run_7g_gene_only_probe import train_flat_region_arm  # noqa: E402

from mbs.paths import DataPaths  # noqa: E402

BASE_CONFIG = ROOT / "configs" / "experiment" / "stage0_12b_cpgpt_nlight_smoke.yaml"
GEN_DIR = ROOT / "configs" / "experiment" / "_generated"
REPORT = ROOT / "reports" / "inspection" / "stage0_12b_cpgpt_arch_sweep"
LEDGER = REPORT / "ledger.json"

# (slug, phi/rho width, n_layers, cpg_hidden_dim, dropout)
# base = current locked N-light shape, carried as the in-sweep control.
VARIANTS: tuple[tuple[str, int, int, int, float], ...] = (
    ("base-64", 64, 2, 64, 0.1),
    ("w128", 128, 2, 128, 0.1),
    ("w256", 256, 2, 128, 0.1),
    ("w128-d3", 128, 3, 128, 0.1),
    ("w256-d3-drop2", 256, 3, 256, 0.2),
    ("w512-drop2", 512, 2, 256, 0.2),
)


def _write_variant_config(slug: str, width: int, layers: int, cpg_hidden: int, dropout: float) -> Path:
    cfg: dict[str, Any] = yaml.safe_load(BASE_CONFIG.read_text(encoding="utf-8"))
    cfg["experiment"]["name"] = f"stage0_12b_cpgpt_arch_{slug.replace('-', '_')}"
    model = cfg["model"]
    model["phi_hidden_dimension"] = width
    model["rho_hidden_dimension"] = width
    model["phi_layers"] = layers
    model["rho_layers"] = layers
    model["encoder"]["cpg_hidden_dim"] = cpg_hidden
    model["encoder"]["dropout"] = dropout
    GEN_DIR.mkdir(parents=True, exist_ok=True)
    out = GEN_DIR / f"stage0_12b_cpgpt_arch_{slug.replace('-', '_')}.yaml"
    out.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return out


def _read_e2e(run_id: str, paths: DataPaths) -> dict[str, Any]:
    """Pull mbs_e2e tissue/age/sex out of a finished run's metrics.json."""
    mpath = paths.artifact_root / "runs" / run_id / "metrics.json"
    if not mpath.is_file():
        return {"error": f"missing {mpath}"}
    m = json.loads(mpath.read_text(encoding="utf-8"))
    e2e = ((m.get("evaluations") or {}).get("mbs_e2e") or {}).get("metrics") or {}
    return {
        "best_epoch": m.get("best_epoch"),
        "tissue_f1": (e2e.get("tissue") or {}).get("macro_f1"),
        "age_mae": (e2e.get("age") or {}).get("mae"),
        "sex_auroc": (e2e.get("sex") or {}).get("auroc"),
    }


def _append_ledger(row: dict[str, Any]) -> None:
    REPORT.mkdir(parents=True, exist_ok=True)
    rows = []
    if LEDGER.is_file():
        rows = json.loads(LEDGER.read_text(encoding="utf-8")).get("rows", [])
    rows = [r for r in rows if r.get("variant") != row["variant"]] + [row]
    LEDGER.write_text(
        json.dumps(
            {
                "sweep": "stage0_12b_cpgpt_arch_sweep",
                "protocol": (
                    "fold 0 hub-nine-pack-3fold-v1, 65536-locus prefix, 16 epochs, "
                    "CpGPT2M on in all arms; ranked on mbs_e2e age MAE (primary)"
                ),
                "ranking_metric": "age_mae (mbs_e2e, lower better)",
                "caveat": "smoke-level, single fold/seed; nested enet not computed inline",
                "rows": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    paths = DataPaths.from_environment()
    REPORT.mkdir(parents=True, exist_ok=True)
    for slug, width, layers, cpg_hidden, dropout in VARIANTS:
        run_prefix = f"stage0-12b-cpgpt-arch-{slug}"
        run_id = f"{run_prefix}-f0"
        if (paths.artifact_root / "runs" / run_id / "metrics.json").is_file():
            print(f"[arch-sweep] skip {slug} (metrics.json exists)", flush=True)
            _append_ledger(
                {
                    "variant": slug,
                    "width": width,
                    "layers": layers,
                    "cpg_hidden": cpg_hidden,
                    "dropout": dropout,
                    "run_id": run_id,
                    "status": "skipped_existing",
                    **_read_e2e(run_id, paths),
                }
            )
            continue
        cfg_path = _write_variant_config(slug, width, layers, cpg_hidden, dropout)
        print(
            f"[arch-sweep] === {slug}: phi/rho={width}x{layers} cpg_hidden={cpg_hidden} "
            f"dropout={dropout} cfg={cfg_path.name} ===",
            flush=True,
        )
        t0 = time.time()
        status = "ok"
        try:
            train_flat_region_arm(
                paths=paths,
                config_path=cfg_path,
                run_prefix=run_prefix,
                device="cuda",
                report_dir=REPORT,
                fold_filter=0,
            )
        except Exception as exc:  # keep the sweep alive; record the failure
            status = f"failed: {type(exc).__name__}: {exc}"
            print(f"[arch-sweep] {slug} FAILED: {status}", flush=True)
        row = {
            "variant": slug,
            "width": width,
            "layers": layers,
            "cpg_hidden": cpg_hidden,
            "dropout": dropout,
            "run_id": run_id,
            "status": status,
            "wall_seconds": round(time.time() - t0, 1),
            **_read_e2e(run_id, paths),
        }
        _append_ledger(row)
        print(f"[arch-sweep] {slug} -> {row}", flush=True)
    print(f"[arch-sweep] done; ledger={LEDGER}", flush=True)


if __name__ == "__main__":
    main()
