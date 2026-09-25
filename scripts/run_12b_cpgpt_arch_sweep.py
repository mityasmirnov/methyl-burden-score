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
  * ranked on ``mbs_enet_nested`` **age MAE** (the product readout, and the
    user-designated primary metric). Nested is fit post-hoc per variant because
    measured restart spread is sd ~0.15 for nested vs ~1.75 for e2e -- ranking
    on e2e would mostly rank noise. Smoke-level, single fold/seed -- not a gate.

Ledger: reports/inspection/stage0_12b_cpgpt_arch_sweep/ledger.json (one row
per variant, appended as each finishes, so a crash keeps prior results).
"""

from __future__ import annotations

import json
import subprocess
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

# Batch size is FIXED across every variant, not `auto`. Peak activation memory
# scales ~ n_edges (57 430) x batch x cpg_hidden, so `auto` would calibrate each
# variant to a different batch and confound capacity with batch size. 128 is the
# largest value that leaves headroom for cpg_hidden=128 on a 49 GB card (the
# first attempt OOM'd at batch 256 with cpg_hidden=128).
MATCHED_BATCH = 128

# (slug, phi/rho width, n_layers, cpg_hidden_dim, dropout)
# `cpg_hidden` is the memory-critical dim (per-edge); phi/rho act on pooled
# gene-level tensors and are comparatively cheap. Capped at 128 for that reason:
# cpg_hidden=256 would need batch 64, and matching that across all arms would
# roughly quadruple wall time. c64-* isolates downstream capacity; c128-* tests
# the 152->cpg_hidden compression that motivated the sweep.
VARIANTS: tuple[tuple[str, int, int, int, float], ...] = (
    ("c64-w64", 64, 2, 64, 0.1),      # in-sweep control (locked N-light shape)
    ("c64-w128", 128, 2, 64, 0.1),
    ("c64-w256", 256, 2, 64, 0.1),
    ("c128-w128", 128, 2, 128, 0.1),  # widens the per-CpG encoder (bottleneck test)
    ("c128-w256", 256, 2, 128, 0.1),
    ("c128-w256-d3", 256, 3, 128, 0.2),
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
    # matched, explicit batch (not "auto") so capacity is not confounded by batch
    cfg.setdefault("training", {})["batch_size"] = MATCHED_BATCH
    GEN_DIR.mkdir(parents=True, exist_ok=True)
    out = GEN_DIR / f"stage0_12b_cpgpt_arch_{slug.replace('-', '_')}.yaml"
    out.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return out


def _fit_nested(run_id: str, cfg_path: Path) -> None:
    """Fit post-hoc nested elastic-net on the saved scores (the product readout).

    Ranking on mbs_e2e alone is unsound here: measured across 6 restarts of the
    same config, e2e age MAE has sd ~1.75 while nested has sd ~0.15 (an 11x
    difference). Single-seed e2e gaps between variants would therefore be mostly
    noise, so the sweep ranks on nested instead.
    """
    subprocess.run(
        [
            "uv", "run", "python", "-u", "scripts/eval_mbs_enet_from_scores.py",
            "--run-id", run_id, "--config", str(cfg_path), "--nested", "--force",
        ],
        cwd=ROOT,
        check=False,  # a failed nested fit must not kill the sweep; recorded as null
    )


def _read_metrics(run_id: str, paths: DataPaths) -> dict[str, Any]:
    """Pull both readouts out of a finished run's metrics.json."""
    mpath = paths.artifact_root / "runs" / run_id / "metrics.json"
    if not mpath.is_file():
        return {"error": f"missing {mpath}"}
    m = json.loads(mpath.read_text(encoding="utf-8"))
    ev = m.get("evaluations") or {}

    def grab(name: str) -> dict[str, Any]:
        g = (ev.get(name) or {}).get("metrics") or {}
        return {
            "tissue_f1": (g.get("tissue") or {}).get("macro_f1"),
            "age_mae": (g.get("age") or {}).get("mae"),
            "sex_auroc": (g.get("sex") or {}).get("auroc"),
        }

    e2e, nested = grab("mbs_e2e"), grab("mbs_enet_nested")
    return {
        "best_epoch": m.get("best_epoch"),
        "nested_tissue_f1": nested["tissue_f1"],
        "nested_age_mae": nested["age_mae"],
        "nested_sex_auroc": nested["sex_auroc"],
        "e2e_tissue_f1": e2e["tissue_f1"],
        "e2e_age_mae": e2e["age_mae"],
        "e2e_sex_auroc": e2e["sex_auroc"],
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
                    f"CpGPT2M on in all arms, batch_size={MATCHED_BATCH} (matched, not auto); "
                    "ranked on mbs_enet_nested age MAE"
                ),
                "ranking_metric": "nested_age_mae (mbs_enet_nested, lower better)",
                "caveat": (
                    "smoke-level, single fold/seed. Ranked on nested, not e2e (e2e "
                    "restart sd ~1.6 vs nested ~0.19). Observed run-to-run spread for a "
                    "FIXED config is ~0.53 MAE wide, so only gaps >~0.5 MAE vs the "
                    "in-sweep control (c64-w64) are credible; smaller orderings are ties."
                ),
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
                    "batch_size": MATCHED_BATCH,
                    "run_id": run_id,
                    "status": "skipped_existing",
                    **_read_metrics(run_id, paths),
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
        if status == "ok":
            _fit_nested(run_id, cfg_path)
        row = {
            "variant": slug,
            "width": width,
            "layers": layers,
            "cpg_hidden": cpg_hidden,
            "dropout": dropout,
            "batch_size": MATCHED_BATCH,
            "run_id": run_id,
            "status": status,
            "wall_seconds": round(time.time() - t0, 1),
            **_read_metrics(run_id, paths),
        }
        _append_ledger(row)
        print(f"[arch-sweep] {slug} -> {row}", flush=True)
    print(f"[arch-sweep] done; ledger={LEDGER}", flush=True)


if __name__ == "__main__":
    main()
