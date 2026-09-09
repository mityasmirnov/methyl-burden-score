#!/usr/bin/env python3
"""10e correct 1-fold S1→S2→S3→S4 smoke (fold 0 only).

Reuses dense S1 fold-0 checkpoint (scalar mean/mean — imperfect vs preferred
vector-dense S1; documented). Runs:

  S2 — freeze+scalarize max/max, no gene hop
  S3 — freeze RBS stack; train gene_rho (region_hidden)
  S4 — unfreeze FT @ 3e-4

Then post-hoc nested enet on S4 scores and prints a gate table vs native
P2-G fold 0.

Do **not** steal GPU 2 from N-light OOF — default CUDA device is GPU 0.
Leave headroom on GPU 0 (only GPU 2 may be fully occupied): configs cap
``batch_size=32`` / ``gpu_share=2``, and the runner sets
``MBS_CASCADE_GPU_RESERVED_MIB`` (default 16384).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "inspection" / "stage0_7h_nine_pack_smoke"
CFG_S2 = ROOT / "configs/experiment/stage0_7h_nine_pack_staged_s2_scalar_max_max.yaml"
CFG_S3 = ROOT / "configs/experiment/stage0_7h_nine_pack_staged_s3_vector_max_max.yaml"
CFG_S4 = ROOT / "configs/experiment/stage0_7h_nine_pack_staged_s4_ft.yaml"
S1_CKPT = ROOT / "artifacts/runs/stage0-7h-nine-pack-dense-stage1-mean-mean/fold_0/best.pt"
P2G_METRICS = ROOT / "artifacts/runs/stage0-7h-nine-pack-P2-G/fold_0/metrics.json"

STAGES = (
    ("s2", CFG_S2, "stage0-7h-nine-pack-staged-s2-fold0"),
    ("s3", CFG_S3, "stage0-7h-nine-pack-staged-s3-fold0"),
    ("s4", CFG_S4, "stage0-7h-nine-pack-staged-s4-fold0"),
)


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _run(cmd: list[str], *, env: dict[str, str] | None = None) -> None:
    print(f"[staged-s1s4] {_utc()} RUN {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)


def _read_metrics(run_id: str) -> dict[str, Any]:
    path = ROOT / "artifacts" / "runs" / run_id / "fold_0" / "metrics.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _triplet(ev: dict[str, Any] | None) -> tuple[float | None, float | None, float | None]:
    if not isinstance(ev, dict):
        return None, None, None
    metrics = ev.get("metrics") if isinstance(ev.get("metrics"), dict) else ev
    tissue = metrics.get("tissue") if isinstance(metrics, dict) else None
    age = metrics.get("age") if isinstance(metrics, dict) else None
    sex = metrics.get("sex") if isinstance(metrics, dict) else None
    f1 = None
    mae = None
    auroc = None
    if isinstance(tissue, dict):
        f1 = tissue.get("macro_f1")
    if isinstance(age, dict):
        mae = age.get("mae")
    if isinstance(sex, dict):
        auroc = sex.get("auroc")
    return (
        None if f1 is None else float(f1),
        None if mae is None else float(mae),
        None if auroc is None else float(auroc),
    )


def _fmt(v: float | None) -> str:
    return "—" if v is None else f"{v:.3f}"


def _gate_row(name: str, trip: tuple[float | None, float | None, float | None]) -> str:
    return f"| {name} | {_fmt(trip[0])} | {_fmt(trip[1])} | {_fmt(trip[2])} |"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--device",
        default="cuda",
        help="Torch device string passed to mbs train cascade (default: cuda)",
    )
    parser.add_argument(
        "--cuda-visible-devices",
        default="0",
        help="CUDA_VISIBLE_DEVICES for the smoke (default: 0; never steal GPU 2)",
    )
    parser.add_argument(
        "--gpu-reserved-mib",
        type=int,
        default=16384,
        help="MBS_CASCADE_GPU_RESERVED_MIB headroom on shared GPU 0 (default: 16384)",
    )
    parser.add_argument(
        "--skip-train",
        action="store_true",
        help="Skip S2–S4 training; only nested enet + gate table",
    )
    parser.add_argument(
        "--skip-enet",
        action="store_true",
        help="Skip post-hoc nested enet on S4",
    )
    args = parser.parse_args()

    if not S1_CKPT.is_file():
        raise SystemExit(f"missing S1 checkpoint: {S1_CKPT}")

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(args.cuda_visible_devices)
    # N-light owns physical GPU 2 (may fill it). GPU 0 must keep spare VRAM.
    env["MBS_CASCADE_GPU_RESERVED_MIB"] = str(max(0, int(args.gpu_reserved_mib)))
    env.setdefault("MBS_CASCADE_GPU_SHARE", "2")
    print(
        f"[staged-s1s4] {_utc()} CUDA_VISIBLE_DEVICES={env['CUDA_VISIBLE_DEVICES']} "
        f"device={args.device} reserved_mib={env['MBS_CASCADE_GPU_RESERVED_MIB']} "
        f"gpu_share={env.get('MBS_CASCADE_GPU_SHARE')} (S1 reuse: {S1_CKPT})",
        flush=True,
    )

    if not args.skip_train:
        for stage, cfg, run_id in STAGES:
            report = REPORT / f"_staging_staged_{stage}_fold0"
            report.mkdir(parents=True, exist_ok=True)
            cmd = [
                "uv",
                "run",
                "mbs",
                "train",
                "cascade",
                "--config",
                str(cfg),
                "--run-id",
                run_id,
                "--device",
                args.device,
                "--max-folds",
                "1",
                "--skip-if-done",
                "--report-dir",
                str(report),
            ]
            _run(cmd, env=env)

    s4_run = STAGES[-1][2]
    if not args.skip_enet:
        cmd = [
            "uv",
            "run",
            "python",
            str(ROOT / "scripts" / "eval_mbs_enet_from_scores.py"),
            "--config",
            str(CFG_S4),
            "--run-id",
            s4_run,
            "--which",
            "both",
            "--nested",
        ]
        _run(cmd, env=env)

    s4 = _read_metrics(s4_run)
    s4_ev = s4.get("evaluations") if isinstance(s4.get("evaluations"), dict) else {}
    p2g = json.loads(P2G_METRICS.read_text(encoding="utf-8"))
    p2g_ev = p2g.get("evaluations") if isinstance(p2g.get("evaluations"), dict) else {}

    rows = [
        ("S4 mbs_e2e", _triplet(s4_ev.get("mbs_e2e"))),
        ("S4 mbs_linear_probe", _triplet(s4_ev.get("mbs_linear_probe"))),
        ("S4 rbs_linear_probe", _triplet(s4_ev.get("rbs_linear_probe"))),
        ("S4 mbs_enet_nested", _triplet(s4_ev.get("mbs_enet_nested"))),
        ("S4 rbs_enet_nested", _triplet(s4_ev.get("rbs_enet_nested"))),
        ("P2-G f0 mbs_e2e", _triplet(p2g_ev.get("mbs_e2e"))),
        ("P2-G f0 mbs_linear_probe", _triplet(p2g_ev.get("mbs_linear_probe"))),
        ("P2-G f0 rbs_linear_probe", _triplet(p2g_ev.get("rbs_linear_probe"))),
        ("P2-G f0 mbs_enet_nested", _triplet(p2g_ev.get("mbs_enet_nested"))),
    ]

    s4_e2e = _triplet(s4_ev.get("mbs_e2e"))
    p2g_e2e = _triplet(p2g_ev.get("mbs_e2e"))
    e2e_beats = (
        s4_e2e[0] is not None
        and p2g_e2e[0] is not None
        and s4_e2e[0] > p2g_e2e[0]
        and s4_e2e[1] is not None
        and p2g_e2e[1] is not None
        and s4_e2e[1] < p2g_e2e[1]
    )
    # Gate also wants probes/enet match-or-beat; report clearly either way.
    gate = "PASS (promote 3-fold candidate)" if e2e_beats else "FAIL / hold native P2-G"
    out = {
        "updated_utc": _utc(),
        "s1_checkpoint": str(S1_CKPT),
        "s1_note": "scalar mean/mean dense S1 fold0 (imperfect vs preferred vector-dense)",
        "stages": {name: run_id for name, _, run_id in STAGES},
        "gate_e2e_tissue_and_age": gate,
        "rows": [
            {
                "name": name,
                "tissue_f1": trip[0],
                "age_mae": trip[1],
                "sex_auroc": trip[2],
            }
            for name, trip in rows
        ],
    }
    out_path = REPORT / "staged_s1s4_fold0_gate.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    print("\n## Fair staged S1→S2→S3→S4 vs P2-G fold 0\n", flush=True)
    print("| Readout | Tissue F1 | Age MAE | Sex AUROC |", flush=True)
    print("|---------|----------:|--------:|----------:|", flush=True)
    for name, trip in rows:
        print(_gate_row(name, trip), flush=True)
    print(f"\nGate (e2e tissue↑ and age↓ vs P2-G f0): **{gate}**", flush=True)
    print(f"Wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
