#!/usr/bin/env python3
"""Fail-safe weekend supervisor for 7G′ matched 16-ep screen + optional seed-mask.

Owns GPU 0 via flock. Does not duplicate live training. Does not launch Stage B,
P5, or final OOF. Prefer completing fewer defensible experiments.

Modes:
  --status   print queue / running / complete summary and exit
  --dry-run  plan actions without launching GPU or writing reports
  (default)  run the state machine until hard stop

Event log: scratch/logs/7g_prime_weekend_events.jsonl
Lock:      scratch/locks/7g_prime_weekend.lock
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts" / "runs"
REPORT = ROOT / "reports" / "inspection" / "stage0_7g_gene_only_probe"
SEED_REPORT = ROOT / "reports" / "inspection" / "stage0_7g_prime_seed_mask"
SCRATCH = ROOT / "scratch"
LOCK_PATH = SCRATCH / "locks" / "7g_prime_weekend.lock"
EVENT_LOG = SCRATCH / "logs" / "7g_prime_weekend_events.jsonl"
PID_PATH = SCRATCH / "logs" / "7g_prime_weekend.pid"
SUMMARY_PATH = SCRATCH / "logs" / "7g_prime_weekend_summary.json"
BLOCK_PATH = SCRATCH / "SEED_MASK_GPU_BLOCKED.txt"
MAIN_CFG = ROOT / "configs" / "experiment" / "stage0_7g_gene_only_probe.yaml"
RUNNER = ROOT / "scripts" / "run_7g_gene_only_probe.py"
SEED_RUNNER = ROOT / "scripts" / "run_7g_prime_seed_mask.py"
WRITE_REPORT = ROOT / "scripts" / "write_7g_gene_only_probe_report.py"
APPLY_DECISION = ROOT / "scripts" / "apply_7g_16ep_decision.py"
ENET = ROOT / "scripts" / "eval_mbs_enet_from_scores.py"

# Matched 16-ep promotion arms (order after any already-running job finishes).
# User priority prefers vector, but never interrupt an in-flight scalar fold.
GPU_ARMS: list[dict[str, Any]] = [
    {
        "arm_id": "N-cascade-scalar-max-mean",
        "run_id": "stage0-7g-gene-probe-scalar-max-mean-16ep",
        "kind": "cascade",
        "min_ceiling": 15,
    },
    {
        "arm_id": "N-cascade-vector-mean-max",
        "run_id": "stage0-7g-gene-probe-vector-mean-max-16ep",
        "kind": "cascade",
        "min_ceiling": 15,
    },
]
# Already-complete references (validated, not retrained).
DONE_ARMS: list[dict[str, Any]] = [
    {"arm_id": "P2-G", "run_id": "stage0-7g-gene-probe-P2-G-explicit", "kind": "cascade", "min_ceiling": 15},
    {
        "arm_id": "N-light-gene-max",
        "run_prefix": "stage0-7g-gene-probe-light-max",
        "kind": "flat",
        "min_ceiling": 15,
    },
    {
        "arm_id": "N-light-gene-mean",
        "run_prefix": "stage0-7g-gene-probe-light-mean-16ep",
        "kind": "flat",
        "min_ceiling": 15,
    },
    {
        "arm_id": "N-cascade-scalar-mean-max",
        "run_id": "stage0-7g-gene-probe-scalar-mean-max-16ep",
        "kind": "cascade",
        "min_ceiling": 15,
    },
]

POLL_S = 120
GPU_ID = 0


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha_short(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:16]


def emit(event: str, **payload: Any) -> None:
    EVENT_LOG.parent.mkdir(parents=True, exist_ok=True)
    row = {"ts": _utc(), "event": event, **payload}
    with EVENT_LOG.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, default=str) + "\n")
    print(f"[{row['ts']}] {event} {json.dumps(payload, default=str)[:300]}", flush=True)


def git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


@dataclass
class FoldStatus:
    arm_id: str
    fold: int
    path: str
    status: str  # ok | incomplete | missing | running
    matched16: bool = False
    tissue_f1: float | None = None
    age_mae: float | None = None
    sex_auroc: float | None = None
    best_epoch: int | None = None
    max_epochs: int | None = None
    eval_split: str | None = None
    reason: str | None = None


def _cascade_fold_dir(run_id: str, fold: int) -> Path:
    return ART / run_id / f"fold_{fold}"


def _flat_run_dir(prefix: str, fold: int) -> Path:
    return ART / f"{prefix}-f{fold}"


def validate_fold(
    *,
    arm_id: str,
    fold: int,
    kind: str,
    run_id: str | None = None,
    run_prefix: str | None = None,
    min_ceiling: int = 15,
) -> FoldStatus:
    if kind == "cascade":
        assert run_id
        root = _cascade_fold_dir(run_id, fold)
    else:
        assert run_prefix
        root = _flat_run_dir(run_prefix, fold)
    mp = root / "metrics.json"
    if not mp.is_file():
        return FoldStatus(arm_id, fold, str(root), "missing", reason="metrics.json missing")
    try:
        blob = json.loads(mp.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return FoldStatus(arm_id, fold, str(root), "incomplete", reason=f"metrics unreadable: {exc}")
    ev = blob.get("evaluations") or {}
    e2e = ev.get("mbs_e2e") or {}
    metrics = e2e.get("metrics") or {}
    tissue = (metrics.get("tissue") or {}).get("macro_f1")
    age = (metrics.get("age") or {}).get("mae")
    sex = (metrics.get("sex") or {}).get("auroc")
    split = e2e.get("eval_split")
    best = blob.get("best_epoch")
    ckpt = blob.get("checkpoint_selection") or {}
    if best is None and isinstance(ckpt, dict):
        best = ckpt.get("best_epoch")
    max_ep = None
    if isinstance(ckpt, dict) and ckpt.get("max_epochs") is not None:
        max_ep = int(ckpt["max_epochs"])
    elif blob.get("max_epochs") is not None:
        max_ep = int(blob["max_epochs"])
    hist = []
    if isinstance(ckpt, dict):
        hist = ckpt.get("val_history") or []
    if not hist:
        hist = blob.get("history") or []
    reached = None
    epochs = [int(h["epoch"]) for h in hist if isinstance(h, dict) and h.get("epoch") is not None]
    if epochs:
        reached = max(epochs)
    matched = (max_ep is not None and max_ep >= min_ceiling) or (
        reached is not None and reached >= min_ceiling
    ) or (best is not None and int(best) >= min_ceiling)
    scores = root / "scores"
    has_mbs = (scores / "mbs.zarr").exists() or (scores / "mbs.npy").is_file()
    pending = bool(blob.get("cpu_probes_pending"))
    finite_ok = True
    for v in (tissue, age):
        if v is None:
            continue
        if not isinstance(v, (int, float)) or v != v or abs(float(v)) > 1e6:
            finite_ok = False
    reasons: list[str] = []
    if split != "test":
        reasons.append(f"eval_split={split!r}")
    if tissue is None:
        reasons.append("missing tissue")
    if age is None:
        reasons.append("missing age")
    if best is None:
        reasons.append("missing best_epoch")
    if not has_mbs:
        reasons.append("missing mbs scores")
    if not matched:
        reasons.append("not matched-16 ceiling")
    if pending:
        reasons.append("cpu_probes_pending")
    if not finite_ok:
        reasons.append("non_finite_or_implausible_metric")
    ok = not reasons
    return FoldStatus(
        arm_id=arm_id,
        fold=fold,
        path=str(root),
        status="ok" if ok else "incomplete",
        matched16=bool(matched),
        tissue_f1=float(tissue) if tissue is not None else None,
        age_mae=float(age) if age is not None else None,
        sex_auroc=float(sex) if sex is not None else None,
        best_epoch=int(best) if best is not None else None,
        max_epochs=max_ep,
        eval_split=str(split) if split is not None else None,
        reason="; ".join(reasons) if reasons else None,
    )


def discover_status() -> dict[str, Any]:
    arms = DONE_ARMS + GPU_ARMS
    out: dict[str, Any] = {"arms": {}, "missing_gpu": [], "running": detect_gpu_jobs()}
    for arm in arms:
        folds = []
        for i in range(3):
            st = validate_fold(
                arm_id=str(arm["arm_id"]),
                fold=i,
                kind=str(arm["kind"]),
                run_id=arm.get("run_id"),
                run_prefix=arm.get("run_prefix"),
                min_ceiling=int(arm.get("min_ceiling", 15)),
            )
            folds.append(asdict(st))
            if st.status != "ok" and arm in GPU_ARMS:
                out["missing_gpu"].append({"arm_id": arm["arm_id"], "fold": i, "status": st.status})
        n_ok = sum(1 for f in folds if f["status"] == "ok")
        out["arms"][arm["arm_id"]] = {
            "n_ok": n_ok,
            "folds": folds,
            "run_id": arm.get("run_id"),
            "run_prefix": arm.get("run_prefix"),
        }
    return out


def detect_gpu_jobs() -> list[dict[str, Any]]:
    """Detect live 7G′ GPU owners on CUDA device 0 without killing anything."""
    jobs: list[dict[str, Any]] = []
    try:
        ps = subprocess.check_output(["ps", "-ef"], text=True)
    except Exception:
        return jobs
    needles = (
        "run_7g_gene_only_probe",
        "mbs train cascade",
        "run_7g_16ep_promotion",
        "run_7g_prime_seed_mask",
        "supervise_7g_prime_weekend",
    )
    for line in ps.splitlines():
        if any(n in line for n in needles) and "rg" not in line and "grep" not in line:
            parts = line.split()
            pid = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else None
            jobs.append({"pid": pid, "cmd": line[line.find(parts[7]) :] if len(parts) > 7 else line})
    return jobs


def gpu0_busy_for_new_train(threshold_mib: int = 6000) -> bool:
    """True if another approved owner still holds the queue (do not start another train).

    Importantly: treat ``run_7g_16ep_promotion_resume.sh`` as the owner even while
    it is between arms (GPU memory may drop during sync/report). Never start a
    duplicate fold while that bash parent is alive.
    """
    jobs = detect_gpu_jobs()
    my_pid = os.getpid()
    foreign = []
    for j in jobs:
        if j.get("pid") == my_pid:
            continue
        cmd = j.get("cmd") or ""
        if "supervise_7g_prime_weekend" in cmd:
            continue
        if "eval_mbs_enet" in cmd:
            continue
        if any(
            x in cmd
            for x in (
                "mbs train cascade",
                "run_7g_gene_only_probe",
                "run_7g_16ep_promotion",
                "run_7g_prime_seed_mask",
            )
        ):
            foreign.append(j)
    if foreign:
        return True
    try:
        out = subprocess.check_output(
            [
                "nvidia-smi",
                f"--id={GPU_ID}",
                "--query-gpu=memory.used",
                "--format=csv,noheader,nounits",
            ],
            text=True,
        ).strip()
        used = int(float(out.splitlines()[0].strip()))
        return used > threshold_mib
    except Exception:
        return False


def run_cmd(
    cmd: list[str],
    *,
    log_path: Path,
    dry_run: bool,
    env: dict[str, str] | None = None,
) -> int:
    emit("cmd_start", cmd=cmd, log=str(log_path), dry_run=dry_run)
    if dry_run:
        return 0
    log_path.parent.mkdir(parents=True, exist_ok=True)
    merged = os.environ.copy()
    if env:
        merged.update(env)
    merged.setdefault("CUDA_VISIBLE_DEVICES", str(GPU_ID))
    merged.setdefault("PYTHONUNBUFFERED", "1")
    merged.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(f"\n===== {_utc()} CMD {' '.join(cmd)} =====\n")
        fh.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            stdout=fh,
            stderr=subprocess.STDOUT,
            env=merged,
        )
        emit("cmd_pid", pid=proc.pid, cmd=cmd)
        rc = proc.wait()
    emit("cmd_end", cmd=cmd, returncode=rc, log=str(log_path))
    return int(rc)


def sync_per_arm_cascade(arm_id: str, run_id: str) -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_7g_gene_only_probe import _mean_metric, _slim_cascade_fold, write_per_arm  # type: ignore

    folds: list[dict[str, Any]] = []
    for i in range(3):
        mp = _cascade_fold_dir(run_id, i) / "metrics.json"
        if mp.is_file():
            folds.append(_slim_cascade_fold(json.loads(mp.read_text(encoding="utf-8"))))
    write_per_arm(
        REPORT,
        arm_id,
        {
            "arm_id": arm_id,
            "kind": "cascade_train",
            "run_id": run_id,
            "folds": folds,
            "mean_mbs_e2e_tissue_f1": _mean_metric(folds, "mbs_e2e.metrics.tissue.macro_f1"),
            "mean_mbs_linear_probe_tissue_f1": _mean_metric(
                folds, "mbs_linear_probe.metrics.tissue.macro_f1"
            ),
            "mean_mbs_enet_tissue_f1": _mean_metric(folds, "mbs_enet.metrics.tissue.macro_f1"),
            "mean_mbs_enet_nested_tissue_f1": _mean_metric(
                folds, "mbs_enet_nested.metrics.tissue.macro_f1"
            ),
            "mean_rbs_linear_probe_tissue_f1": _mean_metric(
                folds, "rbs_linear_probe.metrics.tissue.macro_f1"
            ),
            "mean_rbs_enet_tissue_f1": _mean_metric(folds, "rbs_enet.metrics.tissue.macro_f1"),
        },
    )


def sync_per_arm_flat(arm_id: str, prefix: str, pool: str) -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_7g_gene_only_probe import _mean_metric, _slim_cascade_fold, write_per_arm  # type: ignore

    folds: list[dict[str, Any]] = []
    for i in range(3):
        mp = _flat_run_dir(prefix, i) / "metrics.json"
        if mp.is_file():
            folds.append(_slim_cascade_fold(json.loads(mp.read_text(encoding="utf-8"))))
    write_per_arm(
        REPORT,
        arm_id,
        {
            "arm_id": arm_id,
            "kind": "flat_region_train",
            "run_prefix": prefix,
            "pool": pool,
            "folds": folds,
            "mean_mbs_e2e_tissue_f1": _mean_metric(folds, "mbs_e2e.metrics.tissue.macro_f1"),
            "mean_mbs_linear_probe_tissue_f1": _mean_metric(
                folds, "mbs_linear_probe.metrics.tissue.macro_f1"
            ),
            "mean_mbs_enet_tissue_f1": _mean_metric(folds, "mbs_enet.metrics.tissue.macro_f1"),
            "mean_mbs_enet_nested_tissue_f1": _mean_metric(
                folds, "mbs_enet_nested.metrics.tissue.macro_f1"
            ),
        },
    )


def refresh_reports(dry_run: bool) -> None:
    if dry_run:
        emit("refresh_reports_skipped", dry_run=True)
        return
    # Keep Tier-1 / promotion arms in sync from artifacts.
    sync_per_arm_flat("N-light-gene-max", "stage0-7g-gene-probe-light-max", "max")
    sync_per_arm_flat("N-light-gene-mean", "stage0-7g-gene-probe-light-mean-16ep", "mean")
    sync_per_arm_cascade("N-cascade-scalar-mean-max", "stage0-7g-gene-probe-scalar-mean-max-16ep")
    sync_per_arm_cascade("N-cascade-scalar-max-mean", "stage0-7g-gene-probe-scalar-max-mean-16ep")
    # Prefer 16-ep vector when present; else leave previous per_arm untouched if empty.
    vec_ok = sum(
        1
        for i in range(3)
        if validate_fold(
            arm_id="N-cascade-vector-mean-max",
            fold=i,
            kind="cascade",
            run_id="stage0-7g-gene-probe-vector-mean-max-16ep",
            min_ceiling=15,
        ).status
        == "ok"
    )
    if vec_ok:
        sync_per_arm_cascade(
            "N-cascade-vector-mean-max", "stage0-7g-gene-probe-vector-mean-max-16ep"
        )
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_cmd(
        ["uv", "run", "python", str(WRITE_REPORT)],
        log_path=SCRATCH / "logs" / f"weekend_write_report_{ts}.log",
        dry_run=False,
    )
    run_cmd(
        ["uv", "run", "python", str(APPLY_DECISION)],
        log_path=SCRATCH / "logs" / f"weekend_apply_decision_{ts}.log",
        dry_run=False,
    )


def train_arm(arm: dict[str, Any], *, dry_run: bool) -> int:
    arm_id = str(arm["arm_id"])
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = SCRATCH / "logs" / f"weekend_train_{arm_id}_{ts}.log"
    cmd = [
        "uv",
        "run",
        "python",
        "-u",
        str(RUNNER),
        "--config",
        str(MAIN_CFG),
        "--device",
        "cuda",
        "--arm",
        arm_id,
    ]
    # Runner uses --skip-if-done inside cascade train.
    return run_cmd(cmd, log_path=log_path, dry_run=dry_run)


def posthoc_enet_for_completed(*, dry_run: bool) -> None:
    """CPU enet after GPU queue; skip if another heavy enet already owns RAM."""
    jobs = detect_gpu_jobs()
    if any("eval_mbs_enet_from_scores" in j["cmd"] for j in jobs):
        emit("posthoc_deferred", reason="enet already running")
        return
    targets = [
        ("--run-id", "stage0-7g-gene-probe-scalar-max-mean-16ep"),
        ("--run-id", "stage0-7g-gene-probe-vector-mean-max-16ep"),
        ("--run-id", "stage0-7g-gene-probe-scalar-mean-max-16ep"),
        ("--run-prefix", "stage0-7g-gene-probe-light-max"),
        ("--run-prefix", "stage0-7g-gene-probe-light-mean-16ep"),
    ]
    for flag, value in targets:
        # Prefer mbs first (faster / safer); nested next. Skip rbs nested overnight
        # unless mbs nested already present (rbs can be multi-hour).
        for nested in (False, True):
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            tag = value.replace("/", "_")
            log_path = SCRATCH / "logs" / f"weekend_enet_{tag}_{'nested' if nested else 'fixed'}_{ts}.log"
            cmd = ["uv", "run", "python", "-u", str(ENET), flag, value, "--which", "mbs"]
            if nested:
                cmd.append("--nested")
            rc = run_cmd(cmd, log_path=log_path, dry_run=dry_run)
            if rc != 0:
                emit("posthoc_failed", target=value, nested=nested, returncode=rc)


def seed_panel_audit_ok() -> tuple[bool, list[str]]:
    problems: list[str] = []
    panel = SEED_REPORT / "seed_panels" / "fold_0" / "seed_panel.json"
    if not panel.is_file():
        return False, ["missing fold_0 seed_panel.json"]
    try:
        blob = json.loads(panel.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return False, [f"seed_panel unreadable: {exc}"]
    if not blob.get("graph_content_hash"):
        problems.append("null graph_content_hash")
    if not blob.get("panel_hash"):
        problems.append("null panel_hash")
    traits = blob.get("traits") or {}
    for trait_name in ("age", "tissue", "sex", "sex_autosome"):
        trait = traits.get(trait_name) or {}
        if not isinstance(trait, dict):
            continue
        n_sex = trait.get("n_sex_chrom_seed_cpgs")
        if n_sex is not None and int(n_sex) != 0:
            problems.append(f"{trait_name}.n_sex_chrom_seed_cpgs={n_sex}")
    age = traits.get("age") or {}
    if isinstance(age, dict) and age.get("ranking_fallback"):
        # Documented fallback is allowed; do not treat as failure.
        emit(
            "seed_age_fallback_noted",
            ranking_fallback=age.get("ranking_fallback"),
            sparsity_ok=age.get("sparsity_ok"),
        )
    g3 = blob.get("g3_matched_random_quality") or {}
    if isinstance(g3, dict):
        emit(
            "seed_g3_quality",
            fraction_exact=g3.get("fraction_exact_cpg_match"),
            abs_err_max=g3.get("cpg_count_abs_err_max"),
            abs_err_p90=g3.get("cpg_count_abs_err_p90"),
            gene_length_bp_used=g3.get("gene_length_bp_used"),
        )
    audit_md = SEED_REPORT / "panel_audit.md"
    if audit_md.is_file() and "ok_for_seed_mask_gpu: false" in audit_md.read_text(
        encoding="utf-8"
    ):
        problems.append("panel_audit.md marks ok_for_seed_mask_gpu: false")
    return (len(problems) == 0), problems


def maybe_unlock_and_run_seed(*, dry_run: bool) -> str:
    decision_path = REPORT / "promotion_decision.json"
    if not decision_path.is_file():
        return "no_decision"
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    if not decision.get("screen_complete"):
        return "screen_incomplete"
    ok, problems = seed_panel_audit_ok()
    if not ok:
        emit("seed_audit_blocked", problems=problems)
        return "seed_audit_blocked"
    next_gate = str(decision.get("next_gate") or "")
    allowed = {
        "age_primary_seed_mask",
        "manual_review",
        "prefer_one_hop_max",
        "retain_pooling_2x2",
    }
    if next_gate not in allowed and next_gate != "typed_rbs_aggregation":
        emit("seed_gate_not_allowed", next_gate=next_gate)
        return f"next_gate={next_gate}"
    if next_gate == "typed_rbs_aggregation":
        # Still allow seed-mask as the planned small experiment; typed-RBS is CPU follow-up.
        emit("seed_note", note="typed_rbs rule fired; still running planned seed-mask grid")
    if BLOCK_PATH.is_file() and not dry_run:
        BLOCK_PATH.unlink()
        emit("seed_block_removed", path=str(BLOCK_PATH), next_gate=next_gate)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = SCRATCH / "logs" / f"weekend_seed_mask_{ts}.log"
    cmd = [
        "uv",
        "run",
        "python",
        "-u",
        str(SEED_RUNNER),
        "--device",
        "cuda",
        "--reuse-panels",
    ]
    rc = run_cmd(cmd, log_path=log_path, dry_run=dry_run)
    return "seed_done" if rc == 0 else f"seed_failed_rc={rc}"


def try_commit_push(message: str, *, dry_run: bool) -> str | None:
    if dry_run:
        return None
    # Only stage intended compact reports + supervisor sources.
    paths = [
        "scripts/supervise_7g_prime_weekend.py",
        "scripts/run_7g_prime_weekend.sh",
        "reports/inspection/stage0_7g_gene_only_probe/analysis.md",
        "reports/inspection/stage0_7g_gene_only_probe/promotion_decision.json",
        "reports/inspection/stage0_7g_gene_only_probe/task_comparison.json",
        "reports/inspection/stage0_7g_gene_only_probe/summary.json",
        "reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json",
    ]
    per_arm = REPORT / "per_arm"
    for name in (
        "N-light-gene-max.json",
        "N-light-gene-mean.json",
        "N-cascade-scalar-mean-max.json",
        "N-cascade-scalar-max-mean.json",
        "N-cascade-vector-mean-max.json",
        "P2-G.json",
    ):
        p = per_arm / name
        # Prefer slim per_arm only (fat history/ROC blobs stay untracked).
        if p.is_file() and p.stat().st_size < 3_000_000:
            paths.append(str(p.relative_to(ROOT)))
        elif p.is_file():
            emit("commit_skip_large_per_arm", path=str(p), bytes=p.stat().st_size)
    for extra in (
        "reports/inspection/stage0_7g_gene_only_probe/weekend_supervisor_summary.md",
        "reports/inspection/stage0_7g_prime_seed_mask/analysis.md",
        "reports/inspection/stage0_7g_prime_seed_mask/summary.json",
        "scripts/apply_7g_16ep_decision.py",
        "tests/unit/test_supervise_7g_prime_weekend.py",
    ):
        if (ROOT / extra).exists():
            paths.append(extra)
    existing = [p for p in paths if (ROOT / p).exists()]
    if not existing:
        return None
    subprocess.run(["git", "add", "--", *existing], cwd=ROOT, check=False)
    # Skip empty commit
    staged = subprocess.run(
        ["git", "diff", "--cached", "--quiet"], cwd=ROOT, check=False
    ).returncode
    if staged == 0:
        emit("commit_skipped", reason="nothing staged")
        return None
    subprocess.run(["git", "commit", "-m", message], cwd=ROOT, check=False)
    subprocess.run(["git", "fetch", "origin"], cwd=ROOT, check=False)
    push = subprocess.run(["git", "push", "origin", "HEAD"], cwd=ROOT, check=False)
    sha = git_head()
    emit("pushed", sha=sha, push_rc=push.returncode, message=message)
    return sha


def write_summary(payload: dict[str, Any]) -> None:
    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_PATH.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    final_md = REPORT / "weekend_supervisor_summary.md"
    lines = [
        "# 7G′ weekend supervisor summary",
        "",
        f"Generated: `{_utc()}`",
        "",
        "```json",
        json.dumps(payload, indent=2, default=str),
        "```",
        "",
    ]
    final_md.write_text("\n".join(lines), encoding="utf-8")


def acquire_lock() -> Any:
    LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    fh = LOCK_PATH.open("w", encoding="utf-8")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        fh.close()
        raise SystemExit(f"another weekend supervisor holds {LOCK_PATH}") from exc
    fh.write(f"pid={os.getpid()}\nstarted={_utc()}\nhead={git_head()}\n")
    fh.flush()
    PID_PATH.parent.mkdir(parents=True, exist_ok=True)
    PID_PATH.write_text(str(os.getpid()) + "\n", encoding="utf-8")
    return fh


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--once", action="store_true", help="Single pass then exit")
    parser.add_argument("--skip-seed", action="store_true", help="Hard-stop after promotion screen")
    parser.add_argument("--poll-seconds", type=int, default=POLL_S)
    args = parser.parse_args()

    status = discover_status()
    if args.status:
        print(json.dumps(status, indent=2, default=str))
        return

    lock_fh = acquire_lock()
    emit(
        "supervisor_start",
        head=git_head(),
        config_hash=_sha_short(MAIN_CFG),
        dry_run=args.dry_run,
        missing=status["missing_gpu"],
        running=status["running"],
    )

    completed: list[str] = []
    skipped: list[str] = []
    failures: list[dict[str, Any]] = []
    retries: dict[str, int] = {}
    pushed: list[str] = []

    try:
        while True:
            status = discover_status()
            missing = status["missing_gpu"]
            emit("loop", missing=missing, running=status["running"])

            if gpu0_busy_for_new_train():
                emit("wait_existing_gpu_owner", running=status["running"])
                if args.once:
                    break
                time.sleep(args.poll_seconds)
                continue

            # Train next missing arm (whole arm via runner; skip-if-done handles done folds).
            next_arm = None
            for arm in GPU_ARMS:
                arm_stat = status["arms"].get(arm["arm_id"]) or {}
                if int(arm_stat.get("n_ok") or 0) < 3:
                    next_arm = arm
                    break

            if next_arm is not None:
                key = str(next_arm["arm_id"])
                rc = train_arm(next_arm, dry_run=args.dry_run)
                # Re-validate
                time.sleep(5)
                status = discover_status()
                n_ok = int((status["arms"].get(key) or {}).get("n_ok") or 0)
                if n_ok >= 3:
                    completed.append(key)
                    refresh_reports(args.dry_run)
                elif rc != 0:
                    retries[key] = retries.get(key, 0) + 1
                    failures.append({"arm": key, "rc": rc, "attempt": retries[key]})
                    emit("train_failed", arm=key, rc=rc, attempt=retries[key])
                    if retries[key] >= 2:
                        emit("train_abandoned", arm=key)
                        # continue to other arms if any remain independent
                if args.once:
                    break
                continue

            # All matched GPU arms present → posthoc + report + decision
            emit("gpu_queue_complete")
            posthoc_enet_for_completed(dry_run=args.dry_run)
            refresh_reports(args.dry_run)
            sha = try_commit_push(
                "reports(7g-prime): complete matched 16-epoch promotion screen",
                dry_run=args.dry_run,
            )
            if sha:
                pushed.append(sha)

            seed_status = "skipped"
            if not args.skip_seed:
                seed_status = maybe_unlock_and_run_seed(dry_run=args.dry_run)
                if seed_status.startswith("seed_done"):
                    sha2 = try_commit_push(
                        "reports(7g-prime): evaluate age-primary seed-gene screen",
                        dry_run=args.dry_run,
                    )
                    if sha2:
                        pushed.append(sha2)

            summary = {
                "started_head": git_head(),
                "discovered_missing_at_start": status.get("missing_gpu"),
                "completed": completed,
                "skipped_already_valid": skipped,
                "failures": failures,
                "retries": retries,
                "seed_status": seed_status,
                "promotion_decision": (
                    json.loads((REPORT / "promotion_decision.json").read_text())
                    if (REPORT / "promotion_decision.json").is_file()
                    else None
                ),
                "report_paths": {
                    "analysis": str(REPORT / "analysis.md"),
                    "decision": str(REPORT / "promotion_decision.json"),
                    "weekend_summary": str(REPORT / "weekend_supervisor_summary.md"),
                },
                "pushed_shas": pushed,
                "hard_stop": True,
                "recommended_next_command": (
                    "Review reports/inspection/stage0_7g_gene_only_probe/analysis.md and "
                    "promotion_decision.json; do not launch Stage B automatically."
                ),
            }
            write_summary(summary)
            emit("hard_stop", summary=summary)
            break
    finally:
        try:
            fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
            lock_fh.close()
        except Exception:
            pass
        if PID_PATH.is_file():
            try:
                PID_PATH.unlink()
            except Exception:
                pass


if __name__ == "__main__":
    main()
