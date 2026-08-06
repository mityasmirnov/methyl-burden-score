"""Lightweight terminal dashboard for live / finished training runs."""

from __future__ import annotations

import json
import os
import re
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml
from rich.console import Group, RenderableType
from rich.live import Live
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from mbs.training.run_artifacts import checkpoint_dir, run_dir

_RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_PS_BIN = "/usr/bin/ps"
_NVIDIA_SMI_BIN = "/usr/bin/nvidia-smi"
_DEFAULT_TB_PORT = 6006
_TB_META_NAME = "tensorboard.json"
_PORT_RE = re.compile(r"(?:--port(?:=|\s+)|-p\s+)(\d+)")


@dataclass(frozen=True)
class EpochMetrics:
    epoch: int
    train_loss: float | None
    val_loss: float | None
    train_mae: float | None
    val_mae: float | None
    train_accuracy: float | None
    val_accuracy: float | None
    train_macro_f1: float | None
    val_macro_f1: float | None
    learning_rate: float | None
    task: str | None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class GpuSnapshot:
    index: int
    name: str
    util_pct: float | None
    mem_used_mib: float | None
    mem_total_mib: float | None


@dataclass(frozen=True)
class TensorBoardServer:
    port: int
    pid: int | None
    logdir: Path
    url: str
    reused: bool
    meta_path: Path


@dataclass
class MonitorSnapshot:
    run_id: str
    run_root: Path
    ckpt_root: Path
    jsonl_path: Path
    history: list[EpochMetrics]
    latest: EpochMetrics | None
    max_epochs: int | None
    status: str  # running | finished | waiting | stalled
    train_pid: int | None
    best_ckpt: Path | None
    last_ckpt: Path | None
    gpus: list[GpuSnapshot]
    eta_seconds: float | None
    seconds_per_epoch: float | None
    jsonl_mtime: float | None
    note: str | None = None
    tensorboard: TensorBoardServer | None = None


def validate_run_id(run_id: str) -> str:
    if not _RUN_ID_RE.fullmatch(run_id):
        raise ValueError(
            "run_id must be a single path segment matching "
            r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$"
        )
    return run_id


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def port_is_free(port: int, *, host: str = "127.0.0.1") -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind((host, port))
        except OSError:
            return False
    return True


def find_free_port(preferred: int = _DEFAULT_TB_PORT, *, attempts: int = 40) -> int:
    if preferred < 1 or preferred > 65535:
        raise ValueError("preferred port out of range")
    for offset in range(attempts):
        port = preferred + offset
        if port > 65535:
            break
        if port_is_free(port):
            return port
    raise RuntimeError(f"no free TCP port near {preferred}")


def _tensorboard_bin() -> Path:
    candidate = Path(sys.executable).resolve().parent / "tensorboard"
    if candidate.is_file():
        return candidate
    raise RuntimeError(
        "tensorboard executable not found next to the active Python. "
        "Install the training extra: uv sync --extra training"
    )


def _parse_port_from_args(args: str) -> int | None:
    match = _PORT_RE.search(args)
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def find_running_tensorboard(logdir: Path) -> TensorBoardServer | None:
    """Return a live TensorBoard process whose cmdline includes ``logdir``."""
    target = str(logdir.resolve())
    try:
        proc = subprocess.run(  # noqa: S603
            [_PS_BIN, "-eo", "pid=,args="],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    for raw_line in proc.stdout.splitlines():
        line = raw_line.strip()
        if not line or "tensorboard" not in line.lower():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        pid_s, args = parts
        if target not in args and str(logdir) not in args:
            continue
        try:
            pid = int(pid_s)
        except ValueError:
            continue
        port = _parse_port_from_args(args) or _DEFAULT_TB_PORT
        return TensorBoardServer(
            port=port,
            pid=pid,
            logdir=logdir.resolve(),
            url=f"http://127.0.0.1:{port}",
            reused=True,
            meta_path=logdir.parent / _TB_META_NAME,
        )
    return None


def read_tensorboard_meta(run_root: Path) -> TensorBoardServer | None:
    path = run_root / _TB_META_NAME
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    server: TensorBoardServer | None = None
    if isinstance(payload, dict):
        try:
            port = int(payload["port"])
            logdir = Path(str(payload["logdir"]))
            pid_raw = payload.get("pid")
            pid = int(pid_raw) if pid_raw is not None else None
            alive = pid is not None and _pid_alive(pid)
            port_held = pid is None and not port_is_free(port)
            if alive or port_held:
                server = TensorBoardServer(
                    port=port,
                    pid=pid,
                    logdir=logdir,
                    url=str(payload.get("url") or f"http://127.0.0.1:{port}"),
                    reused=True,
                    meta_path=path,
                )
        except (KeyError, TypeError, ValueError):
            server = None
    return server


def write_tensorboard_meta(run_root: Path, server: TensorBoardServer) -> Path:
    path = run_root / _TB_META_NAME
    payload = {
        "port": server.port,
        "pid": server.pid,
        "logdir": str(server.logdir),
        "url": server.url,
        "ssh_tunnel": ssh_tunnel_hint(server.port),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def ssh_tunnel_hint(port: int) -> str:
    return f"ssh -L {port}:localhost:{port} <user>@<power-horse-host>"


def ensure_tensorboard(
    *,
    run_root: Path,
    logdir: Path | None = None,
    preferred_port: int = _DEFAULT_TB_PORT,
) -> TensorBoardServer:
    """Start or reuse TensorBoard for a run's ``tb/`` directory."""
    run_root = run_root.resolve()
    tb_dir = (logdir or (run_root / "tb")).resolve()
    tb_dir.mkdir(parents=True, exist_ok=True)

    existing = read_tensorboard_meta(run_root) or find_running_tensorboard(tb_dir)
    if existing is not None:
        server = TensorBoardServer(
            port=existing.port,
            pid=existing.pid,
            logdir=tb_dir,
            url=f"http://127.0.0.1:{existing.port}",
            reused=True,
            meta_path=run_root / _TB_META_NAME,
        )
        write_tensorboard_meta(run_root, server)
        return server

    port = preferred_port if port_is_free(preferred_port) else find_free_port(preferred_port)
    tb_bin = _tensorboard_bin()
    log_path = run_root / "tensorboard.log"
    with log_path.open("a", encoding="utf-8") as log_handle:
        proc = subprocess.Popen(  # noqa: S603
            [
                str(tb_bin),
                "--logdir",
                str(tb_dir),
                "--host",
                "0.0.0.0",  # noqa: S104 — bind for SSH tunnel from laptop
                "--port",
                str(port),
            ],
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    time.sleep(0.4)
    if proc.poll() is not None:
        raise RuntimeError(
            f"TensorBoard exited immediately (code {proc.returncode}). See {log_path}"
        )
    server = TensorBoardServer(
        port=port,
        pid=proc.pid,
        logdir=tb_dir,
        url=f"http://127.0.0.1:{port}",
        reused=False,
        meta_path=run_root / _TB_META_NAME,
    )
    write_tensorboard_meta(run_root, server)
    return server


def _opt_float(row: dict[str, Any], *keys: str) -> float | None:
    for key in keys:
        if key in row and row[key] is not None:
            try:
                return float(row[key])
            except (TypeError, ValueError):
                return None
    return None


def parse_metrics_row(row: dict[str, Any]) -> EpochMetrics | None:
    epoch_raw = row.get("epoch")
    if epoch_raw is None:
        return None
    try:
        epoch = int(epoch_raw)
    except (TypeError, ValueError):
        return None
    return EpochMetrics(
        epoch=epoch,
        train_loss=_opt_float(row, "train_loss"),
        val_loss=_opt_float(row, "val_loss"),
        train_mae=_opt_float(row, "train_mae"),
        val_mae=_opt_float(row, "val_mae"),
        train_accuracy=_opt_float(row, "train_accuracy"),
        val_accuracy=_opt_float(row, "val_accuracy"),
        train_macro_f1=_opt_float(row, "train_macro_f1", "train_f1", "train_macro_F1"),
        val_macro_f1=_opt_float(row, "val_macro_f1", "val_f1", "val_macro_F1"),
        learning_rate=_opt_float(row, "learning_rate", "lr"),
        task=str(row["task"]) if row.get("task") is not None else None,
        raw=row,
    )


def read_metrics_jsonl(path: Path) -> list[EpochMetrics]:
    if not path.is_file():
        return []
    out: list[EpochMetrics] = []
    text = path.read_text(encoding="utf-8")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(row, dict):
            continue
        parsed = parse_metrics_row(row)
        if parsed is not None:
            out.append(parsed)
    return out


def resolve_max_epochs(
    *,
    run_root: Path,
    config_path: Path | None,
    max_epochs_override: int | None,
) -> int | None:
    if max_epochs_override is not None and max_epochs_override >= 1:
        return int(max_epochs_override)
    for candidate in (
        run_root / "resolved_config.yaml",
        config_path,
    ):
        if candidate is None or not candidate.is_file():
            continue
        try:
            payload = yaml.safe_load(candidate.read_text(encoding="utf-8"))
        except (OSError, yaml.YAMLError):
            continue
        if not isinstance(payload, dict):
            continue
        training = payload.get("training")
        if isinstance(training, dict) and training.get("max_epochs") is not None:
            try:
                value = int(training["max_epochs"])
            except (TypeError, ValueError):
                continue
            if value >= 1:
                return value
    return None


def find_train_pid(run_id: str) -> int | None:
    """Return PID of a live `mbs train` process that mentions run_id, if any."""
    try:
        proc = subprocess.run(  # noqa: S603
            [_PS_BIN, "-eo", "pid=,args="],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    needle = run_id
    for raw_line in proc.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        pid_s, args = parts
        if "mbs" not in args or "train" not in args:
            continue
        if needle not in args:
            continue
        if " monitor" in f" {args}":
            continue
        try:
            return int(pid_s)
        except ValueError:
            continue
    return None


def query_gpus() -> list[GpuSnapshot]:
    try:
        proc = subprocess.run(  # noqa: S603
            [
                _NVIDIA_SMI_BIN,
                "--query-gpu=index,name,utilization.gpu,memory.used,memory.total",
                "--format=csv,noheader,nounits",
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    if proc.returncode != 0:
        return []
    gpus: list[GpuSnapshot] = []
    for line in proc.stdout.splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 5:
            continue
        try:
            index = int(parts[0])
        except ValueError:
            continue

        def _f(value: str) -> float | None:
            try:
                return float(value)
            except ValueError:
                return None

        gpus.append(
            GpuSnapshot(
                index=index,
                name=parts[1],
                util_pct=_f(parts[2]),
                mem_used_mib=_f(parts[3]),
                mem_total_mib=_f(parts[4]),
            )
        )
    return gpus


def _ckpt_if_exists(path: Path) -> Path | None:
    return path if path.is_file() else None


def estimate_eta(
    *,
    history: list[EpochMetrics],
    max_epochs: int | None,
    epoch_timestamps: dict[int, float],
) -> tuple[float | None, float | None]:
    """Return (eta_seconds, seconds_per_epoch) from observed epoch wall times."""
    if len(epoch_timestamps) < 2:
        return None, None
    ordered = sorted(epoch_timestamps.items())
    deltas: list[float] = []
    for (e0, t0), (e1, t1) in zip(ordered, ordered[1:], strict=False):
        if e1 <= e0:
            continue
        dt = t1 - t0
        if dt > 0:
            deltas.append(dt / (e1 - e0))
    if not deltas:
        return None, None
    window = deltas[-5:]
    sec_per = sum(window) / len(window)
    if max_epochs is None or not history:
        return None, sec_per
    remaining = max_epochs - history[-1].epoch
    if remaining <= 0:
        return 0.0, sec_per
    return remaining * sec_per, sec_per


def collect_snapshot(
    *,
    run_id: str,
    artifact_root: Path,
    config_path: Path | None = None,
    max_epochs_override: int | None = None,
    epoch_timestamps: dict[int, float] | None = None,
    tensorboard: TensorBoardServer | None = None,
) -> MonitorSnapshot:
    run_id = validate_run_id(run_id)
    run_root = run_dir(artifact_root, run_id)
    ckpt_root = checkpoint_dir(artifact_root, run_id)
    jsonl_path = run_root / "metrics.jsonl"
    history = read_metrics_jsonl(jsonl_path)
    latest = history[-1] if history else None
    max_epochs = resolve_max_epochs(
        run_root=run_root,
        config_path=config_path,
        max_epochs_override=max_epochs_override,
    )
    stamps = dict(epoch_timestamps or {})
    if latest is not None and latest.epoch not in stamps:
        stamps[latest.epoch] = time.time()
    eta, spe = estimate_eta(history=history, max_epochs=max_epochs, epoch_timestamps=stamps)

    train_pid = find_train_pid(run_id)
    finished = (run_root / "metrics.json").is_file()
    note: str | None = None
    if finished:
        status = "finished"
    elif train_pid is not None:
        status = "running"
    elif history:
        status = "stalled"
        note = "metrics.jsonl present but no matching train process"
    elif run_root.is_dir():
        status = "waiting"
        note = "run directory exists; waiting for first metrics.jsonl row"
    else:
        status = "waiting"
        note = f"run directory missing: {run_root}"

    tb = tensorboard or read_tensorboard_meta(run_root)
    jsonl_mtime = jsonl_path.stat().st_mtime if jsonl_path.is_file() else None
    return MonitorSnapshot(
        run_id=run_id,
        run_root=run_root,
        ckpt_root=ckpt_root,
        jsonl_path=jsonl_path,
        history=history,
        latest=latest,
        max_epochs=max_epochs,
        status=status,
        train_pid=train_pid,
        best_ckpt=_ckpt_if_exists(ckpt_root / "best.pt"),
        last_ckpt=_ckpt_if_exists(ckpt_root / "last.pt"),
        gpus=query_gpus(),
        eta_seconds=eta,
        seconds_per_epoch=spe,
        jsonl_mtime=jsonl_mtime,
        note=note,
        tensorboard=tb,
    )


def _fmt(value: float | None, *, pct: bool = False, digits: int = 4) -> str:
    if value is None:
        return "—"
    if pct:
        return f"{100.0 * value:6.2f}%"
    return f"{value:.{digits}f}"


def _fmt_duration(seconds: float | None) -> str:
    if seconds is None:
        return "—"
    if seconds < 0:
        seconds = 0.0
    total = round(seconds)
    h, rem = divmod(total, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h:d}h {m:02d}m {s:02d}s"
    if m:
        return f"{m:d}m {s:02d}s"
    return f"{s:d}s"


def render_snapshot(snap: MonitorSnapshot) -> RenderableType:
    status_style = {
        "running": "bold green",
        "finished": "bold cyan",
        "waiting": "yellow",
        "stalled": "bold red",
    }.get(snap.status, "white")

    header = Table.grid(expand=True)
    header.add_column(ratio=3)
    header.add_column(justify="right", ratio=2)
    left = Text()
    left.append("run ", style="dim")
    left.append(snap.run_id, style="bold")
    if snap.latest and snap.latest.task:
        left.append(f"  ·  {snap.latest.task}", style="dim")
    right = Text(snap.status.upper(), style=status_style)
    if snap.train_pid is not None:
        right.append(f"  pid={snap.train_pid}", style="dim")
    header.add_row(left, right)

    metrics = Table(show_header=True, header_style="bold", expand=True, box=None)
    metrics.add_column("metric")
    metrics.add_column("train", justify="right")
    metrics.add_column("val", justify="right")

    latest = snap.latest
    epoch_label = "—"
    if latest is not None:
        if snap.max_epochs is not None:
            epoch_label = f"{latest.epoch} / {snap.max_epochs}"
        else:
            epoch_label = str(latest.epoch)

    metrics.add_row("epoch", epoch_label, "")
    metrics.add_row(
        "loss",
        _fmt(latest.train_loss if latest else None),
        _fmt(latest.val_loss if latest else None),
    )
    metrics.add_row(
        "age MAE",
        _fmt(latest.train_mae if latest else None),
        _fmt(latest.val_mae if latest else None),
    )
    metrics.add_row(
        "tissue acc",
        _fmt(latest.train_accuracy if latest else None, pct=True),
        _fmt(latest.val_accuracy if latest else None, pct=True),
    )
    metrics.add_row(
        "tissue macro-F1",
        _fmt(latest.train_macro_f1 if latest else None),
        _fmt(latest.val_macro_f1 if latest else None),
    )
    if latest and latest.learning_rate is not None:
        metrics.add_row("lr", _fmt(latest.learning_rate, digits=6), "")

    meta = Table(show_header=False, expand=True, box=None)
    meta.add_column(style="dim")
    meta.add_column()
    meta.add_row("sec / epoch", _fmt_duration(snap.seconds_per_epoch))
    meta.add_row("ETA", _fmt_duration(snap.eta_seconds))
    meta.add_row("best.pt", str(snap.best_ckpt) if snap.best_ckpt else "—")
    meta.add_row("last.pt", str(snap.last_ckpt) if snap.last_ckpt else "—")
    meta.add_row("metrics.jsonl", str(snap.jsonl_path))
    if snap.tensorboard is not None:
        meta.add_row("TensorBoard", snap.tensorboard.url)
        meta.add_row("SSH tunnel", ssh_tunnel_hint(snap.tensorboard.port))
    else:
        meta.add_row("TensorBoard", "—")

    gpu_table = Table(show_header=True, header_style="bold", expand=True, box=None)
    gpu_table.add_column("#", justify="right")
    gpu_table.add_column("name")
    gpu_table.add_column("util", justify="right")
    gpu_table.add_column("memory", justify="right")
    if snap.gpus:
        for gpu in snap.gpus:
            util = "—" if gpu.util_pct is None else f"{gpu.util_pct:.0f}%"
            if gpu.mem_used_mib is None or gpu.mem_total_mib is None:
                mem = "—"
            else:
                mem = f"{gpu.mem_used_mib:,.0f} / {gpu.mem_total_mib:,.0f} MiB"
            gpu_table.add_row(str(gpu.index), gpu.name, util, mem)
    else:
        gpu_table.add_row("—", "nvidia-smi unavailable", "—", "—")

    body: list[RenderableType] = [
        header,
        Text(""),
        Panel(metrics, title="metrics", border_style="blue"),
        Panel(meta, title="progress", border_style="blue"),
        Panel(gpu_table, title="GPU", border_style="blue"),
    ]
    if snap.note:
        body.append(Text(snap.note, style="yellow"))
    if (
        latest
        and latest.val_accuracy == 0.0
        and (latest.task or "")
        in {
            "multitask",
            "multiclass",
        }
    ):
        body.append(
            Text(
                "note: val tissue acc 0% can be expected for disjoint study-holdout classes",
                style="dim",
            )
        )
    return Group(*body)


def run_monitor(
    *,
    run_id: str,
    artifact_root: Path,
    config_path: Path | None = None,
    max_epochs: int | None = None,
    interval_s: float = 2.0,
    once: bool = False,
    start_tensorboard: bool = True,
    tb_port: int = _DEFAULT_TB_PORT,
) -> MonitorSnapshot:
    """Start TensorBoard (default) and render a live Rich TUI until finished."""
    if interval_s < 0.2:
        raise ValueError("interval_s must be >= 0.2")

    run_id = validate_run_id(run_id)
    run_root = run_dir(artifact_root, run_id)
    run_root.mkdir(parents=True, exist_ok=True)
    tb_server: TensorBoardServer | None = None
    if start_tensorboard:
        tb_server = ensure_tensorboard(
            run_root=run_root,
            preferred_port=tb_port,
        )

    epoch_timestamps: dict[int, float] = {}
    snap = collect_snapshot(
        run_id=run_id,
        artifact_root=artifact_root,
        config_path=config_path,
        max_epochs_override=max_epochs,
        epoch_timestamps=epoch_timestamps,
        tensorboard=tb_server,
    )
    if snap.latest is not None:
        epoch_timestamps[snap.latest.epoch] = time.time()

    with Live(
        render_snapshot(snap),
        refresh_per_second=max(1, int(1 / min(interval_s, 1.0))),
        screen=False,
    ) as live:
        while True:
            snap = collect_snapshot(
                run_id=run_id,
                artifact_root=artifact_root,
                config_path=config_path,
                max_epochs_override=max_epochs,
                epoch_timestamps=epoch_timestamps,
                tensorboard=tb_server,
            )
            if snap.latest is not None and snap.latest.epoch not in epoch_timestamps:
                epoch_timestamps[snap.latest.epoch] = time.time()
            eta, spe = estimate_eta(
                history=snap.history,
                max_epochs=snap.max_epochs,
                epoch_timestamps=epoch_timestamps,
            )
            snap.eta_seconds = eta
            snap.seconds_per_epoch = spe
            live.update(render_snapshot(snap))
            if once or snap.status == "finished":
                break
            time.sleep(interval_s)
    return snap
