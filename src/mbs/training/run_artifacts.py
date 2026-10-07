"""Write resolved configs, metrics, and checkpoint manifests under artifact roots."""

from __future__ import annotations

import hashlib
import json
import os
import platform
from pathlib import Path
from typing import Any

import torch
import yaml

from mbs.annotation.manifest import write_json


class ResumeError(ValueError):
    """Fail-closed resume: truncated, mismatched, or incomplete checkpoint."""


def run_dir(artifact_root: Path, run_id: str) -> Path:
    return artifact_root / "runs" / run_id


def checkpoint_dir(artifact_root: Path, run_id: str) -> Path:
    return artifact_root / "checkpoints" / run_id


def write_resolved_config(path: Path, config: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def collect_environment(*, device: torch.device) -> dict[str, Any]:
    info: dict[str, Any] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "torch_version": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "device": str(device),
        "model_public_name": "deepMAT",
        "package_name": "methyl-burden-score",
        "cli_entrypoint": "mbs",
    }
    if torch.cuda.is_available():
        info["cuda_device_count"] = torch.cuda.device_count()
        # After CUDA_VISIBLE_DEVICES remapping we always train on logical device 0.
        idx = 0
        info["cuda_device_name"] = torch.cuda.get_device_name(idx)
        info["cuda_device_index"] = idx
    return info


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_checkpoint(
    path: Path,
    *,
    model_state: dict[str, Any],
    head_state: dict[str, Any],
    optimizer_state: dict[str, Any],
    epoch: int,
    metrics: dict[str, Any],
    config_hash: str,
    resume_state: dict[str, Any] | None = None,
) -> str:
    """Atomically write a checkpoint (temp file + ``os.replace``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "model_state": model_state,
        "head_state": head_state,
        "optimizer_state": optimizer_state,
        "epoch": epoch,
        "metrics": metrics,
        "config_hash": config_hash,
    }
    if resume_state is not None:
        payload["resume_state"] = resume_state
    # Unique temp beside the target so replace stays on the same filesystem.
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        torch.save(payload, tmp)
        os.replace(tmp, path)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    return sha256_file(path)


def resolve_resume_path(
    resume_cfg: Any,
    *,
    ckpt_root: Path,
) -> Path | None:
    """Return checkpoint path for resume, or None when resume is disabled.

    ``auto`` → ``ckpt_root/last.pt`` when that file exists (else fresh start).
    A concrete path string must exist.
    """
    if resume_cfg in (None, False, "", "false", "False", "off", "OFF", "none", "None"):
        return None
    if isinstance(resume_cfg, str) and resume_cfg.strip().lower() == "auto":
        candidate = ckpt_root / "last.pt"
        return candidate if candidate.is_file() else None
    path = Path(str(resume_cfg)).expanduser()
    if not path.is_file():
        raise ResumeError(f"training.resume path does not exist: {path}")
    return path


def load_resume_checkpoint(
    path: Path,
    *,
    expected_config_hash: str,
    max_epochs: int,
) -> dict[str, Any]:
    """Load ``last.pt`` for resume; raise :class:`ResumeError` on any mismatch."""
    try:
        payload = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:  # noqa: BLE001 — fail closed on truncated/corrupt
        raise ResumeError(f"cannot load checkpoint {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ResumeError(f"checkpoint {path} is not a dict payload")
    if "resume_state" not in payload:
        raise ResumeError(
            f"checkpoint {path} has no resume_state (pre-§3.6 or incomplete write); "
            "refusing to resume"
        )
    if not isinstance(payload["resume_state"], dict):
        raise ResumeError(f"checkpoint {path} resume_state must be a dict")
    got_hash = payload.get("config_hash")
    if got_hash != expected_config_hash:
        raise ResumeError(
            f"checkpoint config_hash mismatch for {path}: "
            f"checkpoint={got_hash!r} current={expected_config_hash!r}"
        )
    try:
        epoch = int(payload["epoch"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ResumeError(f"checkpoint {path} missing/invalid epoch") from exc
    if epoch < 1:
        raise ResumeError(f"checkpoint {path} has invalid epoch={epoch}")
    if epoch >= int(max_epochs):
        raise ResumeError(
            f"checkpoint {path} epoch={epoch} >= max_epochs={max_epochs}; nothing to resume"
        )
    for key in ("model_state", "head_state", "optimizer_state"):
        if key not in payload:
            raise ResumeError(f"checkpoint {path} missing {key}")
    return payload


def write_run_artifacts(
    *,
    run_root: Path,
    ckpt_root: Path,
    config: dict[str, Any],
    environment: dict[str, Any],
    metrics: dict[str, Any],
    split: dict[str, Any],
    checkpoint_hashes: dict[str, str],
    config_hash: str,
) -> None:
    run_root.mkdir(parents=True, exist_ok=True)
    ckpt_root.mkdir(parents=True, exist_ok=True)
    write_resolved_config(run_root / "resolved_config.yaml", config)
    write_json(run_root / "environment.json", environment)
    write_json(run_root / "metrics.json", metrics)
    write_json(run_root / "split.json", split)
    write_json(
        run_root / "checksums.json",
        {"config_hash": config_hash, "checkpoints": checkpoint_hashes},
    )
    write_json(
        ckpt_root / "checkpoint_manifest.json",
        {
            "config_hash": config_hash,
            "checkpoints": checkpoint_hashes,
        },
    )


def config_sha256(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()
