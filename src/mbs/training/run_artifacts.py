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
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "model_state": model_state,
        "head_state": head_state,
        "optimizer_state": optimizer_state,
        "epoch": epoch,
        "metrics": metrics,
        "config_hash": config_hash,
    }
    torch.save(payload, path)
    return sha256_file(path)


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
