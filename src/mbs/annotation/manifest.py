"""Annotation and graph manifest helpers."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    commit = result.stdout.strip()
    if not _COMMIT_RE.fullmatch(commit):
        raise RuntimeError(f"unexpected git commit: {commit!r}")
    return commit


def source_file_entry(
    name: str,
    path: Path,
    *,
    version: str,
    uri: str | None = None,
    license_note: str | None = None,
) -> dict[str, Any]:
    return {
        "name": name,
        "version": version,
        "uri": uri,
        "sha256": sha256_file(path),
        "license_note": license_note,
    }


def validate_graph_manifest(manifest: dict[str, Any]) -> None:
    """Validate required fields against ``schemas/graph_manifest.schema.json`` rules."""
    required = [
        "artifact_version",
        "graph_id",
        "genome_build",
        "builder_commit",
        "source_files",
        "region_policy",
        "genes_path",
        "regions_path",
        "locus_region_edges_path",
        "regions_bed_path",
        "created_at",
    ]
    missing = [key for key in required if key not in manifest]
    if missing:
        raise ValueError(f"graph manifest missing keys: {missing}")
    if manifest["genome_build"] != "GRCh38":
        raise ValueError("genome_build must be GRCh38")
    if not _COMMIT_RE.fullmatch(str(manifest["builder_commit"])):
        raise ValueError("builder_commit must be a 40-char lowercase hex SHA")
    if not isinstance(manifest["source_files"], list) or not manifest["source_files"]:
        raise ValueError("source_files must be a non-empty list")
    for entry in manifest["source_files"]:
        for key in ("name", "version", "sha256"):
            if key not in entry:
                raise ValueError(f"source_files entry missing {key}")
        if not _SHA256_RE.fullmatch(str(entry["sha256"])):
            raise ValueError(f"invalid sha256 for {entry.get('name')}")
    policy = manifest["region_policy"]
    for key in ("region_types", "role_precedence"):
        if key not in policy or not isinstance(policy[key], list) or not policy[key]:
            raise ValueError(f"region_policy.{key} must be a non-empty list")
    for key in (
        "genes_path",
        "regions_path",
        "locus_region_edges_path",
        "regions_bed_path",
    ):
        path = str(manifest[key])
        if not path.startswith("/data/"):
            raise ValueError(f"{key} must be an absolute /data path, got {path}")
    optional = manifest.get("region_gene_edges_path")
    if optional is not None and not str(optional).startswith("/data/"):
        raise ValueError("region_gene_edges_path must be an absolute /data path or null")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def utc_now_iso() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
