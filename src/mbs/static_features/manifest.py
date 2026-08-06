"""Static feature artifact manifest helpers."""

from __future__ import annotations

import re
from typing import Any

from mbs.annotation.manifest import utc_now_iso, write_json

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
_DATA_PATH_RE = re.compile(r"^/data/")

ARTIFACT_VERSION = "static-feature-manifest-v1"
SOURCE_REPOSITORY = "https://github.com/lucascamillomd/CpGPT"


def validate_static_feature_manifest(manifest: dict[str, Any]) -> None:
    """Validate required fields against ``schemas/static_feature_manifest.schema.json``."""
    required = [
        "artifact_version",
        "feature_set_id",
        "source_model",
        "source_repository",
        "source_commit",
        "checkpoint_sha256",
        "genome_build",
        "output_dimension",
        "storage_dtype",
        "locus_table_sha256",
        "embedding_path",
        "locus_index_path",
        "export_command",
        "created_at",
    ]
    missing = [key for key in required if key not in manifest]
    if missing:
        raise ValueError(f"static feature manifest missing keys: {missing}")
    if not str(manifest["feature_set_id"]).strip():
        raise ValueError("feature_set_id must be non-empty")
    if not str(manifest["source_model"]).strip():
        raise ValueError("source_model must be non-empty")
    if not str(manifest["source_repository"]).startswith(("http://", "https://")):
        raise ValueError("source_repository must be an absolute URI")
    if not _COMMIT_RE.fullmatch(str(manifest["source_commit"])):
        raise ValueError("source_commit must be a 40-char lowercase hex SHA")
    if not _SHA256_RE.fullmatch(str(manifest["checkpoint_sha256"])):
        raise ValueError("checkpoint_sha256 must be a 64-char lowercase hex SHA")
    config_hash = manifest.get("configuration_sha256")
    if config_hash is not None and not _SHA256_RE.fullmatch(str(config_hash)):
        raise ValueError("configuration_sha256 must be a 64-char lowercase hex SHA or null")
    vocab_hash = manifest.get("vocabulary_sha256")
    if vocab_hash is not None and not _SHA256_RE.fullmatch(str(vocab_hash)):
        raise ValueError("vocabulary_sha256 must be a 64-char lowercase hex SHA or null")
    if manifest["genome_build"] != "GRCh38":
        raise ValueError("genome_build must be GRCh38")
    if not isinstance(manifest["output_dimension"], int) or manifest["output_dimension"] < 1:
        raise ValueError("output_dimension must be a positive integer")
    if manifest["storage_dtype"] not in {"float16", "bfloat16", "float32"}:
        raise ValueError(f"unsupported storage_dtype: {manifest['storage_dtype']}")
    if not _SHA256_RE.fullmatch(str(manifest["locus_table_sha256"])):
        raise ValueError("locus_table_sha256 must be a 64-char lowercase hex SHA")
    for key in ("embedding_path", "locus_index_path"):
        path = str(manifest[key])
        if not _DATA_PATH_RE.match(path):
            raise ValueError(f"{key} must be an absolute /data path, got {path}")
    if not str(manifest["export_command"]).strip():
        raise ValueError("export_command must be non-empty")
    if not str(manifest["created_at"]).strip():
        raise ValueError("created_at must be non-empty")
    context = manifest.get("context_length")
    if context is not None and (not isinstance(context, int) or context < 1):
        raise ValueError("context_length must be a positive integer or null")
    input_dim = manifest.get("input_dimension")
    if input_dim is not None and (not isinstance(input_dim, int) or input_dim < 1):
        raise ValueError("input_dimension must be a positive integer or null")
    n_loci = manifest.get("n_loci")
    if n_loci is not None and (not isinstance(n_loci, int) or n_loci < 0):
        raise ValueError("n_loci must be a non-negative integer or absent")


def write_static_feature_manifest(path: Any, manifest: dict[str, Any]) -> None:
    validate_static_feature_manifest(manifest)
    write_json(path, manifest)


__all__ = [
    "ARTIFACT_VERSION",
    "SOURCE_REPOSITORY",
    "utc_now_iso",
    "validate_static_feature_manifest",
    "write_static_feature_manifest",
]
