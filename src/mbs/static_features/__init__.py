"""Offline static locus feature export and artifact I/O.

Training must load immutable arrays from disk and must not import CpGPT.
"""

from __future__ import annotations

from mbs.static_features.coordinates import mbs_locus_to_cpgpt_location
from mbs.static_features.export_cpgpt import DEFAULT_FEATURE_SET_ID, export_cpgpt_adapter
from mbs.static_features.manifest import validate_static_feature_manifest
from mbs.static_features.store import StaticFeatureStorePaths, static_feature_store_paths

__all__ = [
    "DEFAULT_FEATURE_SET_ID",
    "StaticFeatureStorePaths",
    "export_cpgpt_adapter",
    "mbs_locus_to_cpgpt_location",
    "static_feature_store_paths",
    "validate_static_feature_manifest",
]
