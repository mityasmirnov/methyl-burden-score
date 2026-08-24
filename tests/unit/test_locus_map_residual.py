"""Unit tests for matrix probe→locus map residual retention."""

from __future__ import annotations

import pandas as pd

from mbs.batch import ANNOTATION_STATUS_UNMAPPED
from mbs.matrix.locus_map import (
    build_probe_locus_map,
    is_residual_canonical_key,
    residual_canonical_key,
)


def test_build_probe_locus_map_retains_unmapped_probes() -> None:
    edges = pd.DataFrame(
        {
            "probe_id": ["cg001", "cg002"],
            "locus_id": [10, 20],
            "canonical_key": ["chr1:10", "chr1:20"],
            "genome_build": ["GRCh38", "GRCh38"],
            "platform_id": ["HM450", "HM450"],
            "is_primary": [True, True],
        }
    )
    observed = ["cg001", "cg002", "cg999"]
    locus_map = build_probe_locus_map(observed, edges, platform_id="HM450")
    assert locus_map.n_mapped_probes == 2
    assert locus_map.n_residual_probes == 1
    assert locus_map.unmapped_probe_ids == ("cg999",)
    assert len(locus_map.locus_ids) == 3
    assert is_residual_canonical_key(str(locus_map.canonical_keys[-1]))
    assert str(locus_map.canonical_keys[-1]) == residual_canonical_key("cg999")
    assert int(locus_map.annotation_status[-1]) == ANNOTATION_STATUS_UNMAPPED
    assert locus_map.contributing_probe_ids[-1] == ("cg999",)
    assert locus_map.collapse_method[-1] == "identity"
