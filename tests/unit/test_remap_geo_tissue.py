"""GEO tissue remapping from source_name / aliases."""

from __future__ import annotations

import pandas as pd

from mbs.geo_metadata import load_geo_tissue_aliases, remap_geo_tissue_frame
from mbs.training.phenotype_table import TissueOntology


def _tiny_ontology() -> TissueOntology:
    labels = (
        "whole blood",
        "peripheral blood mononuclear cell",
        "leukocyte",
        "brain",
        "brain - tumor",
        "breast",
        "breast - tumor",
        "CD14+ monocyte",
    )
    label_to_id = {lab: i for i, lab in enumerate(labels)}
    return TissueOntology(labels=labels, min_n=1, label_to_id=label_to_id)


def test_remap_recovers_source_name_pbmc() -> None:
    df = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "study_id": "GSE1",
                "source_name": "peripheral blood mononuclear cells",
                "tissue_raw": None,
                "tissue": None,
                "tissue_ontology_id": None,
                "tissue_map_status": "empty",
            },
            {
                "sample_id": "GSM2",
                "study_id": "GSE1",
                "source_name": "brain tumor",
                "tissue_raw": "brain tumor",
                "tissue": "brain tumor",
                "tissue_ontology_id": None,
                "tissue_map_status": "unmapped",
            },
        ]
    )
    aliases = load_geo_tissue_aliases()
    out, stats = remap_geo_tissue_frame(df, ontology=_tiny_ontology(), aliases=aliases)
    assert out.loc[0, "tissue_map_status"] == "mapped"
    assert out.loc[0, "tissue"] == "peripheral blood mononuclear cell"
    assert out.loc[1, "tissue_map_status"] == "mapped"
    assert out.loc[1, "tissue"] == "brain - tumor"
    assert stats["tissue_map_status_after"]["mapped"] == 2
