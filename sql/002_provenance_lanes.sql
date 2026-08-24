-- Provenance lanes for raw Stage 0 trees.
-- Do not merge Atlas / DataHub / CpGCorpus / manifests into one ingest path.

CREATE TABLE IF NOT EXISTS provenance_lane (
    source_system VARCHAR PRIMARY KEY,
    display_name VARCHAR NOT NULL,
    raw_root_relpath VARCHAR NOT NULL,
    role VARCHAR NOT NULL,
    default_file_formats VARCHAR NOT NULL,
    notes VARCHAR
);

INSERT OR REPLACE INTO provenance_lane VALUES
    (
        'cpgcorpus',
        'CpGCorpus',
        'raw/cpgcorpus',
        'training_matrices',
        'arrow',
        'Stage 0 GSEs only under raw/cpgcorpus/{GSE}/. Partial full-sync leftovers live in raw/cpgcorpus/_partial_fullsync/ and are not ingested by default.'
    ),
    (
        'ewas_atlas',
        'EWAS Atlas',
        'raw/ewas_atlas',
        'association_knowledge',
        'tsv',
        'Association / study / cohort / probe annotation batch exports. Never mixed with DataHub betas.'
    ),
    (
        'ewas_datahub_db',
        'EWAS DataHub All Data',
        'raw/ewas_datahub/EWAS_db',
        'sample_beta_text',
        'probe_beta_txt',
        'Per-study probe_id<TAB>beta text files. Separate release from baseline zip packs.'
    ),
    (
        'ewas_datahub_baseline',
        'EWAS DataHub Baseline',
        'raw/ewas_datahub/download',
        'baseline_trait_packs',
        'zip|zip_r',
        'Trait packs (*_methylation_v1.zip) and sample_*.zip (often R serialized). All packs live under download/ only — no flat ewas_datahub/*.zip.'
    ),
    (
        'epicv2_manifest',
        'EPICv2 reannotated manifest',
        'raw/manifests/epicv2',
        'array_annotation',
        'csv_gz',
        'Probe annotation only; not a methylation matrix source.'
    );

-- Existing catalogs created before source_system existed.
ALTER TABLE source_release ADD COLUMN IF NOT EXISTS source_system VARCHAR;

-- Milestone 7A columns for catalogs created before retrieved_at existed.
ALTER TABLE study ADD COLUMN IF NOT EXISTS retrieved_at TIMESTAMPTZ;
