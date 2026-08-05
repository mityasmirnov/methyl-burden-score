-- Core analytical catalog schema.
-- Large sample-by-locus matrices remain in Zarr or HDF5 and are referenced by manifests.

CREATE TABLE IF NOT EXISTS source_release (
    source_release_id VARCHAR PRIMARY KEY,
    source_name VARCHAR NOT NULL,
    source_version VARCHAR,
    retrieved_at TIMESTAMPTZ,
    source_uri VARCHAR,
    license_note VARCHAR,
    manifest_sha256 VARCHAR
);

CREATE TABLE IF NOT EXISTS platform (
    platform_id VARCHAR PRIMARY KEY,
    platform_name VARCHAR NOT NULL,
    manufacturer VARCHAR,
    manifest_version VARCHAR,
    nominal_probe_count BIGINT,
    genome_build VARCHAR
);

CREATE TABLE IF NOT EXISTS study (
    study_id VARCHAR PRIMARY KEY,
    source_release_id VARCHAR,
    gse_id VARCHAR,
    cohort_id VARCHAR,
    platform_id VARCHAR,
    processing_level VARCHAR,
    genome_build VARCHAR,
    metadata_json JSON,
    FOREIGN KEY (source_release_id) REFERENCES source_release(source_release_id),
    FOREIGN KEY (platform_id) REFERENCES platform(platform_id)
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id VARCHAR PRIMARY KEY,
    study_id VARCHAR NOT NULL,
    source_sample_id VARCHAR,
    donor_id VARCHAR,
    replicate_group VARCHAR,
    age DOUBLE,
    sex VARCHAR,
    tissue_raw VARCHAR,
    tissue_ontology_id VARCHAR,
    case_control VARCHAR,
    metadata_json JSON,
    FOREIGN KEY (study_id) REFERENCES study(study_id)
);

CREATE TABLE IF NOT EXISTS assay_file (
    assay_file_id VARCHAR PRIMARY KEY,
    study_id VARCHAR NOT NULL,
    path VARCHAR NOT NULL,
    format VARCHAR NOT NULL,
    sha256 VARCHAR NOT NULL,
    byte_size UBIGINT,
    n_rows UBIGINT,
    n_columns UBIGINT,
    matrix_orientation VARCHAR,
    schema_hash VARCHAR,
    processing_level VARCHAR,
    FOREIGN KEY (study_id) REFERENCES study(study_id)
);

CREATE TABLE IF NOT EXISTS phenotype (
    phenotype_id VARCHAR PRIMARY KEY,
    phenotype_name VARCHAR NOT NULL,
    phenotype_type VARCHAR NOT NULL,
    ontology_id VARCHAR,
    unit VARCHAR
);

CREATE TABLE IF NOT EXISTS sample_phenotype (
    sample_id VARCHAR NOT NULL,
    phenotype_id VARCHAR NOT NULL,
    numeric_value DOUBLE,
    categorical_value VARCHAR,
    is_observed BOOLEAN NOT NULL DEFAULT TRUE,
    PRIMARY KEY (sample_id, phenotype_id),
    FOREIGN KEY (sample_id) REFERENCES sample(sample_id),
    FOREIGN KEY (phenotype_id) REFERENCES phenotype(phenotype_id)
);

CREATE TABLE IF NOT EXISTS probe (
    probe_id VARCHAR NOT NULL,
    platform_id VARCHAR NOT NULL,
    probe_design VARCHAR,
    manifest_version VARCHAR,
    quality_flags JSON,
    PRIMARY KEY (probe_id, platform_id),
    FOREIGN KEY (platform_id) REFERENCES platform(platform_id)
);

CREATE TABLE IF NOT EXISTS locus (
    locus_id UBIGINT PRIMARY KEY,
    genome_build VARCHAR NOT NULL,
    chromosome VARCHAR NOT NULL,
    position UBIGINT NOT NULL,
    canonical_key VARCHAR NOT NULL UNIQUE,
    mapping_status VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS probe_locus_edge (
    probe_id VARCHAR NOT NULL,
    platform_id VARCHAR NOT NULL,
    locus_id UBIGINT NOT NULL,
    mapping_source VARCHAR NOT NULL,
    mapping_confidence DOUBLE,
    is_primary BOOLEAN NOT NULL DEFAULT TRUE,
    PRIMARY KEY (probe_id, platform_id, locus_id),
    FOREIGN KEY (probe_id, platform_id) REFERENCES probe(probe_id, platform_id),
    FOREIGN KEY (locus_id) REFERENCES locus(locus_id)
);

CREATE TABLE IF NOT EXISTS gene (
    gene_id VARCHAR PRIMARY KEY,
    gene_name VARCHAR,
    chromosome VARCHAR NOT NULL,
    start UBIGINT NOT NULL,
    "end" UBIGINT NOT NULL,
    strand VARCHAR NOT NULL,
    gene_type VARCHAR,
    source_version VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS region (
    region_id VARCHAR PRIMARY KEY,
    gene_id VARCHAR NOT NULL,
    region_type VARCHAR NOT NULL,
    chromosome VARCHAR NOT NULL,
    start UBIGINT NOT NULL,
    "end" UBIGINT NOT NULL,
    strand VARCHAR NOT NULL,
    source VARCHAR NOT NULL,
    source_version VARCHAR NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES gene(gene_id)
);

CREATE TABLE IF NOT EXISTS locus_region_edge (
    locus_id UBIGINT NOT NULL,
    region_id VARCHAR NOT NULL,
    edge_weight DOUBLE NOT NULL DEFAULT 1.0,
    evidence_type VARCHAR NOT NULL,
    primary_gene_role BOOLEAN NOT NULL DEFAULT TRUE,
    PRIMARY KEY (locus_id, region_id),
    FOREIGN KEY (locus_id) REFERENCES locus(locus_id),
    FOREIGN KEY (region_id) REFERENCES region(region_id)
);

CREATE TABLE IF NOT EXISTS capsule (
    capsule_id VARCHAR PRIMARY KEY,
    capsule_type VARCHAR NOT NULL,
    parent_capsule_id VARCHAR,
    label VARCHAR NOT NULL,
    source VARCHAR NOT NULL,
    source_version VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS locus_capsule_edge (
    locus_id UBIGINT NOT NULL,
    capsule_id VARCHAR NOT NULL,
    edge_weight DOUBLE NOT NULL DEFAULT 1.0,
    evidence_type VARCHAR NOT NULL,
    PRIMARY KEY (locus_id, capsule_id),
    FOREIGN KEY (locus_id) REFERENCES locus(locus_id),
    FOREIGN KEY (capsule_id) REFERENCES capsule(capsule_id)
);

CREATE TABLE IF NOT EXISTS fold_assignment (
    split_id VARCHAR NOT NULL,
    sample_id VARCHAR NOT NULL,
    outer_fold INTEGER NOT NULL,
    inner_fold INTEGER,
    group_id VARCHAR NOT NULL,
    role VARCHAR NOT NULL,
    PRIMARY KEY (split_id, sample_id),
    FOREIGN KEY (sample_id) REFERENCES sample(sample_id)
);

CREATE TABLE IF NOT EXISTS artifact (
    artifact_id VARCHAR PRIMARY KEY,
    artifact_type VARCHAR NOT NULL,
    path VARCHAR NOT NULL,
    sha256 VARCHAR,
    manifest_json JSON NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT current_timestamp
);

CREATE TABLE IF NOT EXISTS experiment (
    experiment_id VARCHAR PRIMARY KEY,
    git_commit VARCHAR NOT NULL,
    resolved_config_path VARCHAR NOT NULL,
    data_release_id VARCHAR NOT NULL,
    graph_artifact_id VARCHAR,
    feature_artifact_ids VARCHAR[],
    split_id VARCHAR,
    status VARCHAR NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT current_timestamp,
    FOREIGN KEY (graph_artifact_id) REFERENCES artifact(artifact_id)
);
