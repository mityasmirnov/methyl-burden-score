-- Atlas association + seed-panel knowledge track (parallel, non-blocking).
--
-- Knowledge only: published CpG↔trait associations, an explicit locus→gene
-- edge view for seed selection, and derived per-trait seed gene panels.
-- No sample×CpG observations live here (ADR 0005: catalog stays independent of
-- the matrix backend). Atlas gene symbols are retained as SOURCE METADATA only
-- and never used as gene allocation (ADR 0006, ADR 0010). Seed panels use
-- explicit, non-nearest-gene edges only (ADR 0010 explicit_only; ADR 0004).

CREATE TABLE IF NOT EXISTS association_study (
    association_study_id VARCHAR PRIMARY KEY,
    source_database VARCHAR NOT NULL,
    source_release VARCHAR,
    pmid VARCHAR,
    doi VARCHAR,
    geo_accessions VARCHAR,       -- CSV or JSON string of GSE/GSM accessions
    cohort_name VARCHAR,
    platform VARCHAR,
    genome_build VARCHAR,
    sample_size BIGINT,
    tissue_cell_type VARCHAR,
    phenotype_definition VARCHAR,
    covariate_model VARCHAR,
    ancestry_summary VARCHAR,
    demographic_summary VARCHAR,
    atlas_study_id VARCHAR         -- ES* when sourced from the EWAS Atlas
);

CREATE TABLE IF NOT EXISTS trait_ontology (
    trait_id VARCHAR PRIMARY KEY,
    raw_trait_name VARCHAR NOT NULL,
    trait_family VARCHAR,
    value_type VARCHAR,            -- continuous|binary|multiclass
    units VARCHAR,
    tissue_conditionality VARCHAR,
    usage_class VARCHAR            -- core|auxiliary|evaluation_only
);

CREATE TABLE IF NOT EXISTS cpg_trait_association (
    association_id VARCHAR PRIMARY KEY,
    association_study_id VARCHAR NOT NULL,
    trait_id VARCHAR NOT NULL,
    probe_id VARCHAR,
    locus_id VARCHAR,
    genomic_chrom VARCHAR,
    genomic_pos BIGINT,
    effect DOUBLE,
    standard_error DOUBLE,
    p_value DOUBLE,
    fdr DOUBLE,
    effect_scale VARCHAR,
    effect_direction VARCHAR,
    sample_size BIGINT,
    discovery_replication VARCHAR,  -- discovery|replication
    atlas_gene_symbol VARCHAR,      -- SOURCE METADATA ONLY (never allocation)
    FOREIGN KEY (association_study_id) REFERENCES association_study(association_study_id),
    FOREIGN KEY (trait_id) REFERENCES trait_ontology(trait_id)
);

-- Explicit locus→gene edges for seed selection. is_nearest_gene must be FALSE
-- for any edge used to build a seed panel (ADR 0004/0010: no nearest-gene).
CREATE TABLE IF NOT EXISTS association_locus_gene_edge (
    locus_id VARCHAR NOT NULL,
    gene_id VARCHAR NOT NULL,
    gene_role VARCHAR NOT NULL,     -- promoter|body|UTR|other
    mapping_evidence VARCHAR,
    mapping_source VARCHAR,
    is_explicit BOOLEAN NOT NULL DEFAULT TRUE,
    is_nearest_gene BOOLEAN NOT NULL DEFAULT FALSE,
    cpg_context VARCHAR,
    PRIMARY KEY (locus_id, gene_id, gene_role)
);

CREATE TABLE IF NOT EXISTS seed_panel (
    seed_panel_id VARCHAR PRIMARY KEY,
    trait_id VARCHAR NOT NULL,
    fold_id VARCHAR,
    selection_source VARCHAR NOT NULL,  -- external_clean|internal_fold|hybrid_fold
    selection_method VARCHAR,
    thresholds_json JSON,
    database_release VARCHAR,
    excluded_study_ids VARCHAR,         -- CSV/JSON of leakage-excluded study IDs
    panel_hash VARCHAR,
    graph_hash VARCHAR,
    created_at TIMESTAMPTZ NOT NULL DEFAULT current_timestamp,
    FOREIGN KEY (trait_id) REFERENCES trait_ontology(trait_id)
);

CREATE TABLE IF NOT EXISTS seed_panel_gene (
    seed_panel_id VARCHAR NOT NULL,
    gene_id VARCHAR NOT NULL,
    rank BIGINT,
    score DOUBLE,
    n_associated_cpgs BIGINT,
    n_independent_studies BIGINT,
    direction_consistency DOUBLE,
    min_fdr DOUBLE,
    promoter_body_coverage DOUBLE,
    n_cpgs_in_matrix BIGINT,
    inclusion_reason VARCHAR,
    PRIMARY KEY (seed_panel_id, gene_id),
    FOREIGN KEY (seed_panel_id) REFERENCES seed_panel(seed_panel_id)
);

CREATE TABLE IF NOT EXISTS seed_panel_locus (
    seed_panel_id VARCHAR NOT NULL,
    gene_id VARCHAR NOT NULL,
    locus_id VARCHAR NOT NULL,
    gene_role VARCHAR,
    source_association_status VARCHAR,
    is_seed_cpg BOOLEAN NOT NULL DEFAULT TRUE,  -- vs gene-sibling enrichment
    PRIMARY KEY (seed_panel_id, gene_id, locus_id),
    FOREIGN KEY (seed_panel_id) REFERENCES seed_panel(seed_panel_id)
);
