-- Study-level EWAS Atlas enrichment (external stratification; not sample labels).

CREATE TABLE IF NOT EXISTS study_atlas_enrichment (
    study_id VARCHAR PRIMARY KEY,
    join_method VARCHAR NOT NULL,
    atlas_study_ids VARCHAR,
    pmid VARCHAR,
    n_atlas_cohorts BIGINT,
    total_sample_size BIGINT,
    tissues VARCHAR,
    cohort_descriptions VARCHAR,
    platforms VARCHAR,
    ancestries VARCHAR,
    atlas_traits VARCHAR,
    FOREIGN KEY (study_id) REFERENCES study(study_id)
);

CREATE OR REPLACE VIEW v_study_atlas_enrichment AS
SELECT
    study_id,
    join_method,
    atlas_study_ids,
    pmid,
    n_atlas_cohorts,
    total_sample_size,
    tissues,
    cohort_descriptions,
    platforms,
    ancestries,
    atlas_traits
FROM study_atlas_enrichment
WHERE join_method <> 'none';
