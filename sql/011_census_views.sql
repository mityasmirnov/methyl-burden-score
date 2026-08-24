-- Phenotype census and eligibility views (Milestone 7A).

CREATE OR REPLACE VIEW v_sample_pack_overlap AS
SELECT
    sample_id,
    count(DISTINCT phenotype_family) AS n_families,
    list(DISTINCT phenotype_family ORDER BY phenotype_family) AS families
FROM sample_source_membership
GROUP BY sample_id;

CREATE OR REPLACE VIEW v_sample_label_conflicts AS
SELECT
    sample_id,
    phenotype_id,
    count(DISTINCT coalesce(cast(numeric_value AS VARCHAR), categorical_value)) AS n_distinct_values,
    list(DISTINCT source_family ORDER BY source_family) AS source_families,
    list(DISTINCT coalesce(cast(numeric_value AS VARCHAR), categorical_value)) AS values
FROM sample_phenotype
WHERE is_observed
GROUP BY sample_id, phenotype_id
HAVING count(DISTINCT coalesce(cast(numeric_value AS VARCHAR), categorical_value)) > 1;

CREATE OR REPLACE VIEW v_phenotype_prevalence AS
SELECT
    phenotype_id,
    source_family,
    count(*) AS n_rows,
    count(*) FILTER (WHERE is_observed) AS n_observed,
    count(DISTINCT sample_id) AS n_samples,
    count(DISTINCT categorical_value) FILTER (
        WHERE categorical_value IS NOT NULL
    ) AS n_categorical_levels
FROM sample_phenotype
GROUP BY phenotype_id, source_family;

CREATE OR REPLACE VIEW v_phenotype_study_support AS
SELECT
    sp.phenotype_id,
    sp.source_family,
    s.study_id,
    count(DISTINCT sp.sample_id) AS n_samples
FROM sample_phenotype AS sp
JOIN sample AS s USING (sample_id)
WHERE sp.is_observed
GROUP BY ALL;

CREATE OR REPLACE VIEW v_phenotype_platform_support AS
SELECT
    sp.phenotype_id,
    sp.source_family,
    coalesce(st.platform_id, 'unknown') AS platform_id,
    count(DISTINCT sp.sample_id) AS n_samples
FROM sample_phenotype AS sp
JOIN sample AS sm USING (sample_id)
LEFT JOIN study AS st USING (study_id)
WHERE sp.is_observed
GROUP BY ALL;

CREATE OR REPLACE VIEW v_tissue_class_distribution AS
SELECT
    categorical_value AS tissue_label,
    source_family,
    count(DISTINCT sample_id) AS n_samples
FROM sample_phenotype
WHERE phenotype_id = 'tissue'
  AND is_observed
  AND categorical_value IS NOT NULL
GROUP BY ALL
ORDER BY n_samples DESC;

CREATE OR REPLACE VIEW v_disease_case_control_distribution AS
SELECT
    sp.source_family,
    sp.label_status,
    count(DISTINCT sp.sample_id) AS n_samples,
    count(DISTINCT sm.study_id) AS n_studies
FROM sample_phenotype AS sp
JOIN sample AS sm USING (sample_id)
WHERE sp.phenotype_id IN ('disease', 'cancer')
GROUP BY ALL;

CREATE OR REPLACE VIEW v_age_distribution_by_study AS
SELECT
    sm.study_id,
    count(DISTINCT sp.sample_id) AS n_samples,
    min(sp.numeric_value) AS age_min,
    max(sp.numeric_value) AS age_max,
    avg(sp.numeric_value) AS age_mean
FROM sample_phenotype AS sp
JOIN sample AS sm USING (sample_id)
WHERE sp.phenotype_id = 'age'
  AND sp.is_observed
  AND sp.numeric_value IS NOT NULL
GROUP BY sm.study_id;

CREATE OR REPLACE VIEW v_trait_missingness AS
SELECT
    phenotype_id,
    source_family,
    count(*) AS n_rows,
    count(*) FILTER (WHERE NOT is_observed OR label_status = 'unknown') AS n_missing_or_unknown,
    count(*) FILTER (WHERE is_observed AND label_status = 'observed') AS n_observed
FROM sample_phenotype
GROUP BY phenotype_id, source_family;

CREATE OR REPLACE VIEW v_trait_training_eligibility AS
SELECT
    phenotype_id,
    phenotype_family,
    task_type,
    n_samples,
    n_cases,
    n_controls,
    n_unknown,
    prevalence,
    n_studies,
    n_platforms,
    eligible_core_task,
    eligible_auxiliary_task,
    eligible_external_evaluation,
    exclusion_reason
FROM trait_eligibility
ORDER BY phenotype_family, phenotype_id;

CREATE OR REPLACE VIEW v_split_trait_balance AS
SELECT
    fa.split_id,
    fa.role,
    sp.phenotype_id,
    count(DISTINCT fa.sample_id) AS n_samples,
    count(DISTINCT sm.study_id) AS n_studies
FROM fold_assignment AS fa
JOIN sample AS sm USING (sample_id)
LEFT JOIN sample_phenotype AS sp
    ON sp.sample_id = fa.sample_id
   AND sp.is_observed
GROUP BY ALL;

CREATE OR REPLACE VIEW v_ewas_db_ingest_status AS
SELECT
    sr.source_release_id,
    sr.source_system,
    count(DISTINCT s.study_id) AS n_studies,
    count(DISTINCT sm.sample_id) AS n_samples,
    count(DISTINCT af.assay_file_id) AS n_assay_files,
    sum(coalesce(af.byte_size, 0)) AS total_bytes
FROM source_release AS sr
LEFT JOIN study AS s USING (source_release_id)
LEFT JOIN sample AS sm USING (study_id)
LEFT JOIN assay_file AS af USING (study_id)
WHERE sr.source_system = 'ewas_datahub_db'
GROUP BY ALL;
