-- Reusable catalog views for source inspection and leakage checks.

CREATE OR REPLACE VIEW v_source_inventory AS
SELECT
    sr.source_name,
    sr.source_version,
    s.study_id,
    s.gse_id,
    s.platform_id,
    s.processing_level,
    count(DISTINCT sm.sample_id) AS n_samples,
    count(DISTINCT af.assay_file_id) AS n_assay_files,
    sum(coalesce(af.byte_size, 0)) AS total_bytes
FROM study AS s
LEFT JOIN source_release AS sr USING (source_release_id)
LEFT JOIN sample AS sm USING (study_id)
LEFT JOIN assay_file AS af USING (study_id)
GROUP BY ALL;

CREATE OR REPLACE VIEW v_study_platform_tissue AS
SELECT
    s.study_id,
    s.gse_id,
    s.platform_id,
    s.processing_level,
    sm.tissue_ontology_id,
    sm.tissue_raw,
    count(*) AS n_samples,
    count(sm.age) AS n_age,
    min(sm.age) AS age_min,
    max(sm.age) AS age_max,
    avg(sm.age) AS age_mean
FROM sample AS sm
JOIN study AS s USING (study_id)
GROUP BY ALL;

CREATE OR REPLACE VIEW v_duplicate_source_samples AS
SELECT
    source_sample_id,
    count(*) AS n_rows,
    count(DISTINCT study_id) AS n_studies,
    list(sample_id ORDER BY sample_id) AS sample_ids
FROM sample
WHERE source_sample_id IS NOT NULL
GROUP BY source_sample_id
HAVING count(*) > 1;

CREATE OR REPLACE VIEW v_replicate_groups AS
SELECT
    replicate_group,
    count(*) AS n_samples,
    count(DISTINCT donor_id) AS n_donors,
    count(DISTINCT study_id) AS n_studies,
    list(sample_id ORDER BY sample_id) AS sample_ids
FROM sample
WHERE replicate_group IS NOT NULL
GROUP BY replicate_group;

CREATE OR REPLACE VIEW v_platform_locus_overlap AS
SELECT
    p.platform_id,
    count(DISTINCT p.probe_id) AS n_probes,
    count(DISTINCT ple.locus_id) AS n_mapped_loci,
    count(*) FILTER (WHERE ple.is_primary) AS n_primary_edges
FROM probe AS p
LEFT JOIN probe_locus_edge AS ple
    ON p.probe_id = ple.probe_id
   AND p.platform_id = ple.platform_id
GROUP BY p.platform_id;

CREATE OR REPLACE VIEW v_gene_platform_coverage AS
SELECT
    r.gene_id,
    r.region_type,
    pr.platform_id,
    count(DISTINCT lre.locus_id) AS n_loci
FROM region AS r
JOIN locus_region_edge AS lre USING (region_id)
JOIN probe_locus_edge AS ple USING (locus_id)
JOIN probe AS pr
    ON pr.probe_id = ple.probe_id
   AND pr.platform_id = ple.platform_id
GROUP BY ALL;

CREATE OR REPLACE VIEW v_fold_balance AS
SELECT
    fa.split_id,
    fa.outer_fold,
    fa.role,
    s.study_id,
    st.platform_id,
    s.tissue_ontology_id,
    count(*) AS n_samples,
    count(s.age) AS n_age,
    min(s.age) AS age_min,
    max(s.age) AS age_max
FROM fold_assignment AS fa
JOIN sample AS s USING (sample_id)
JOIN study AS st USING (study_id)
GROUP BY ALL;

CREATE OR REPLACE VIEW v_artifact_lineage AS
SELECT
    e.experiment_id,
    e.git_commit,
    e.data_release_id,
    e.split_id,
    e.status,
    e.graph_artifact_id,
    e.feature_artifact_ids,
    e.created_at
FROM experiment AS e;
