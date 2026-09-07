-- Hub baseline vs EWAS_db lane flags (membership census; not training labels).

CREATE TABLE IF NOT EXISTS sample_lane_flags (
    sample_id VARCHAR PRIMARY KEY,
    study_id VARCHAR,
    in_hub_baseline BOOLEAN NOT NULL,
    in_ewas_db BOOLEAN NOT NULL,
    hub_families VARCHAR,
    FOREIGN KEY (sample_id) REFERENCES sample(sample_id)
);

CREATE TABLE IF NOT EXISTS study_lane_flags (
    study_id VARCHAR PRIMARY KEY,
    in_hub_baseline BOOLEAN NOT NULL,
    in_ewas_db BOOLEAN NOT NULL,
    n_hub_samples BIGINT,
    n_ewas_db_samples BIGINT,
    n_samples BIGINT,
    FOREIGN KEY (study_id) REFERENCES study(study_id)
);

CREATE OR REPLACE VIEW v_sample_lane_flags AS
SELECT
    sample_id,
    study_id,
    in_hub_baseline,
    in_ewas_db,
    hub_families,
    CASE
        WHEN in_hub_baseline AND in_ewas_db THEN 'both'
        WHEN in_hub_baseline THEN 'hub_only'
        WHEN in_ewas_db THEN 'ewas_db_only'
        ELSE 'neither'
    END AS lane_class
FROM sample_lane_flags;

CREATE OR REPLACE VIEW v_study_lane_flags AS
SELECT
    study_id,
    in_hub_baseline,
    in_ewas_db,
    n_hub_samples,
    n_ewas_db_samples,
    n_samples,
    CASE
        WHEN in_hub_baseline AND in_ewas_db THEN 'both'
        WHEN in_hub_baseline THEN 'hub_only'
        WHEN in_ewas_db THEN 'ewas_db_only'
        ELSE 'neither'
    END AS lane_class
FROM study_lane_flags;
