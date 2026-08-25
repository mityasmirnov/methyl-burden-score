# Milestone 7D — Hub Level-1 A/B smoke

- **Created:** 2026-08-25T10:46:53.526531+00:00
- **Matrix:** `matrix-hub-age-tissue-sex-full-v1`
- **Cohort:** deeprvat_hub age/tissue/sex
- **Device:** `cuda:0`

## Runs

| Channel | run_id | robust_deviation | fold_norm | input_dim | n_train |
|---------|--------|------------------|-----------|-----------|---------|
| A | `stage0-7d-level1-a` | False | False | 131 | 9794 |
| B | `stage0-7d-level1-b` | True | True | 133 | 9794 |

## Level-1 (channel B)

- n_train_samples: `9794`
- n_estimated: `512`
- n_unestimated: `0`
- mu_sha256: `1a40231ba13eee0e56c5bc9578e5de6afe502d1627094e9354c5d8025b729b98`
- sigma_sha256: `6021ca818f327c6b52cb31a17fe54dcf747dbac33574e2bbe24421707985d5d3`
- locus_ids_sha256: `6be6d4b54c4fa574f19d1544249d9a8f86c974c693d4b363ecca2c1db270902c`
- sigma_min: `1e-06`
- formula: `z=(M-median)/max(1.4826*MAD,sigma_min); train-fold only`

## Acceptance

- Channel A has no `fold_norm/`: **True**
- Channel B has schema-valid `fold_norm/`: **True**
- Identical study-grouped split IDs: **True**
- input_dim B > A (z + norm_present channels): **True**
- GMQN canonical betas untouched (no writes under matrices): **True**

## Non-claims

- Short smoke (`max_loci` + few epochs); not phenotype SOTA.
- Comprehensive RBS/TBS + graph-v2 remain 7E prerequisites, not 7D.
