# Milestone 7C: Supervised architecture corrections

Status: **done** (fixture acceptance under
`reports/inspection/stage0_7c_architecture/` + `tests/unit/test_stage0_7c.py`).
Normative ADRs: [0006](../adr/0006-multipath-noncoding-scores.md),
[0008](../adr/0008-score-identifiability.md).
Checklist: [`TODO_PIPELINE.md`](../TODO_PIPELINE.md).
Programme context: [`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).

## Scope and acceptance

Phase A (age/tissue/sex; does not wait on 7B): trainer shuffle + token-budget
sampler; centered heads; constraint-aware splits; honest donor/replicate;
full metric suite helpers; controls; `static_present`; polarity anchor schema;
matched encoder.

Phase B: graph builder id `graph-grch38-gencode38-cgi-tile-v2` (CGI RBS +
CpG-count tiles); elastic-net direct branch; independently tagged
`gene|rbs|tbs|direct` arms on fixtures; masked disease/cancer heads on
fixtures (`hub_longform_ready` + sidecar gates Hub join).

**Done when (met):** unit tests + fixture inspection report.
Do not retrain v0.1. Hub-scale 7E CV is out of scope.

## Locked decisions

| Choice | Decision |
|--------|----------|
| Orientation | Polarity flip vs signed gene-mean M-value; persist `score_polarity` |
| RBS | UCSC CGI/shore on gene-unassigned loci (no SCREEN cCRE required) |
| TBS | Greedy 50-CpG tiles on leftover mapped loci |
| Direct \(z\) | Train-fold centered M until 7D MAD |
| Donor | `null` when unknown; never impersonate `study_id` |
| Metadata-only | Separate linear model; IDs never enter the encoder |
| Branch ablations | Independent trains on identical fixture folds |

## Schemas / contracts

- [`schemas/score_manifest.schema.json`](../../schemas/score_manifest.schema.json)
- Graph `region_policy.region_systems`: `gene`, `rbs`, `tbs`

## Remaining follow-ons (not blocking 7C Done when)

Fixture acceptance and residual code polish (orientation train-path + long-form
loader) are closed. Split leftovers by whether they need **7B Hub matrices**.

### Blocked on 7B (convert still running)

| Item | Status | Notes |
|------|--------|-------|
| Hub disease multilabel smoke | Waiting on `matrix-hub-disease-full-v1` | Config ready: `configs/experiment/stage0_flat_hub_disease_multilabel.yaml`. Needs `sample_index.parquet` + `sample_phenotypes.parquet` (`hub_longform_ready`). |
| Hub cancer multilabel smoke | Waiting on `matrix-hub-cancer-full-v1` | Same join path; cancer convert in progress under `$MBS_DATA_ROOT/canonical/matrices/`. |
| AUROC / AUPRC / ECE in trainer JSON | Helpers in `metrics.py`; not emitted from holdout path | Prefer wiring when Hub disease/cancer (or binary sex) eval runs |

Do **not** treat `matrix-hub-disease-from-agepack-v1` (unique-GSM multiclass) as
the multi-label solution.

### Not blocked on 7B

| Item | Status | Notes |
|------|--------|-------|
| `apply_orientation` on real gene-mean M | **done** (flat + hier train) | Train-fold gene-mean MBS vs signed gene-mean M; flips head weights + rewrites ckpts; `score_manifest.json` |
| Disease/cancer long-form join (code) | **done** (fixtures) | `load_longform_multilabel` + masked BCE; Hub smoke above |
| Full-genome graph-v2 artifact | Builder + fixture path done; do not rebuild genome in CI | Run `mbs graph build --graph-id graph-grch38-gencode38-cgi-tile-v2` under `$MBS_*` when needed |
| Multi-system hier index (RBS/TBS DeepSet combo) | `locus_region_gene` still filters `region_system==gene` for v0.1-compatible hier | Shared CpG encoder + per-system region embeddings is 7E topology work |
| Branch arms `rbs`/`tbs` | CLI/run-dir + fixture overfit; still gene FlatDeepSet features | True region-system masks need graph-v2 train index |

## Non-goals

v0.1 overwrite; 7D MAD; 7E 3×2; dual hyper/hypo channels; LOEUF/constraint scores;
SCREEN cCRE download.
