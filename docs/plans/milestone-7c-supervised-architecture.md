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
loader) are closed. **7B** Hub full matrices are done.

### Hub join (after 7B)

| Item | Status | Notes |
|------|--------|-------|
| Disease/cancer long-form join | **done** (fixtures + Hub sidecar verify) | `hub_longform_ready` on `matrix-hub-*-full-v1`; multi-hot ≠ control |
| Optional Hub train smoke | **done** | `stage0_flat_hub_disease_multilabel_smoke.yaml`; report `reports/inspection/stage0_7c_hub_disease_smoke/` |
| AUROC / AUPRC / ECE in trainer JSON | **done** | Holdout emits binary sex/tissue AUROC/ECE; multilabel `disease_auroc` when both classes observed under mask |

Do **not** treat `matrix-hub-disease-from-agepack-v1` (unique-GSM multiclass) as
the multi-label solution.

### Ops / later topology

| Item | Status | Notes |
|------|--------|-------|
| `apply_orientation` on real gene-mean M | **done** (flat + hier train) | Train-fold gene-mean MBS vs signed gene-mean M; flips head weights + rewrites ckpts; `score_manifest.json` |
| Full-genome graph-v2 artifact | **on disk** | `$MBS_DATA_ROOT/canonical/graphs/graph-grch38-gencode38-cgi-tile-v2/`; `reports/inspection/annotation_graph_cgi_tile_v2/`; plan [`milestone-7c-graph-v2-topology.md`](milestone-7c-graph-v2-topology.md) |
| Multi-system hier index (RBS/TBS) | **done** | `region_systems` on `build_locus_region_gene_index`; default gene-only |
| Branch arms `rbs`/`tbs` | **train-time masks done** | System filter at index build; Hub-scale 7E still pending |

## Non-goals

v0.1 overwrite; 7D MAD; 7E 3×2; dual hyper/hypo channels; LOEUF/constraint scores;
SCREEN cCRE download.
