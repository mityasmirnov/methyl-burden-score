# Milestone 7E′: Hub multitask + analysis hygiene

Status: **pending** (extra steps beside **7E**; required before Milestone **7**).
Does not reopen 7A–7D or the frozen ATS matrix.
Checklist: [`TODO_PIPELINE.md`](../TODO_PIPELINE.md).
7E bake-off: [`milestone-7e-development-cv.md`](milestone-7e-development-cv.md).
Graph-v2 (done): [`milestone-7c-graph-v2-topology.md`](milestone-7c-graph-v2-topology.md).

## What ATS is

**ATS** = **A**ge / **T**issue / **S**ex. It is the frozen Hub **GSM-union** of
the three baseline packs:

```text
matrix-hub-age-tissue-sex-full-v1   13 548 samples × 482 379 loci
freeze name: deepmat-data-age-tissue-sex-v1
split: hub-age-tissue-sex-full-auto-v1 (9489 / 2074 / 1985)
```

It is **not** “all downloaded methylation.” Unique Hub GSM across nine packs is
**34 234**. EWAS_db listing is **92 971** GSM files and is a different lane.

## Can ATS grow because EWAS_db downloaded more?

**No.** More `EWAS_db/{GSE}/` dirs do not enlarge the Hub age/tissue/sex zips.
Those packs are already complete. Do **not** overwrite the ATS freeze.

**Yes, a larger training cohort already exists on disk:** union (or virtual
multi-store) of the nine Hub full matrices, with masked heads for every
available trait — including disease and cancer. That is this milestone, not a
bigger ATS freeze.

EPIC/850K is also not in ATS: every Hub `*_sample_info.parquet` row is
`platform=450K`. EPIC lives in the annotation graph and in EWAS_db
`add_txt_850`, not in the nine baseline zips.

## Disease/cancer “core” vs train-on-everything

Eligibility `eligible_core_task=False` for disease/cancer means: **do not treat
missing labels as controls** and **do not claim case–control epidemiology**
until documented controls exist (ADR / DATA_CONTRACT). It does **not** mean
skip those heads.

**Locked here:** train age, tissue, sex, **and** disease (cancer too) with
**masked** loss (`unknown ≠ 0`). Pack-wide blood `cell_component` stays unused
(~1.1% populated).

## Scope and acceptance

| Deliverable | Done when |
|-------------|-----------|
| Hub multitask | Age + tissue + sex + disease (+ cancer) heads on Hub packs already converted; masked unknown≠control; new split **not** overwriting ATS freeze |
| Metadata-only | Study/platform/tissue → phenotype on the **same** 7E folds; confounding ceiling in the 7E report |
| Catalog string | 7B `platform_id=450K` aliases to `HM450` on next `refresh-release` (probe map already HM450) |
| Census follow-ons | Donor/replicate IDs when present; within-study age/BMI ranges |
| Tests | Census/eligibility CLI tests write a **temp** `--report-dir` (no hyphen-path clobber) |
| Git | `*.RData` ignored; no Hub sample blobs under `reports/inspection/ewas_datahub_samples/` |
| Blood | Config/docs: do not use `cell_component` as a pack-wide head |
| Matched encoders | 7E YAML: same width/GELU/dropout/LN for flat vs hier (also 7E Done when) |

Graph-v2 independent RBS/TBS arms are **7E**, not this list (artifact is on disk).

## Non-goals

- Overwriting ATS / v0.1; waiting on EWAS_db completeness; converting EPIC
  `add_txt_850` into Hub-scale matrices; treating unlabeled disease rows as
  controls.

## Open questions

None blocking 7E ATS bake-off. Hub-union matrix layout (physical concat vs
virtual multi-store) can follow the 7B index.
