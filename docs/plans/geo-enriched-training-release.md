# GEO-enriched training release (design)

**Status:** design locked; **not built** (no matrix/phenotype release artifact yet)  
**Parent:** [`geo-metadata-backfill-pre-scale.md`](geo-metadata-backfill-pre-scale.md)  
**Depends on:** audited GEO batch (≥50 GSE), trait eligibility by study, 7G′
architecture decisions for when GEO enters training

## Scope and acceptance

**Goal:** After GEO QC, derive a **new** phenotype/matrix development release that
can include GEO-enriched EWAS_db samples — without mutating frozen ATS
(`matrix-hub-age-tissue-sex-full-v1`) or Hub nine-pack packs.

**Done when:**

1. Release id and artifact layout documented (this file) and referenced from
   [`DATA_CONTRACT.md`](../DATA_CONTRACT.md) / data-infra plan.
2. Eligibility report exists by trait × study for `geo_metadata_backfill` (+ Hub).
3. Four comparison arms are defined in configs (or experiment YAML stubs):

   | Arm | Membership | Role |
   |-----|------------|------|
   | Hub-only | Existing Hub pack samples | Baseline |
   | GEO-only | EWAS_db GSM with ≥1 GEO phenotype, no Hub membership | Coverage / leakage |
   | Hub + GEO | Union of eligible labels | Combined development |
   | Metadata-only confounding | Study/platform features only (7E′ style) | Alarm |

4. No silent overwrite of `deepmat-data-v1` frozen pack semantics; new release id
   e.g. `deepmat-data-geo-dev-v1`.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Mutate ATS? | **No** | Frozen training SoT; GEO is optional enrichment |
| Encoder features from GEO SOFT? | **No** | Phenotypes only; batch/GPL/study-id stay out |
| Disease/cancer from GEO | Only explicit case/control rows; unknown ≠ control | ADR / DATA_CONTRACT |
| When to train | After 7G′ Stage A/B gates + eligibility cutoffs | Ordering in pre-scale plan |
| Release id | `deepmat-data-geo-dev-v1` (provisional) | Distinct from `deepmat-data-v1` |

## Artifact sketch

```text
$MBS_DATA_ROOT/canonical/releases/deepmat-data-geo-dev-v1/
  manifest.json
  catalog/           # optional view or symlink policy TBD at build time
  phenotypes/
    sample_info_hub_only.parquet
    sample_info_geo_only.parquet
    sample_info_hub_geo.parquet
  matrices/          # only if convert paths are approved; else phenotype-only first
```

**Do not build until** repaired-pilot audit is accepted, an expanded batch is
QC'd, and trait eligibility shows acceptable labels for the intended heads.
**Do not** use GEO disease/cancer while cases or controls are missing for
eligibility. **Do not** mutate `matrix-hub-age-tissue-sex-full-v1`.

## Non-goals

- Full EWAS_db crawl as a training gate
- ComBat / batch as encoder inputs
- Treating Atlas cohort fields as sample labels
- Milestone **7** OOF before 7G′ Stage A/B

## Open (resolve at build time)

- Exact eligibility cutoffs for GEO tissue/age aux heads
- Whether GEO-only matrix convert uses EWAS_db beta txt paths only
- Seed-gene selection: Hub-only vs Hub+GEO (default Hub-only until audited)
