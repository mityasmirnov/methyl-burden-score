# GEO-enriched training release

**Status:** **phenotype release built** (2026-09-07) — matrices not converted;
training launch still gated  
**Parent:** [`geo-metadata-backfill-pre-scale.md`](geo-metadata-backfill-pre-scale.md)  
**Depends on:** audited GEO batch (≥50 GSE), trait eligibility by study,
Milestone **10/11** architecture decisions for when GEO enters training

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
| When to train | After Milestone **10/11** gates + eligibility cutoffs | Ordering in pre-scale plan |
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

**Built (phenotype-only):** `scripts/build_geo_dev_release.py` →
`$MBS_DATA_ROOT/canonical/releases/deepmat-data-geo-dev-v1/` +
`reports/inspection/deepmat_data_geo_dev_v1/summary.md` +
`configs/experiment/geo_dev/*.yaml` stubs.

| Arm (2026-09-07) | n_samples |
|------------------|---------:|
| `hub_only` | 34 090 |
| `geo_only` | 34 990 |
| `hub_geo` | 69 080 |
| `metadata_only` | 69 080 |

**Do not** use GEO disease/cancer for training while per-study case/control is
unbalanced (global `trait_eligibility` can look green — see
`geo_backfill_batch/eligibility_by_study.md`). **Do not** mutate
`matrix-hub-age-tissue-sex-full-v1`. **Do not** launch age/tissue training or
seed-gene selection on GEO until Milestone **10/11** gates allow.

## Non-goals

- Full EWAS_db crawl / larger GEO crawl as a training gate
- ComBat / batch as encoder inputs
- Treating Atlas cohort fields as sample labels
- Milestone **12** OOF before architecture finalists lock
- Wiring GEO into live `mbs train` configs in this release

## Open (resolve before training launch)

- Exact eligibility cutoffs for GEO tissue/age aux heads
- Whether GEO-only matrix convert uses EWAS_db beta txt paths only
- Seed-gene selection: Hub-only vs Hub+GEO (default Hub-only until audited)
