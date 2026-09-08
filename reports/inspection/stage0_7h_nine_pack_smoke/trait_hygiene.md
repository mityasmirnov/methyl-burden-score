# Milestone 10c — trait expansion hygiene

Policy + census. Nine-pack samples: **34234**.

BMI/ancestry labels are now **joined** into
`sample_phenotype_table_hub_nine_pack_v1.parquet` with real masks; heads are
wired for DeepRVAT freeze-and-reuse (GPU not queued). Brain region labels are
honest; **no `brain_head`**. Blood still deferred.

## Locked policy

1. **Disease/cancer case≠control:** use Hub `sample_type` mapped through
   `mbs.datahub_census.SAMPLE_TYPE_CASE_CONTROL`. Train masked heads only on
   rows with `label_status ∈ {case, control}`. Pack-membership
   `disease_mask`/`cancer_mask` false means **unknown / not in pack**, **not**
   control.
2. **BMI:** real numeric Hub labels (`bmi` + `bmi_mask`); regression head
   ready for freeze-encoder reuse. GPU gated on M10 review.
3. **Ancestry:** real `race` labels collapsed rare→`other` (`min_n=10`);
   CE head ready. No class currently meets ≥1k-per-arm GPU bar.
4. **Brain:** region labels joined (`brain_label` + `brain_mask`); **defer
   head** (catalogue, 100% `sample_type=control`).
5. **Blood:** **defer** head (sparse `cell_component`; pack ≠ trait).
6. **Tissue:** do **not** expand the tissue CE head — only `whole blood` ≥1k.

## Disease / cancer census (sample_type → label_status)

| Family | pack samples | case | control | adjacent_normal | ≥1k case+control? | nine-pack mask_true |
|--------|-------------:|-----:|--------:|----------------:|:-----------------:|--------------------:|
| disease | 12218 | 5264 | 6930 | 24 | yes | 12218 |
| cancer | 10101 | 7157 | 1920 | 1024 | yes | 10101 |

## BMI / ancestry / brain (post label-prep)

| Family | Hub pack n | nine-pack mask_true | Head | Note |
|--------|-----------:|--------------------:|------|------|
| bmi | 2070 | 2070 | `bmi_head` wired | continuous; freeze-reuse |
| ancestry | 1380 | 1380 | `ancestry_head` wired | 14 classes after collapse; all &lt;1k |
| brain | 1997 | 1997 | **defer_head** | region catalogue |
| blood | 3402 | 0 | **defer_head** | mask still false |

## Tissue collapse

- labeled=13457; labels=64
- ≥1k: `{'whole blood': 2070}`
- decision: **do_not_expand_tissue_ce**

## Next

- After **10e** encoder: freeze-reuse probes on RBS/MBS for ATS, BMI,
  ancestry, cancer types, and individual diseases (AD 945 / PD 333 /
  stroke 204 cases — census matched controls before GPU).
- BMI/ancestry head-only configs stubbed — **do not auto-queue** until
  ≥1k policy exception or collapse ADR + GPU free.
- Full-catalog sample overview: `sample_overview_hub_geo_v1` +
  `docs/plans/label-prep-bmi-ancestry-sample-overview.md`.
