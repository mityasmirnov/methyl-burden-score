# Data and infrastructure improvements (Stage 0)

Binding priorities for **better scientific results** from the existing MBS
pipeline — not optional product layers (ComBat, ClickHouse, TileDB) deferred
per [`TODO_PIPELINE.md`](../TODO_PIPELINE.md).

**Snapshot:** 2026-09-02. Catalog `deepmat-data-v1`: **149 244** samples,
**1 584** studies, **132 289** EWAS_db GSM files (**1 353**/1 989 studies;
`mirror_complete=false`). Hub nine-pack matrices and ATS union (**13 548** GSM)
are training-ready. Gate: **7G′** gene-only architecture selection.

## 1. Download reliability (EWAS_db)

**Problem:** ~79% through the study list; **4 224** GSM files still missing
after logged `wget` failures; **~1 474** bogus `(.+?)` filenames from HTML href
parser leakage ([`ewas_db_download_failures.md`](../../reports/inspection/deepmat_data_v1/ewas_db_download_failures.md)).

| Action | Impact | Effort |
|--------|--------|--------|
| Filter `list_hrefs` to drop regex artifacts and non-`GSM*.txt` names | Cuts false failures; faster retries | Small |
| Port disease-pack resilient `wget` flags (`--retry-connrefused`, longer timeouts) into `mirror_ewas_db` | Fewer transient CNCB drops | Small |
| Run `make retry-ewas-db-failures` after mirror pass | Recovers real missing GSM without full re-crawl | Ops |
| Post-download hook (done): summarize + `catalog refresh-release` | Catalog stays aligned with disk | Done |
| Optional: parallel wget per study (cap concurrency) | Faster mirror without hammering CNCB | Medium |

EWAS_db is **not** a 7A/7B/7G gate ([ADR 0007](../adr/0007-crossfit-prerequisites.md))
but improves external-test diversity and EWAS_db-only convert paths.

## 2. Catalog and phenotype SoT

**Problem:** Census reports can lag one refresh behind ingest; donor/replicate
columns empty; EWAS_db-only samples lack Hub phenotype rows.

| Action | Impact | Effort |
|--------|--------|--------|
| Run `make catalog-refresh-release` after large ingest (hook now calls subset) | Consistent census + eligibility | Done (partial) |
| Extend post hook to `validate-release` + `phenotype-census` + `trait-eligibility` | Single command = full 7A reports | Small |
| Ingest EWAS_db sample metadata where available (Atlas study/cohort joins) | Labels for ~113k EWAS_db-only GSM | Medium |
| Registry `sample_count` from matrix sample indexes | Honest N in `phenotype_registry.yaml` | Small |
| Normalize `450K` → `HM450` on catalog refresh | Fewer platform string splits | Small |

**Use now:** `trait_eligibility` + `v_sample_pack_overlap` to restrict heads
(disease/cancer need cases+controls; tissue needs class support). Do not treat
missing disease/cancer as controls ([`DATA_CONTRACT.md`](../DATA_CONTRACT.md)).

## 3. Training matrices and compute

**What exists:**

- ATS union: `matrix-hub-age-tissue-sex-full-v1` (13 548 × 482k, HM450)
- Nine-pack virtual: `matrix-hub-nine-pack-virtual-v1` (34 234 routed GSM)
- Graph v2: `graph-grch38-gencode38-cgi-tile-v2` (RBS/TBS topology for 7F/7G)

| Action | Impact | Effort |
|--------|--------|--------|
| **7G′ Stage A** on matched `gene_cols` (Cascade vs Flat vs C-mvalue-G) | Picks encoder before OOF | In progress |
| **7G′ Stage B** fold-selected `C-mvalue-enetS` + `direct_cpg.zarr` | Full multipath scores for Milestone 7 | Pending |
| Train on GPU (`--device cuda`); dense betas prefix in RAM for cascade | 10–50× wall-clock vs CPU | Ops |
| Keep `cv_budget.max_loci` explicit per arm (65k ATS vs 8k smoke) | Fair architecture comparison | Discipline |
| Study-grouped splits only; no sample-ID features in encoder | Prevents leakage inflation | Invariant |

Do **not** expand to full 482k loci dense union across nine packs (~61 GB+ per
tensor) — virtual index + per-pack Zarr is the correct scale path.

## 4. Evaluation and reporting hygiene

| Action | Impact | Effort |
|--------|--------|--------|
| Comparable ranking (`mbs_e2e`, matched probes) — 7G′ | Apples-to-apples neural vs classical | In progress |
| Arm glossary + inspection JSON slimming | Auditable reports without multi-MB ROC blobs | Done |
| Validation vs external tissue F1 separation in tissue probe | Stops optimistic checkpoint selection | Done (7G) |
| Metadata-only baseline as leakage alarm only (7E′) | Guards against study/platform shortcuts | Done |

Before Milestone **7** OOF: score-orientation anchor ([ADR 0008](../adr/0008-score-identifiability.md)),
fold-fitted Level-1 normalization (7D), and locked panel from 7G′ Stage B.

## 5. Scientific leverage (data → better models)

Ordered by expected lift on current gate:

1. **Tissue ontology + class weights** — 72 raw labels → coarse head + conditional fine (brain/blood packs); already wired in cascade loss.
2. **Masked multitask** — age/tissue/sex primary; disease/cancer only where eligibility passes; unknown stays unknown.
3. **Late fusion solvers** (`balanced_logistic`, PCA) — cheap lift when MBS features are linearly separable for tissue (7G tissue probe evidence).
4. **Graph-v2 multipath** — gene-aggregated RBS + orphan RBS + direct CpGs (7F topology); compare gene-only vs full paths in 7G′.
5. **EWAS Atlas** — external validation of association direction, not training labels ([ADR 0002](../adr/0002-ewas-datahub-primary-source.md)).
6. **CpGPT static features** — optional locus prior; manifest-gated; not default until ablation shows gain.

## 6. Deferred (explicit non-goals here)

- ComBat / batch harmonization across studies (§8 programme)
- ClickHouse / TileDB migration ([ADR 0005](../adr/0005-catalog-matrix-independence.md))
- Final 5×6 OOF before 7G′ completes ([ADR 0007](../adr/0007-crossfit-prerequisites.md))
- Retraining frozen v0.1 references

## Commands (maintenance)

```bash
source scripts/activate_data_environment.sh

# Full catalog + reports after ingest
make catalog-refresh-release

# EWAS_db failure audit + retry
make summarize-ewas-db-failures
make retry-ewas-db-failures

# Resume mirror (auto post-hook on exit)
nohup bash scripts/download_ewas_datahub.sh EWAS_db \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_datahub_EWAS_db.log" 2>&1 &
```

## Acceptance

Track in [`TODO_PIPELINE.md`](../TODO_PIPELINE.md) ops notes; no milestone
status change from this plan alone. Individual items graduate to ADRs or
milestone plans when they become binding implementation work.
