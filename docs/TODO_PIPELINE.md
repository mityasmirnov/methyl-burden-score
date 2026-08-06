# Stage 0 scientific pipeline TODO

Authoritative milestone list for coding agents. Update status here when a
milestone is **truly done** (acceptance criteria met), not when scaffolding
exists. Cursor agents must read this file at session start and after finishing
work; see `.cursor/rules/pipeline-todo.mdc`.

Status values: `done` | `in_progress` | `pending` | `deferred`

True next milestone after bootstrap:

> one source ingested cleanly → one graph built → one canonical matrix written →
> one baseline trained → phenotype registry + multi-pack eval → Hub metadata
> contracts → real Hub pack matrices + study-grouped eval → **multitask shared
> encoder (5c — start now)** → hierarchical → one cross-fitted score matrix.

**Current gate:** Milestone **5c**. Prerequisites 5b / 5b′ / 5b″ are `done`.
Do **not** wait on incomplete disease/cancer profile zips to begin 5c.

Primary open data source going forward: **EWAS Data Hub**
([ADR 0002](adr/0002-ewas-datahub-primary-source.md);
[strategic plan](STRATEGIC_PLAN.md)). Milestone 1 evidence remains the
historical CpGCorpus inspection and is not re-opened.

---

## Ops notes (plumbing)

### Done

- Project-local `/data` roots (`MBS_*`), DuckDB schema via `mbs catalog init`,
  shallow `mbs inspect source`, CpGCorpus download scripts with log tee +
  `scripts/download_cpgcorpus_background.sh` (nohup; no auto-download).
- Path / catalog / CLI / no-home tests green; vendor kept as gitlink submodules.
- ADR: `docs/adr/0001-project-local-workspace.md`.
- ADR: `docs/adr/0002-ewas-datahub-primary-source.md` (EWAS Data Hub primary).
- Raw download inventory (sizes, schemas, top rows, cleanup checklist):
  `reports/inspection/raw_inventory/summary.md` (+ `summary.json`).
- EWAS Atlas small tables + Hub sample-info structure contracts (**done; required
  reading for 5c**): `mbs inspect ewas-metadata` →
  `reports/inspection/ewas_metadata_structure/`; durable doc
  [`EWAS_METADATA.md`](EWAS_METADATA.md); plan
  [`plans/ewas-metadata-structure.md`](plans/ewas-metadata-structure.md).
  Sample-info export prefers unpacked
  `reports/inspection/ewas_datahub_samples/*.txt` when zips are absent.
- EWAS Atlas + Data Hub download scripts / Makefile targets
  (`docs/EWAS_DATA.md`).

### Blocked / caveats (do not treat as scientific blockers)

- **uv location:** Bootstrap installs/copies `uv` under `$MBS_ROOT/.tools/uv/bin`.
  Older runs may still have a copy in `$HOME/.local/bin`. Caches remain under
  `$MBS_CACHE_ROOT` (correct). Prefer
  `source scripts/activate_data_environment.sh` (puts `.tools` on `PATH` first).
- **CUDA:** Host driver 570.x / CUDA 12.8 (`nvidia-smi`). Main `.venv` pins
  Linux/Windows torch to the PyTorch **cu128** index
  (`torch==2.11.0+cu128` at last sync) via `[[tool.uv.index]]` /
  `[tool.uv.sources]` in `pyproject.toml`. Verified:
  `torch.cuda.is_available()` is true on all three GPUs. Do not switch back to
  PyPI's default Linux wheel (`cu130`) without a driver that supports CUDA 13.0.
  macOS still resolves torch from PyPI (CPU/MPS).
- **CpGPT (optional extra):** `uv sync --all-groups --extra cpgpt` installs the
  vendored pin; `download_cpgpt(model="small", species="human",
  cache_dir=$HF_HOME)` materializes weights under
  `$MBS_CACHE_ROOT/huggingface` (see `docs/STATIC_FEATURES.md`). Not a CI /
  default sync dep. Full `CpGPTInferencer` import currently breaks on the MBS
  torch / torchtune / torchao pin set; offline export uses
  `mbs.static_features.cpgpt_adapter` (dna_encoder weights only) — see
  `docs/STATIC_FEATURES.md` § “CpGPT vs MBS torch pin”.
- **MethylGPT (dedicated env):** cannot share main `.venv` (torchtext requires
  torch ~2.1–2.4). Use `make setup-methylgpt` → `.venv-methylgpt` and
  `make download-methylgpt` → `$MBS_DATA_ROOT/raw/methylgpt/` (medium
  checkpoint + type3 probe IDs). See `docs/STATIC_FEATURES.md`. Ablation-only
  for Stage 0 static features.
- **Intentional non-goals until later milestones:** Parquet population / ingest
  into catalog tables; committing the shallow whole-corpus inventory under
  `reports/inspection/cpgcorpus/` (GSE/GPL report for milestone 1 is enough).
  MethylGPT token-prior export remains ablation-only (not Stage 0 default).
- **Incomplete Hub profile zips (do not block 5c MVP):**
  `disease_methylation_v1.zip` / `cancer_methylation_v1.zip` may still be
  downloading (`wget -c`). Age/tissue/blood/brain packs +
  `matrix-hub-{age,tissue,blood,brain}-studyholdout-v1` are enough to start 5c.
  Sample-info Parquets for disease/cancer labels already exist under
  `$MBS_DATA_ROOT/canonical/phenotypes/` but betas require the profile zips.

### Useful commands (already safe to re-run)

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
uv run mbs doctor --create-directories
uv run mbs catalog init
uv run mbs inspect source --source-id cpgcorpus
uv run mbs inspect cpgcorpus-gpl --gse GSE125367 --gpl GPL21145
uv run mbs inspect ewas-metadata
make download-ewas-study STUDY=GSE35069
uv run mbs matrix convert --study-id GSE35069 --platform-id HM450 --verify
# optional foundation-model tooling:
# uv sync --all-groups --extra cpgpt
# uv run --extra cpgpt mbs features export-cpgpt
# make setup-methylgpt && make download-methylgpt
```

See [`EWAS_METADATA.md`](EWAS_METADATA.md) before multitask / Hub phenotype joins
(Milestone 5c).
### Primary downloads (EWAS Open Platform)

```bash
make download-ewas-datahub
make download-ewas-atlas
make download-ewas-study STUDY=GSE35069
make download-manifests
```

See [`EWAS_DATA.md`](EWAS_DATA.md) for HTTP trees (`EWAS_db/`, `download/`,
Atlas TSVs) and background `nohup` patterns.

### Optional alternate (CpGCorpus; requester-pays / large)

```bash
bash scripts/download_cpgcorpus_gse.sh
# or: bash scripts/download_cpgcorpus_background.sh gse
```

Not required for milestones 2–7. See [`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md).

---

## 0. Bootstrap / scaffold

- **Status:** `done`
- **Done when:** `/data` paths, DuckDB catalog schema, CLI (`doctor` / `catalog` /
  shallow `inspect`), segment ops, flat + hierarchical model modules, unit tests,
  and experiment YAML exist and pass CI.
- **Evidence:** Stage 0 scaffold on `main`; `make lint typecheck test` green.

---

## 1. Download and inspect one small source

- **Status:** `done`
- **Done when:** One tiny CpGCorpus subset or single GSE/GPL is on disk under
  `$MBS_DATA_ROOT/raw`, and a sanitized report under `reports/inspection/`
  validates file layout, sample alignment, beta ranges, missingness, and platform
  metadata. Catalog + inspector proven on real data (not only fixtures).
- **Evidence:** `GSE125367` / `GPL21145` inspected via
  `mbs inspect cpgcorpus-gpl`; report at
  `reports/inspection/GSE125367_GPL21145/` (44 samples, 865919 probes, perfect
  GSM_ID alignment, betas in `[0,1]`, platform `GPL21145`). Unit coverage in
  `tests/unit/test_inspect_cpgcorpus.py`.
- **Note:** Primary ongoing open source is EWAS Data Hub (ADR 0002). This
  milestone’s CpGCorpus evidence stands; do not re-open milestone 1 to switch
  sources.
- **Next action:** Milestone 6 — hierarchical region model.

---

## 2. Build the canonical annotation graph

- **Status:** `done`
- **Done when:** Stable locus registry and first
  `probe → locus → region → gene` mapping exist. Keep simple: promoter, body,
  UTR, and a few annotation flags. Do **not** expand to full MethylGPT /
  MethylCapsNet topology yet.
- **Evidence:** Graph release
  `graph-grch38-gencode38-five-role-v1` under
  `$MBS_DATA_ROOT/canonical/graphs/graph-grch38-gencode38-five-role-v1/`
  (genes/regions/edges/BED + `graph_manifest.json`); locus registry under
  `$MBS_DATA_ROOT/canonical/annotations/` (`loci`/`probes`/`probe_locus_edges`
  + `annotations_manifest.json`). Validation report:
  `reports/inspection/annotation_graph_v1/` (~1.08M loci, 19937 protein-coding
  genes, five roles, island/QC flags). Built via `mbs graph build` from
  InfiniumAnnotation (HM450/EPIC/EPICv2) + GENCODE v38 + UCSC CpG islands.
  Unit tests: `tests/unit/test_annotation_graph.py`. Build plan:
  [`plans/milestone-2-canonical-annotation-graph.md`](plans/milestone-2-canonical-annotation-graph.md).
- **Depends on:** (1) at least partially (platform/probe IDs known).
- **Leverage (read-only vendor references; do not runtime-import):**
  - [`vendor/infinium_annotation`](https://github.com/zhou-lab/InfiniumAnnotation)
    (Zhou-lab InfiniumAnnotation): hg38 probe ordering, genomic coordinates,
    mapping/quality masks, KnowYourCG feature sets for HM450/EPIC/EPICv2/….
  - [`vendor/epicv2_manifest`](https://github.com/bethan-mallabar-rimmer/EPICv2_manifest):
    EPICv2 reannotation code; Zenodo tables under `$MBS_DATA_ROOT/raw/manifests`.
  - [`vendor/methylcapsnet`](https://github.com/Christensen-Lab-Dartmouth/MethylCapsNet):
    regulatory capsule / typed-region grouping approach to adapt for a simple
    Stage 0 region taxonomy (not the full CapsNet model).
  Convert or export needed tables into `$MBS_DATA_ROOT/canonical/annotations`
  (or graphs); keep bulky vendor blobs out of the Python runtime path.
- **Next action:** Milestone 6 — hierarchical region model.

---

## 3. Export static locus features

- **Status:** `done`
- **Done when:** Offline CpGPT sequence-adapter embedding artifact is the default
  static feature, with a complete static-feature manifest (commit, checkpoint
  hash, vocabulary/locus-table hash, dims, dtype, genome build, export command).
  MethylGPT token priors remain ablation-only.
- **Evidence:** Feature set `cpgpt2m_adapter_128_v1` under
  `$MBS_DATA_ROOT/canonical/static_features/cpgpt2m_adapter_128_v1/`
  (`embeddings.zarr` 128-d float16, `loci.parquet`, schema-valid `artifact.json`).
  Inspection report: `reports/inspection/static_features_cpgpt2m_v1/`
  (~1.076M / 1.082M loci mapped, rate ≈0.994). Export via
  `mbs features export-cpgpt` (optional `--extra cpgpt` for downloads only;
  adapter forward uses `mbs.static_features.cpgpt_adapter` to avoid torchtune /
  torchao breakage on the MBS torch pin). Unit tests:
  `tests/unit/test_static_features.py`. Build plan:
  [`plans/milestone-3-static-locus-features.md`](plans/milestone-3-static-locus-features.md).
- **Depends on:** (2) locus registry.
- **Next action:** Milestone 6 — hierarchical region model.

---

## 4. Convert one pilot matrix into canonical storage

- **Status:** `done`
- **Done when:** One pilot source from **EWAS Data Hub** is written in
  project-local canonical matrix format; slices round-trip correctly from raw
  file to matrix store (checksum / equality checks in tests or inspection
  report). Prefer a labeling GSE already under `EWAS_db/` (see
  [`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md) / [`EWAS_DATA.md`](EWAS_DATA.md))
  or a small baseline subset under `download/`. Do **not** default the pilot to
  CpGCorpus Arrow.
- **Evidence:** Hub labeling study `GSE35069` (60 `GSM*.txt`) under
  `$MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db/GSE35069/`; canonical store
  `matrix-gse35069-ewasdb-v1` at
  `$MBS_DATA_ROOT/canonical/matrices/matrix-gse35069-ewasdb-v1/`
  (`betas.zarr` shape `[60, 485470]` float32, sample/locus parquet indices,
  schema-valid `matrix_manifest.json`). Round-trip PASS in
  `reports/inspection/GSE35069_ewas_db/` (max_abs_diff=0). Convert via
  `mbs matrix convert`; targeted download via
  `make download-ewas-study STUDY=GSE35069`. Unit tests:
  `tests/unit/test_matrix_store.py`.
- **Depends on:** (1), (2).
- **Next action:** Milestone 6 — hierarchical region model.

---

## 5. Train the flat DeepRVAT-style baseline

- **Status:** `done`
- **Done when:** Exact CpG-to-gene max-pooling baseline runs end to end: overfits
  a tiny fixture, then trains on the pilot source. Checkpoints + resolved config
  under `$MBS_ARTIFACT_ROOT`.
- **Evidence:** `mbs train flat --overfit-fixture` reaches accuracy 1.0
  (`artifacts/runs/stage0-flat-overfit-fixture/`). Pilot cell-type CE on
  GSE35069 (`matrix-gse35069-ewasdb-v1` + CpGCorpus metadata labels, 10 classes,
  donor-grouped 4/2 split) via `CUDA_VISIBLE_DEVICES=0 mbs train flat --config
  configs/experiment/stage0_flat_pilot.yaml --run-id stage0-flat-gse35069-v1`
  on a single Ada GPU (`cuda:0`). Artifacts:
  `$MBS_ARTIFACT_ROOT/runs/stage0-flat-gse35069-v1/` (resolved config, metrics,
  split, environment) and
  `$MBS_ARTIFACT_ROOT/checkpoints/stage0-flat-gse35069-v1/` (`best.pt`,
  `last.pt`, checksums). Unit tests: `tests/unit/test_training_flat.py`. Plan:
  [`plans/milestone-5-flat-deeprvat-baseline.md`](plans/milestone-5-flat-deeprvat-baseline.md).
- **Depends on:** (4). Model module scaffolding alone is not sufficient.
- **Next action:** Milestone 5c — multitask shared encoder (5b / 5b′ / 5b″ done;
  do not wait on disease/cancer zips).

---

## 5b. Phenotype registry and multi-pack eval

- **Status:** `done`
- **Done when:** Versioned phenotype/source registry exists; wave-1 Hub packs
  (age, tissue, disease) have family download tooling and sample-info Parquet
  export; evaluation metrics + study-grouped splits are coded and documented;
  flat training logs to TensorBoard/JSONL; first benchmark report covers
  GSE35069 smoke plus age and tissue study-holdout runs (disease subset
  allowed if documented). Atlas remains validation-only.
- **Evidence:** Registry `configs/data/phenotype_registry.yaml` + schema;
  `scripts/download_ewas_phenotype_family.sh` / `make download-ewas-family`;
  sample-info Parquet under `$MBS_DATA_ROOT/canonical/phenotypes/` for age /
  tissue / disease; `src/mbs/evaluation` + unit tests; TensorBoard/JSONL on
  flat loop; fixtures
  `artifacts/runs/stage0-5b-{tissue,age}-holdout-fixture/` plus existing
  GSE35069 pilot; report
  `reports/inspection/stage0_5b_benchmark/`. Plan:
  [`plans/milestone-5b-phenotype-registry-eval.md`](plans/milestone-5b-phenotype-registry-eval.md);
  [ADR 0003](adr/0003-milestone-5b-phenotype-registry.md).
- **Depends on:** (5).
- **Next action:** Milestone 5c (5b′ / 5b″ already done).

---

## 5b′. Hub / Atlas metadata structure (pre-5c gate)

- **Status:** `done`
- **Done when:** Atlas small tables (studies, cohorts, trait×trait) and unpacked
  DataHub `sample_*.txt` packs are profiled with parse recipes, join keys,
  family→column map, and cross-pack ID overlap; sample-info Parquet export works
  from unpacked `.txt` when zips are missing; docs link the report.
- **Evidence:** `uv run mbs inspect ewas-metadata` →
  `reports/inspection/ewas_metadata_structure/`;
  [`EWAS_METADATA.md`](EWAS_METADATA.md);
  [`plans/ewas-metadata-structure.md`](plans/ewas-metadata-structure.md);
  `tests/unit/test_ewas_metadata.py`; nine unpacked families under
  `reports/inspection/ewas_datahub_samples/` (incl. ancestry=`sample_race.txt`, bmi, cancer).
- **Depends on:** (5b) sample-info / registry path.
- **Next action:** **Required reading** at 5c start (do not invent Hub columns;
  do not join Atlas `study_ID` to Hub `project_id` by raw equality). No further
  metadata milestone work before coding 5c.

---

## 5b″. Real Hub pack → matrix + study-grouped eval

- **Status:** `done`
- **Done when:** Downloaded Hub profile packs convert to canonical matrices for
  feasible families; registry lists study IDs / platform / label type / split
  role / matrix paths; study-grouped train/val/test runs on **real** Hub
  matrices (not synthetic fixtures) with TensorBoard + JSONL; per-family and
  combined reports exist. Disease/cancer may remain registered-but-blocked if
  packs are incomplete, with explicit notes.
- **Evidence:** `mbs matrix convert-pack`; matrices
  `matrix-hub-{age,tissue,blood,brain}-studyholdout-v1` under
  `$MBS_DATA_ROOT/canonical/matrices/`; registry benchmark entries in
  `configs/data/phenotype_registry.yaml`; runs
  `artifacts/runs/stage0-hub-*-studyholdout-v1/` (TB + `metrics.jsonl`,
  `model_public_name: deepMAT`); reports
  `reports/inspection/stage0_hub_real_benchmark/`. Scripts:
  `scripts/convert_hub_pack_subsets.sh`, `scripts/train_hub_real_benchmarks.sh`.
  **Caveats (do not block 5c):** disease/cancer profile zips still incomplete on
  disk; single-tissue study-holdout yields disjoint CE classes (0% holdout
  accuracy expected — plumbing gate, not biology). Age external-test MAE logged
  in years.
- **Depends on:** (5b), (5b′).
- **Next action:** Start Milestone 5c on the **ready** age + tissue assets below.
  Re-convert disease/cancer and redesign tissue holdouts for shared classes when
  packs finish — as follow-ons, not 5c start gates.

---

## 5c. Multitask shared encoder (Hub packs)

- **Status:** `pending` — **ready to start with assets already on disk**
- **Done when (MVP — sufficient to mark `done`):**
  - Unified `sample_phenotype_table.parquet` joins Hub sample-info →
    `sample_id` / `study_id` / per-task masks using
    [`EWAS_METADATA.md`](EWAS_METADATA.md) column contracts
  - One shared flat `FlatDeepSet` trains with **linear** age (MSE/Huber) +
    tissue (CE) heads; unlabeled heads masked out per sample
  - Study-grouped train/val/test (no study in more than one role) on **real**
    Hub-derived matrices (`matrix-hub-age-*` + `matrix-hub-tissue-*` or a
    merged multitask matrix build from those packs)
  - Checkpoints + resolved multitask config + TensorBoard/JSONL under
    `$MBS_ARTIFACT_ROOT`; `model_public_name: deepMAT`
  - Unit tests on synthetic multitask fixtures (no full Hub matrices in CI)
- **Explicitly not required for MVP `done`:** complete disease/cancer profile
  zips; disease/cancer aux heads; blood/brain as tissue CE classes; Atlas
  joins; biological tissue accuracy on disjoint single-tissue holdouts.
- **Optional follow-ons (same milestone or soon after, when downloads finish):**
  masked disease/cancer aux heads; blood/brain as **domain aux** after a tissue
  ontology (do not dump into primary tissue CE); shared-class tissue holdouts.
- **Ready inventory (do not wait):**

  | Asset | Status |
  |-------|--------|
  | `age_methylation_v1.zip` + age sample-info parquet | Ready |
  | `tissue_methylation_v1.zip` + tissue sample-info parquet | Ready |
  | `matrix-hub-age-studyholdout-v1` / `matrix-hub-tissue-studyholdout-v1` | Ready |
  | `blood` / `brain` packs + matrices + sample-info | Ready (domain aux only) |
  | Sample-info parquet for disease/cancer/ancestry/bmi | Ready (labels) |
  | `disease_methylation_v1.zip` / `cancer_methylation_v1.zip` | Incomplete (wget in progress) — **skip for MVP** |
  | [`EWAS_METADATA.md`](EWAS_METADATA.md) + ewas_metadata_structure report | Ready |

- **Depends on:** (5b), (5b′), (5b″) — all `done`. Read
  [`EWAS_METADATA.md`](EWAS_METADATA.md) before inventing phenotype columns.
- **Plan:** [`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md);
  draft config [`../configs/experiment/stage0_flat_multitask.yaml`](../configs/experiment/stage0_flat_multitask.yaml);
  schema [`../schemas/sample_phenotype_table.schema.json`](../schemas/sample_phenotype_table.schema.json).
- **Next action:** (1) Build `canonical/phenotypes/sample_phenotype_table.parquet`
  for age+tissue (masks; `450K`→`HM450`); (2) tissue ontology / min-n class
  filter; (3) `multitask.py` masked loss + heads on `FlatDeepSet`; (4) train
  study-grouped age+tissue holdout on real Hub matrices; then Milestone 6.
  Wire disease/cancer aux only after those profile zips verify complete.

---

## 6. Add the hierarchical model

- **Status:** `pending`
- **Done when:** Region layer is trained after the flat baseline is stable;
  promoter/body (and related roles) can be compared to the flat model on the
  same multitask / pilot folds.
- **Depends on:** (5c) preferred; (5b) at minimum if 5c deferred by ADR.

---

## 7. Run study-grouped cross-fitting

- **Status:** `pending`
- **Done when:** Out-of-fold MBS scores, age predictions, and tissue predictions
  are generated with leakage controls (no sample/donor/replicate/held-out study
  scored by a model that saw it). Score matrix + fold-assignment hash stored.
- **Depends on:** (5c) at minimum; (6) preferred for hierarchical OOF scores.

---

## 8. Optional layers (after core pipeline is stable)

- **Status:** `deferred`
- **Rule:** Do not start these until milestones 1–7 produce a real model
  pipeline. Full vision: [`STRATEGIC_PLAN.md`](STRATEGIC_PLAN.md).

### Stage 1+ roadmap (deferred candidates)

| Candidate | Intent / acceptance hint |
|-----------|--------------------------|
| PROTRIDER-style conditional AE | Student-t NLL + missingness mask; two-sided CpG tail probabilities as optional features |
| ComBat-met (rpy2) | Beta-regression batch correction for user/custom IDAT or uncorrected cohorts; assert corrected betas stay in `[0, 1]`. Not required for Data Hub GMQN baselines |
| EPICv2 IlmnID collapse | Strip technical suffixes, groupby core CpG ID, mean beta before matrix build |
| REGENIE association | Export MBS×2 as pseudodosages (BGEN/VCF); Step 1 + Step 2; Firth for imbalanced binary traits; QQ / λ near 1 |
| EWAS Atlas enrichment | Compare significant gene–trait hits to Atlas curated associations / pathway enrichment |
| Richer annotations | MethylGPT priors, MethylCapsNet-inspired topology, intergenic tiles |
| Epivariants / episignatures | Explicit epivariant calling and clinical episignature work |

---

## Agent checklist

Before claiming a milestone `done`:

1. Acceptance criteria above are met with evidence (tests, report paths, artifact
   manifests).
2. Status in this file is updated in the same change set.
3. Required checks from `AGENTS.md` pass for code changes.
4. Do not mark complete because a stub module or download script exists.
