# Stage 0 scientific pipeline TODO

Authoritative milestone list for coding agents. Update status here when a
milestone is **truly done** (acceptance criteria met), not when scaffolding
exists. Cursor agents must read this file at session start and after finishing
work; see `.cursor/rules/pipeline-todo.mdc`.

Status values: `done` | `in_progress` | `pending` | `blocked` | `deferred`

**Numbering:** after **7F**, use integers **8, 9, 10, …** (optional `a/b/c`
sub-tracks). Historical names `7G` / `7G′` / `7H` / “Milestone 7 OOF” are
**aliases only** — see [`plans/MILESTONE_INDEX.md`](plans/MILESTONE_INDEX.md).
On-disk `stage0_7g_*` / `run_7g_*` / `stage0_7h_*` IDs stay unchanged.

## Roadmap (one screen)

```text
1–6          done     bootstrap → flat/hierarchical v0.1 freezes
7A–7F, 7E′   done     release → nine-pack → architecture → CV → cascade (no TBS)
8            done     methylation-only full eval + tissue probe          (was 7G)
9            done     gene-only architecture on ATS (9a–9d screens)      (was 7G′)
10           ← NOW    pretrained MBS/RBS scale / nine-pack campaign      (was 7H)
  10a        done     5-combo grid; P2-G scalar max/max LOCKED cascade finalist
  10a-cap    done     N-light rho_hidden=64 adopted (mbs_e2e 0.308/14.7/0.852)
  10a-agecov rejected age head tissue+sex conditioning — all 3 metrics worse
  10a-warm   done     both vector warms done; max→max 0.342/13.41/0.851 (≤ P2-G tissue)
  10a++      done     one-hop smokes (G0≫G1; seeds 42/43/44 distinct)
  10b        done     ATS seed-43 pooling (P2-G still leads)
  10c        partial  label_status + BMI/ancestry; freeze-reuse n≥600
  10d        pending  reference checkpoints after staged recipe + enet readout
  10e        in_prog  cascade S1–S4 still gated; N-light nested enet landed
11           deferred fold-selected panel Stage B — parallel, not a hard gate on 12
12           in_prog  N-light 5×6 first (GPU2 after S1); cascade 5×6 still gated
13           deferred expression aux after OOF
14           deferred optional a–f
```

Live board: [`plans/milestone-10-pretrained-mbs-rbs.md`](plans/milestone-10-pretrained-mbs-rbs.md)
→ [`plans/milestone-7h-pretrained-mbs-rbs-campaign.md`](plans/milestone-7h-pretrained-mbs-rbs-campaign.md)
→ [`../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md)
→ [`../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`](../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md).

### Next steps (ordered)

**GPU-2 policy:** max-load N-light 5×6 after dense S1 (skip naive transplants).
Milestone **11** Stage B GPU is optional / parallel.

**Priority:**

1. ~~**5-combo pooling grid**~~ — done. **P2-G** locked on joint e2e (0.355 / 13.431 / 0.853).
2. ~~**N-light rho=64**~~ — done as encoder; **nested enet 3-fold 0.368 / 9.88 / 0.803**.
3. ~~**Age covariates**~~ — rejected.
4. ~~**Vector warm-starts**~~ — **both done**. Do not replace P2-G on e2e; re-rank under enet.
5. ~~**One-hop correctness smokes**~~ — done. G0 ≫ G1; seeds 42/43/44 distinct.
6. **NOW — Milestone 12 N-light 5×6** on GPU 2 (max VRAM) after dense S1.
   Split `hub-nine-pack-5fold-v1`. Product readout `mbs_enet_nested`.
   Runner: `scripts/run_12_nlight_oof.sh`. Cascade 5×6 still waits on 10e.
7. **CPU:** finish P2-G nested enet 3-fold (re-rank vs N-light under the same metric).
8. **Freeze-reuse** traits with **n≥600** (BMI, ancestry, cancer pack, Alzheimer’s 945).
   PD/stroke stay below the case bar. Catalog 173k is **not** free extra data
   (only ~110k HM450); encoder must stay gene-invariant for EPIC/ONT later.
9. **Then cascade 5×6** after 10e S1–S4 smoke. Milestone **11** parallel.
10. **Deferred post-OOF:** CpGPT/positional, sex-chr imputation, clock age imputation.

**Still blocked auto:** cascade 5×6, pack-mask trait heads, blood/brain heads,
mixing EPIC/ONT into this HM450 OOF.

### Trustworthy ATS numbers (`explicit_only`, 51 375 gene-linked CpGs, test)

| Question | Best arm | Tissue macro-F1 | Caveat |
|----------|----------|----------------:|--------|
| Classical on gene panel | `C-mvalue-enet-G` | **0.388** (±0.018) | Still the tissue leader |
| Best MBS readout | `P2-G` **`mbs_enet`** | 0.385 (±0.053) | Same encoder as e2e; elastic-net heads |
| Best cascade product path | `P2-G` **`mbs_e2e`** | 0.373 (±0.038) | max/max, 15 epochs; P4 mean 0.370 |
| One-hop mean (16-ep) | `N-light-gene-mean` | **~0.378** e2e | Near P2; age MAE ~17 |
| Ablation | `m_only` preferred vs `full` | — | Matched 16-ep confirmed |
| Seed-mask | `G0` beats G1–G3 | — | **Not adopted** |
| Invalid (do not cite) | pre-fix `mbs_e2e` on P*-G | ~0.67–0.70 | train+val+test leak |

**Refs:** table above is the ATS-scale (13,548-sample) reference; P2-G there
predates the scale decision below. **At nine-pack scale (34,234 samples),
P2-G scalar max/max is the LOCKED cascade finalist** (5-combo grid + age-cov
reject). Cold vector loses; **max→max warm-start** improves to 0.342 / 13.406 /
0.851 but still trails P2-G tissue (−0.013) — **not a finalist change**.
**N-light@64** nested enet 3-fold is **0.368 / 9.88 / 0.803** vs e2e
**0.308 / 14.73 / 0.852**. Milestone **12 N-light 5×6 is starting**; cascade
5×6 still waits on 10e. Architecture stays gene-invariant for later EPIC/ONT.
Live numbers: [`analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md).
**Platform (this OOF):** HM450 nine-pack — no mixed-platform claim yet.

Frozen freezes (do not overwrite): **deepMAT-flat-v0.1** /
**deepMAT-hierarchical-v0.1** / **deepmat-data-age-tissue-sex-v1**.

Hub census (refresh 2026-09-07): **173 076** samples, **1 763** studies,
**21** matrix artifacts including nine Hub full packs — see
`reports/inspection/deepmat_data_v1/`. EWAS_db ingest incomplete
(`mirror_complete=false`) and **not** a Milestone 10/13 gate.

Programme: [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).
ADRs: [0007](adr/0007-crossfit-prerequisites.md) (OOF prerequisites; OOF =
Milestone **12**), [0008](adr/0008-score-identifiability.md),
[0009](adr/0009-drop-tbs-scores.md), [0010](adr/0010-gene-allocation-policy.md).


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
- **Intentional non-goals until later milestones:** committing the shallow
  whole-corpus inventory under `reports/inspection/cpgcorpus/` (GSE/GPL report
  for milestone 1 is enough). MethylGPT token-prior export remains ablation-only
  (not Stage 0 default). **Parquet → DuckDB catalog population is Milestone 7A**
  (no longer deferred): see [ADR 0005](adr/0005-catalog-matrix-independence.md).
- **Hub disease profile zip:** completed 2026-08-11 via
  `scripts/download_disease_pack_resilient.sh` (exact remote size + EOCD).
  Earlier failures were CNCB connection drops / bogus HTTP 416, not disk
  space. Cancer pack was already OK. **EWAS_db** All-Data download still
  running (`download_ewas_datahub.sh EWAS_db`; ~1 582/1 989 studies visited).
  Post-download hook + failure audit:
  [`plans/data-infrastructure-improvements.md`](plans/data-infrastructure-improvements.md).
  **GEO sample backfill:** repaired pilot re-validated; **batch-50 merged**
  (45 762 GSM → catalog `geo_metadata_backfill`; report
  `geo_backfill_batch/`). Eligibility-by-study audit done
  (`eligibility_by_study.md`). Immutable phenotype release
  `deepmat-data-geo-dev-v1` **built** (Hub/GEO arms; no matrix convert yet) —
  [`plans/geo-enriched-training-release.md`](plans/geo-enriched-training-release.md).
  Do **not** mutate frozen ATS; do **not** wire GEO into age/tissue training or
  enlarge the crawl until Milestone **10/11** gates allow. **Wave-1 label fix**
  (2026-09-07): tissue mapped **14 663→43 886** via `remap_geo_tissue` /
  source_name fallback; catalog GEO tissue **36 861** (0 `nan`); Atlas still
  **182**/1 763 matched (PMID ceiling) —
  [`plans/improve-labels-study-context.md`](plans/improve-labels-study-context.md).
  Plans:
  [`plans/geo-metadata-backfill-pre-scale.md`](plans/geo-metadata-backfill-pre-scale.md).
  Inventory: [`EWAS_DATA.md`](EWAS_DATA.md), [`DATA_CATALOG.md`](DATA_CATALOG.md),
  `reports/inspection/raw_inventory/`.

### Useful commands (already safe to re-run)

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
uv run mbs doctor --create-directories
uv run mbs catalog init
uv run mbs catalog refresh-release
uv run mbs catalog validate-release
uv run mbs catalog phenotype-census
uv run mbs catalog trait-eligibility
# or: make catalog-refresh-release
make summarize-ewas-db-failures
make retry-ewas-db-failures   # filters (.+?) artifacts; background-safe
# Multitask flat training checks trait_eligibility for disease/cancer heads
# (override: training.check_trait_eligibility: false)
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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

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

- **Status:** `done`
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
- **Evidence:** `mbs phenotypes build-multitask-table` →
  `matrix-hub-age-tissue-multitask-v1` (299 samples, 16 GSM deduped) +
  `canonical/phenotypes/sample_phenotype_table.parquet` +
  `tissue_ontology.yaml` (5 classes). Train:
  `CUDA_VISIBLE_DEVICES=0 mbs train flat --config
  configs/experiment/stage0_flat_multitask.yaml --run-id
  stage0-flat-multitask-age-tissue-v1` →
  `$MBS_ARTIFACT_ROOT/runs/stage0-flat-multitask-age-tissue-v1/` (TB,
  `metrics.jsonl`, resolved config, split) and
  `checkpoints/stage0-flat-multitask-age-tissue-v1/` (`best.pt` epoch 15,
  `model_public_name: deepMAT`, task=`multitask`). Units:
  `tests/unit/test_multitask.py`, `tests/unit/test_phenotype_table.py`.
  Plan: [`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md).
- **Explicitly not required for MVP `done`:** complete disease/cancer profile
  zips; disease/cancer aux heads; blood/brain as tissue CE classes; Atlas
  joins; biological tissue accuracy on disjoint single-tissue holdouts
  (external-test tissue accuracy 0% expected under current single-tissue
  studyholdout design — same caveat as 5b″).
- **Optional follow-ons:** masked disease/cancer aux heads; blood/brain as
  **domain aux** after ontology; shared-class tissue holdouts.
- **Depends on:** (5b), (5b′), (5b″) — all `done`.
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

---

## 5d. Max-N flat DeepRVAT baseline (age / tissue / sex)

- **Status:** `done`
- **Done when:** Full Hub age/tissue/sex packs convert without study/sample
  caps; GSM-union cohort; **shared** flat deepMAT + **decoupled phenotype
  modules** (age/tissue/sex) with **masked per-trait loss** (DeepRVAT pattern —
  not dynamic head-switching); study-grouped auto split + checkpoints +
  inspection under `reports/inspection/stage0_5d_max_n/`.
- **Evidence:** Uncapped
  `matrix-hub-{age,tissue,sex}-full-v1` (8374 / 5323 / 2978 samples) via
  `scripts/convert_hub_full_packs.sh`; GSM-union
  `matrix-hub-age-tissue-sex-full-v1` (13548 samples) +
  `sample_phenotype_table_age_tissue_sex_full_v1.parquet` (masks age=10002 /
  tissue=7866 / sex=12445); config
  `configs/experiment/stage0_flat_deeprvat_full.yaml`; run
  `stage0-flat-deeprvat-age-tissue-sex-full-v1` (train/val/test 9489/2074/1985;
  47 tissue classes; external tissue acc ~0.666, age MAE ~22 y, sex acc ~0.931;
  `model_public_name: deepMAT`); checkpoints under
  `$MBS_ARTIFACT_ROOT/checkpoints/stage0-flat-deeprvat-age-tissue-sex-full-v1/`;
  report `reports/inspection/stage0_5d_max_n/` via
  `scripts/write_stage0_5d_report.py`.
- **Depends on:** (5c).
- **Plan:** [`plans/milestone-5d-max-n-flat-baseline.md`](plans/milestone-5d-max-n-flat-baseline.md).
- **Next action:** Milestone **7E** (gene-only arms may start; full 3×2 waits
  on graph-v2).
- **Freeze name:** `deepMAT-flat-v0.1` (do not overwrite this run).

---

## 6. Add the hierarchical model

- **Status:** `done`
- **Done when:** Region layer is trained after the flat baseline is stable;
  promoter/body (and related roles) can be compared to the flat model on the
  same multitask / pilot folds; unmapped loci are retained on a residual path
  (not `__unassigned__` gene pooling) with mapped vs residual eval slices.
- **Depends on:** (5d) preferred; (5c) at minimum if 5d deferred by ADR.
- **Plan:** [`plans/milestone-6-hierarchical-region-model.md`](plans/milestone-6-hierarchical-region-model.md).
- **Evidence:** Residual retention policy in batch/matrix/models/hier train
  (ADR 0004); annotation-status masks; mapped vs residual eval slices.
  Uncapped run `stage0-hier-deeprvat-age-tissue-sex-full-v1` on
  `matrix-hub-age-tissue-sex-full-v1` with reused 5d split (9489/2074/1985);
  best epoch 13; topology 19554 genes + residual slot, 108070 residual cols,
  five GENCODE roles only. External vs flat 5d: tissue acc 0.598 vs 0.666,
  age MAE 27.8 vs 22.0 y, sex 0.934 vs 0.931. Ablations show mapped≈full and
  residual_only near chance for tissue/sex. Report:
  `reports/inspection/stage0_6_hierarchical/` via
  `scripts/write_stage0_6_report.py`. One-scalar residual is a **bottleneck
  ablation**, not evidence that noncoding CpGs are uninformative
  ([ADR 0006](adr/0006-multipath-noncoding-scores.md)).
- **Freeze name:** `deepMAT-hierarchical-v0.1` (do not overwrite; not the
  preferred phenotype model vs flat v0.1).
- **Next action:** Milestone **7E** (graph-v2 on disk; multi-path unblocked).

---

## 7A. Harmonized data release and phenotype census

- **Status:** `done`
- **Done when:**
  - Versioned release `$MBS_DATA_ROOT/canonical/releases/deepmat-data-v1/` with
    `release_manifest.json` (source checksums, retrieval dates, preprocessing,
    probe universe, genome build, graph/static versions, phenotype families,
    dedup decisions, code commit)
  - DuckDB catalog **populated** from Parquet/manifests (not schema-only):
    `source_release`, `study`, `platform`, `sample`, `sample_source_membership`,
    `assay_file`, `phenotype`, `sample_phenotype` (long-form + `source_family`),
    matrix inventory tables, `fold_assignment`, `artifact`, `experiment`
  - Census views and reports answer unique GSM vs pack-row sum, studies,
    platforms, phenotype prevalence, pack overlap, label conflicts, confounding
  - `trait_eligibility` table with initial cutoffs from
    [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md)
  - CLI: `mbs catalog refresh-release`, `validate-release`, `phenotype-census`,
    `trait-eligibility`
- **Refresh follow-on (does not reopen 7A):** remaining census fields — within-study
  age/BMI ranges, documented controls, donor/replicate IDs, metadata-only
  predictability, full platform×tissue support — listed in the programme brief.
- **Evidence:** Release
  `$MBS_DATA_ROOT/canonical/releases/deepmat-data-v1/` (schema-valid
  `release_manifest.json`, populated `catalog/catalog.duckdb` +
  `catalog/tables/*.parquet`, phenotype long-form, matrix pointers, ingested
  5d `split.json`). Live refresh 2026-08-25T11:15Z: **121 931** unique GSM
  (Hub pack members **34 234**; rest EWAS_db-only), **1 325** studies,
  **216 476** phenotype rows; EWAS_db listing **924** local studies /
  **92 971** GSM files (`mirror_complete: false`, advertised 1989). Census +
  eligibility match this refresh:
  `reports/inspection/deepmat_data_v1/`. Unit tests:
  `tests/unit/test_catalog_release.py`. Plan:
  [`plans/milestone-7a-harmonized-release.md`](plans/milestone-7a-harmonized-release.md).
- **Depends on:** (6); EWAS_db mirror **not** required.
- **Plan:** [`plans/milestone-7a-harmonized-release.md`](plans/milestone-7a-harmonized-release.md);
  [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md);
  [ADR 0005](adr/0005-catalog-matrix-independence.md),
  [ADR 0007](adr/0007-crossfit-prerequisites.md).
- **Next action:** Milestone **7E**. EWAS_db refresh is optional and not a
  7E gate. Do not cite the hyphen inspection dir.

---

## 7B. Complete canonical Hub matrices

- **Status:** `done`
- **Done when:** Disease, cancer, blood, brain, BMI, and ancestry packs convert
  to canonical full matrices; BMI/ancestry supported in pack converter maps;
  **stream probe chunks directly** to compressed Zarr (no full dense RAM stack);
  **per-sample** platform provenance (not one HM450 map for merged unions);
  probe-collapse records **all** contributing probe IDs (mean/robust mean, not
  lexicographic-first); disease/cancer multi-label via long-form (no
  `dict[gsm]=row` overwrite); **content** checksums (not filename/size);
  overlapping GSM betas **verified** (do not silently take the first pack);
  deduplicated union or virtual multi-store index documented;
  evidence under `reports/inspection/stage0_7b_hub_matrices/`.
- **Evidence:** six full matrices under `$MBS_DATA_ROOT/canonical/matrices/`;
  `hub_pack_matrix_index.parquet`; overlap `concordant` (0 discordant);
  `reports/inspection/stage0_7b_hub_matrices/summary.{md,json}`; unit tests for
  converter path. Disease `12218×482387` (14501 phenotype rows); cancer
  `10101×482387`. Watcher finalize ~2026-08-24T17:22Z.
- **Depends on:** (7A).
- **Plan:** [`plans/milestone-7b-complete-hub-matrices.md`](plans/milestone-7b-complete-hub-matrices.md);
  [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).
- **Next action:** Milestone **7E** (graph-v2 on disk; full 3×2 unblocked).

---

## 7C. Supervised architecture corrections

- **Status:** `done` (fixture acceptance)
- **Done when:**
  - Trainer P0: deterministic epoch shuffle; **used** `batch_token_budget`;
    task/study-balanced sampling; real donor/replicate identifiers
  - All phenotype heads center/mask present scores consistently
  - Constraint-aware study-grouped splits (tissue-class, task-mask, age quantile,
    platform, cases/controls, donor/replicate; not sample-count only)
  - Emit macro-F1, balanced accuracy, RMSE, R², correlations, AUROC/AUPRC,
    calibration; study/platform/tissue-stratified reports
  - Controls: static-only, coverage-only, **metadata-only**, label permutation
  - Mapped loci missing CpGPT kept with `static_present=False` (not dropped);
    residual zeros carry a missingness flag
  - **Score orientation anchor** ([ADR 0008](adr/0008-score-identifiability.md))
    before any OOF average; predictive MBS ≠ constraint/LOEUF score
  - Graph v2: RBS + TBS; first direct branch sparse \(D_k=\sum w_{k,c}z_{s,c}\)
    (elastic-net / group sparsity); independently trained branch ablations
  - Parameter-matched width/activation/dropout/norm when comparing flat vs hier
- **Evidence:** `tests/unit/test_stage0_7c.py`;
  `reports/inspection/stage0_7c_architecture/`.
- **Residual follow-ons** (do not reopen Done when; track in
  [`plans/milestone-7c-supervised-architecture.md`](plans/milestone-7c-supervised-architecture.md)):

  **Hub join (after 7B):**
  - Long-form multi-hot verified on `matrix-hub-disease-full-v1` /
    `matrix-hub-cancer-full-v1` (`hub_longform_ready` + sidecar load).
  - Short train smoke **done** (`stage0_flat_hub_disease_multilabel_smoke.yaml`;
    `reports/inspection/stage0_7c_hub_disease_smoke/`).
  - Holdout AUROC/AUPRC/ECE emission **done** (binary sex/tissue; multilabel
    when both classes observed under mask).

  **Topology residual (closed; not Done when reopen):**
  - Full-genome `graph-grch38-gencode38-cgi-tile-v2` on disk (RBS n≈18 356 ≫ 72;
    TBS n≈5 446; inspection `annotation_graph_cgi_tile_v2/`)
  - Multi-system hier index (`region_systems`; default still gene-only)
  - Train-time region-system masks for `rbs`/`tbs` arms
  - Plan: [`plans/milestone-7c-graph-v2-topology.md`](plans/milestone-7c-graph-v2-topology.md)

  **Landed in residual polish (fixtures + Hub join):** train-path
  `apply_orientation` + honest `score_manifest.json`;
  `load_longform_multilabel` + masked BCE + `lambda_disease`/`lambda_cancer`.
- **Depends on:** (7B) done for Hub data; trainer closed on age/tissue/sex
  fixtures (still **before 7E**).
- **ADRs:** [0006](adr/0006-multipath-noncoding-scores.md),
  [0008](adr/0008-score-identifiability.md).
- **Plan:** [`plans/milestone-7c-supervised-architecture.md`](plans/milestone-7c-supervised-architecture.md);
  [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).

---

## 7D. Fold-fitted normalization ablation

- **Status:** `done`
- **Done when:** Level-1 study-balanced median + **1.4826×MAD** on train-fold
  M-values; persist \(\mu,\sigma\) and hashes; novel loci `z=0` +
  `norm_present=False` (not discarded); A (beta+M) vs B (A + robust z) on
  identical folds (fixtures **and** DeepRVAT Hub ATS smoke). Do **not**
  overwrite Hub GMQN canonical betas. Levels 2 (bounded residual MLP +
  LayerNorm/RMSNorm) and 3 (fold-isolated masked AE) documented as later
  ablations; select on phenotype/stability, not reconstruction loss.
- **Evidence:** `tests/unit/test_stage0_7d.py`; Hub smoke runs
  `stage0-7d-level1-{a,b}`; `reports/inspection/stage0_7d_level1/`.
- **Depends on:** (7C) at least for shared train path; Level-1 can land with 7C.
- **Plan:** [`plans/milestone-7d-fold-fitted-normalization.md`](plans/milestone-7d-fold-fitted-normalization.md);
  [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).


---

## 7E. Development cross-validation (architecture selection)

- **Status:** `done` (winner `N-multipath-l1a`; readable report
  `reports/inspection/stage0_7e_dev_cv/analysis.md`)
- **Prep (done):**
  - `$MBS_DATA_ROOT/canonical/graphs/graph-grch38-gencode38-cgi-tile-v2/`
  - `reports/inspection/annotation_graph_cgi_tile_v2/`
  - Train-time RBS/TBS feature masks (eval-time masking is not an ablation)
- **May start now:** parameter-matched flat/hier on frozen **ATS** (age/tissue/sex
  GSM-union, 13 548; freeze `deepmat-data-age-tissue-sex-v1`) + Level-1 A/B +
  mean/elastic-net + CpGPT ablation; independently trained RBS/TBS/direct arms
  on graph-v2. Disease/cancer heads with **masked** unknown≠control belong on
  Hub packs in **7E′** (do not skip them; do not treat unlabeled rows as
  controls). More EWAS_db files do **not** enlarge ATS.
- **Done when:** 3 outer study-grouped folds × 2 restarts compare independently
  trained arms: **transparent gene/region mean and elastic-net**;
  **parameter-matched** flat gene-only; **parameter-matched** hierarchical
  gene-only; gene + direct CpG; gene + RBS + TBS + direct; each neural arm
  with/without Level-1 robust-z; **CpGPT inclusion as a separate ablation**.
  Report selects architecture for later milestones (now **9→13**). Eval-time branch masking and
  ordered-prefix holdout eval are not sufficient. Encoder width/GELU/dropout/LN
  must match for flat vs hier.
- **Depends on:** (7C), (7D Level-1); graph-v2 on disk for multi-path arms.
- **Plan:** [`plans/milestone-7e-development-cv.md`](plans/milestone-7e-development-cv.md);
  [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).
- **Next action:** 7E closed. Do **not** ship `N-multipath-l1a` as-is (TBS
  dropped; late fusion was region means; 2-epoch / 8 192-locus budget).
  Proceed to **7F**, then **7G**. Do not start Milestone **7** 5×6 OOF.

---

## 7E′. Hub multitask (age/tissue/sex/disease/cancer) + analysis hygiene

- **Status:** `done`
- **Why a separate step:** Frozen ATS is only 13 548 GSM from three packs.
  Unique Hub GSM already on disk is **34 234**. EWAS_db download progress does
  not add rows to Hub packs. Eligibility `core=False` for disease/cancer is
  about **unknown ≠ control**, not about skipping those heads.
- **Done when:**
  - Train **age, tissue, sex, and disease** (cancer too) with masked BCE /
    masked continuous heads on Hub matrices already converted; unlabeled
    disease/cancer is unknown, never a control
  - New study-grouped split for that Hub-wide cohort; **do not** overwrite
    `deepmat-data-age-tissue-sex-v1` or v0.1
  - Metadata-only control (study/platform/tissue → phenotype) on the **same
    folds as 7E**, reported as a confounding ceiling
  - 7A census follow-ons: donor/replicate IDs when present; within-study
    age/BMI ranges
  - Catalog: alias Hub sample-info `450K` → `HM450` on next
    `mbs catalog refresh-release` (same 450K universe; EPIC is not in these
    nine zips — see plan)
  - Census/eligibility tests use a temp `--report-dir` (hyphen CLI default must
    not be fixture-clobbered)
  - `*.RData` gitignored; no Hub sample blobs committed under
    `reports/inspection/ewas_datahub_samples/`
  - Blood `cell_component` not used as a pack-wide head (~1.1% populated)
  - Parameter-matched flat vs hier YAML (shared with 7E Done when)
- **Not in this list (already done):** graph-v2 + train-time RBS/TBS masks
  ([`plans/milestone-7c-graph-v2-topology.md`](plans/milestone-7c-graph-v2-topology.md)).
- **Depends on:** (7B) matrices; (7E) may run in parallel.
- **Plan:** [`plans/milestone-7e-prime-analysis-hygiene.md`](plans/milestone-7e-prime-analysis-hygiene.md).

---

## 7F. RBS→gene cascade + direct leftover (no tile scores)

- **Status:** `done`
- **Plan (impl):** [`plans/milestone-7f-rbs-gene-direct.md`](plans/milestone-7f-rbs-gene-direct.md);
  [ADR 0009](adr/0009-drop-tbs-scores.md).
- **Report:** [`reports/inspection/stage0_7f_rbs_gene_direct/`](../reports/inspection/stage0_7f_rbs_gene_direct/)
- **Why:** 7E’s winner fused gene + RBS + **TBS** + direct as linear models on
  **presence-aware region means**, not saved neural scores. **Tiles are
  dropped** because 50-CpG bins randomly aggregate leftover loci into a
  feature that is not a typed region. Unassigned CpGs stay **direct**, not
  tiles. **Nearest-gene is allowed** as the gene allocation of an already
  typed RBS (how that region score becomes MBS). ADR 0004 still forbids
  collapsing an unmapped *CpG* into a nearest-gene proxy instead of keeping
  it on the direct path. Graph-v2 tile nodes may stay on disk unused.
- **Locked topology (late fusion):**
  ```text
  CpG → cCRE / enhancer / CGI / DMR / ChromHMM / similar / typed gene region
        → RBS (one score per typed region)
  then three families:
    Direct  — CpG with no region assignment (not tiled)
    RBS     — region score with no gene allocation
    MBS     — gene score = pool of RBS allocated to that gene
              (typed gene region and/or nearest-gene)
  ```
  1. **RBS** — observed CpGs that hit a typed region (cCRE / enhancer / CGI /
     DMR / ChromHMM / similar / gene-body roles) → one score per region.
  2. **MBS** — those RBS allocated to a gene (including nearest-gene) are
     pooled to gene-level scores.
  3. **Direct CpGs** — CpGs with **no** region assignment, in parallel
     (elastic-net / group-sparse or a neural per-locus term on
     fold-normalized z). Not tiled; not forced into a gene as a CpG.
  4. **Late fusion** concatenates **saved** orphan RBS, MBS (gene-aggregated
     RBS), and direct contributions, then a linear (or boosted) phenotype
     head. Fusion of region-mean tables is **not** sufficient.
- **Done when:**
  - Assignment + trainer implement the cascade above on frozen 7E folds
    (`hub-ats-7e-3fold-v1`); fixture tests cover leftover→direct and
    RBS→gene aggregation (nearest-gene allocation of RBS is allowed)
  - Per-sample RBS / gene-RBS / direct score matrices are written and fused
    (not region-mean linear stand-ins)
  - No TBS arm in the 7F model matrix
  - Inspection report under `reports/inspection/stage0_7f_rbs_gene_direct/`
- **Depends on:** (7E) folds + report; (7C) graph-v2; (7D) Level-1.
- **Plan:** [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md)
  (§7F).
- **Next:** Milestone **8** (historical 7G).

---

## 8. Methylation-only full evaluation (alias: 7G)

- **Status:** `done`
- **Plan:** [`plans/milestone-8-methylation-eval.md`](plans/milestone-8-methylation-eval.md)
  → detail [`plans/milestone-7g-methylation-eval.md`](plans/milestone-7g-methylation-eval.md)
- **Why:** Close 7E evaluation gaps (budget, classical methylation-only
  comparators, sex metrics). Do **not** conclude “trees beat Deep Sets” from
  the 7E 2-epoch / 8 192-locus table.
- **Evidence:** `reports/inspection/stage0_7g_methylation_eval/` — ranking
  winner **`C-mvalue-enet`**; cascade weak at first budget. Tissue probe
  P0–P3 historical in `stage0_7g_cascade_tissue_probe/` (P2 late-fusion ~0.38
  — not a clean MBS-only lock).
- **Depends on:** (7F).
- **Next:** Milestone **9**.

---

## 9. Gene-only architecture selection on ATS (alias: 7G′ Stage A)

- **Status:** `done` (no architecture lock; 2×2 pooling retained; seed-mask
  not adopted)
- **Plan:** [`plans/milestone-9-gene-only-architecture.md`](plans/milestone-9-gene-only-architecture.md)
- **Master brief:** [`plans/milestone-7g-prime-matched-probe-lightweight.md`](plans/milestone-7g-prime-matched-probe-lightweight.md)
- **Runners (historical names):** `scripts/run_7g_gene_only_probe.py`,
  `scripts/run_7g_16ep_promotion_resume.sh`, `scripts/run_7g_prime_seed_mask.py`

### 9a — DeepRVAT Tier-1 / gene-only probe grid

- **Status:** `done`
- **Plan:** [`plans/milestone-7g-prime-stage-a-deeprvat-screen.md`](plans/milestone-7g-prime-stage-a-deeprvat-screen.md)
- **Evidence:** `P2-G` best cascade `mbs_e2e` **0.373**; classical enet-G
  **0.388**; vector ≤ scalar; `m_only` best annotation under short budget.
  Report: `reports/inspection/stage0_7g_gene_only_probe/`.

### 9b — Matched 16-epoch promotion

- **Status:** `done` — `next_gate: retain_pooling_2x2` (no pooling lock)
- **Plan:** [`plans/milestone-7g-prime-16ep-promotion.md`](plans/milestone-7g-prime-16ep-promotion.md)
- **Evidence:** full 2×2 cascade pooling within noise; one-hop
  `N-light-gene-mean` near P2 tissue (~0.378); `m_only` still beats `full`
  at matched 16 ep.

### 9c — Age-primary seed-mask

- **Status:** `done` — **seed-masking not adopted** (`G0` beats G1–G3)
- **Plan:** [`plans/milestone-7g-prime-age-seed-mask.md`](plans/milestone-7g-prime-age-seed-mask.md)
- **Bugs fixed en route:** `learning_rate` threading, seed-offset/`-s2`,
  classical tissue array KeyError; gradient clipping; epoch budget raised.
- **Report:** `reports/inspection/stage0_7g_prime_seed_mask/`.

### 9d — Typed-RBS R0–R5 (CPU)

- **Status:** `done` — neural typed aggregator **not** promoted (shuffle Δ≪1 y)
- **Plan:** [`plans/milestone-7g-prime-pre-stage-b.md`](plans/milestone-7g-prime-pre-stage-b.md)

**Note:** Fold-selected panel GPU is **Milestone 11**, not part of 9.

---

## 10. Pretrained MBS/RBS scale campaign (alias: 7H)

- **Status:** `in_progress` ← **current gate**
- **Plan:** [`plans/milestone-10-pretrained-mbs-rbs.md`](plans/milestone-10-pretrained-mbs-rbs.md)
- **Running log / A→B→C→D board:**
  [`plans/milestone-7h-pretrained-mbs-rbs-campaign.md`](plans/milestone-7h-pretrained-mbs-rbs-campaign.md)
- **Goal:** best gene-level (MBS) and region-level (RBS) architecture(s) for a
  pretrained scoring framework; scale from ATS (13 548) to nine-pack
  (`matrix-hub-nine-pack-virtual-v1`, 34 234) when loader/refs are honest.
- **Depends on:** (9) reference arms; (7B) matrices; virtual store wiring.

### 10a — Nine-pack virtual loader + reference arms

- **Status:** `done` — **P2-G scalar max/max locked as cascade finalist**
- **Done:** `RoutedBetas` dense `[:, :n]`; cascade/classical/flat openers;
  alignment tests; split `hub-nine-pack-3fold-v1`; full 5-combo pooling
  grid (3-fold, 15 ep): P2-G scalar max/max wins outright (tissue F1
  **0.355**, age MAE **13.431**, sex AUROC **0.853**) over scalar mean/max,
  scalar max/mean, vector mean→max, vector max→max — all worse on F1 and
  MAE. N-light capacity diagnostic (`rho_hidden` 10→64), read correctly via
  the primary `mbs_e2e` evaluation (an earlier read of this used
  `external_test`, which under-reports tissue F1 for flat runs and gave a
  false "worse" verdict), **adopted as the N-light default**: tissue F1
  0.273→0.308, age MAE 15.057→14.727, sex AUROC 0.742→0.852, all improved.
  Age head tissue+sex conditioning built, tested, **rejected** (full
  3-fold ablation regressed all 3 metrics). Full detail: `vector_vs_scalar.md`.
- **Next:** one-hop smokes; CPU enet on existing checkpoints; **10e** staged
  recipe before any 5×6. Freeze-reuse 10c on frozen scores (not pack masks).

### 10a-warm — Vector warm-start (LP-FT)

- **Status:** `done` — **max→max and mean→max both complete** (2026-09-08)
- **max→max warm:** 0.342 / 13.406 / 0.851 (vs cold 0.333 / 15.408 / 0.834; vs P2-G 0.355 / 13.431 / 0.853)
- **mean→max warm:** 0.345 / 14.578 / 0.862 (vs cold 0.335 / 16.448 / 0.780)
- **Verdict:** warm-start beats cold vector; **does not replace P2-G** as cascade finalist
- **Write-up:** `reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`

### 10a++ — One-hop correctness smokes

- **Status:** `in_progress` (GPU 0) — G0 done (~0.335 / 17.3 / 0.75); G1 + multi-seed running
- **Runner:** `scripts/run_7h_onehop_correctness_smokes.py`
- **Done when:** report under `reports/inspection/stage0_7h_onehop_correctness/`

### 10c — Freeze-and-reuse for disease/cancer (+ BMI/ancestry heads)

- **Status:** `partial` — disease/cancer `label_status` ready; BMI/ancestry
  **joined + heads wired** (freeze-reuse); brain labels honest / head
  deferred; GPU arms not launched.
- **Standing philosophy:** freeze-and-reuse (DeepRVAT) — train a
  gene-invariant encoder once, freeze it, reuse frozen scores for new
  traits via lightweight probes / new linear heads rather than joint
  retraining.
- **Ready:** `label_status` sidecar
  (`sample_phenotype_table_hub_nine_pack_v1.label_status.parquet` +
  [`label_status_census.md`](../reports/inspection/stage0_7h_nine_pack_smoke/label_status_census.md))
  — disease 12,194 / cancer 9,077 usable case+control samples.
  BMI/ancestry: real Hub labels in nine-pack table + `bmi_head` /
  `ancestry_head`; configs stubbed, **not** queued on GPU 0.
- **Plan (disease/cancer GPU):** score once with finalized frozen encoder;
  fit per-fold logistic probes; prefer column selection over seed-mask
  weight gating for sparse disease panels.
- **BMI/ancestry / overview plan:**
  [`plans/label-prep-bmi-ancestry-sample-overview.md`](plans/label-prep-bmi-ancestry-sample-overview.md)
- **Still deferred:** `brain_head` / `blood_head` (blood sparse; brain =
  region catalogue). GPU launch gated on ≥1k-per-arm + M10 review.
- **CPU report:**
  [`reports/inspection/stage0_7h_nine_pack_smoke/trait_hygiene.md`](../reports/inspection/stage0_7h_nine_pack_smoke/trait_hygiene.md)

### 10e — Staged RBS→MBS training + frozen enet (blocks 12)

- **Status:** `pending` — **do not launch Milestone 12 5×6 until a 1-fold
  smoke of this recipe beats P2-G joint e2e and reports RBS/MBS enet**
- **Plan:** [`plans/milestone-10e-staged-rbs-mbs-training.md`](plans/milestone-10e-staged-rbs-mbs-training.md)
- **Recipe:** S1 dense vector RBS (mean pool) → S2 freeze, scalar/burden RBS
  → S3 freeze RBS, learned gene MBS hop → S4 unfreeze fine-tune. Then
  frozen **`rbs_enet` / `mbs_enet(_nested)`** as co-primary readout.
- **Why:** RBS probes and nested enet beat joint `mbs_e2e` (especially
  N-light age). LP-FT already showed freeze-then-train recovers vector.
- **Data:** encoder on nine-pack **34k**; Hub **173k** + extra packs =
  freeze-reuse (ATS, BMI, ancestry, cancer types, AD/PD/stroke, …).

### 10d — Reference checkpoint deliverable

- **Status:** `pending` (after 10e, not after joint-e2e P2-G)
- **Done when:** documented pretrained checkpoint(s) (MBS ± RBS), input/score
  contract, short association-testing note (CpG→gene multiple-testing
  reduction). See campaign Phase 4.

### Deferred to after Milestone 12 OOF

CpGPT embeddings/positional features (`use_cpgpt_static_features`, off
everywhere — a different, untried form of annotation enrichment; region-
level annotation richness was already tried and underperformed), sex-
chromosome-based sex imputation for the small `sex` pack, epigenetic-
clock-based (Horvath/CpGPT/MethylGPT) age imputation for age-unlabeled
samples. Real ideas, deliberately sequenced after the OOF architecture
decision, not before.

**Hard stop:** do not auto-launch Milestone **11** or **12** from this campaign.

---

## 11. Fold-selected panel + full model (alias: 7G′ Stage B)

- **Status:** `deferred` (parallel track — **not** a hard gate on Milestone **12**)
- **Plan:** [`plans/milestone-11-fold-selected-panel.md`](plans/milestone-11-fold-selected-panel.md)
- **Runner:** `scripts/run_7g_prime_stage_b.py`
- **Role after Milestone 10:** sparsity / fold-safe **panel product** and classical
  `C-mvalue-enetS` comparator — **not** architecture selection (that locked in **10**).
- **CPU prep:** `--panels-only` / `--classical-only` / `--folds` may continue anytime;
  orphan census under
  `reports/inspection/stage0_7g_prime_matched_probe/orphan_rbs_census.md`.
- **Preferred thinner GPU matrix (when scheduled):** panels + `C-mvalue-enetS` +
  finalists-on-S (`N-cascade-S` = P2-G params, `N-light` on same panel) ± fusion
  ablations — not the full historical Stage A arm soup.
- **Done when:** fold-safe panels, matched classical/finalist-on-S report,
  `direct_cpg.zarr` when direct loci exist; report under
  `reports/inspection/stage0_7g_prime_matched_probe/`.
- **Depends on:** honest **10** scale/architecture decisions (cleared: P2-G +
  N-light@64). Seed-mask gate cleared (9c; not adopted).
- **Does not block:** Milestone **12** finalist OOF on the gene-linked path.
- **Hard stop:** do not auto-launch full Stage B GPU from the Milestone 10 keeper;
  schedule explicitly after review.

---

## 12. Final study-grouped OOF (alias: historical Milestone 7)

- **Status:** `in_progress` — **N-light 5×6 first** (GPU 2 after dense S1);
  cascade 5×6 still waits on **10e** S1–S4 + nested-enet ranking; **does not wait on Milestone 11**
- **Plan:** [`plans/milestone-12-final-oof.md`](plans/milestone-12-final-oof.md)
  (alias stub: [`milestone-13-final-oof.md`](plans/milestone-13-final-oof.md))
- **Arms policy:** **N-light 5×6 first**, then cascade after 10e. Do **not**
  5×6 joint `mbs_e2e` as the product score.
  1. **Light (now)** — N-light@64 + **`mbs_enet_nested`**. Split
     `hub-nine-pack-5fold-v1`. GPU 2 max VRAM.
  2. **Cascade (later)** — P2-G topology trained **S1–S4**; rank vs N-light
     under nested enet (vector `gene_rho` still in play).
  Extra traits with **n≥600** (BMI, ancestry, cancer pack, AD 945) are
  **frozen-score probes**, not extra OOF arms.
- **Depends on:** N-light: (7A)–(7F), **8**, **9**, **10** N-light@64, one-hop
  smokes, 3-fold nested enet. Cascade: **10e** S1–S4 smoke. Milestone **11**
  is parallel/optional.
- **Panel for this OOF:** same **gene-linked** setting used in Milestone **10**
  refs (not the Milestone **11** fold-selected sparse panel unless a separate
  sparse-panel OOF is explicitly scoped later).
- **Done when:** OOF gene-aggregated RBS / MBS (+ orphan RBS + direct as
  applicable), age and tissue predictions, leakage controls, orientation-aligned
  scores (ADR 0008), no TBS (ADR 0009). Protocol: **5** outer folds × up to
  **6** restarts **per finalist**.
- **Depends on:** (7A)–(7F), **8**, **9**, **10** topology screen, **10e**
  staged-recipe smoke, **10a++** one-hop smokes. Milestone **11** is
  parallel/optional for this OOF.
- **Note:** 3-fold / 1-restart plumbing smoke allowed; must not overwrite v0.1
  freezes ([ADR 0007](adr/0007-crossfit-prerequisites.md)).
- **Next:** Milestone **13** expression continue/finetune (optional path after
  OOF scores exist).

---

## 13. Expression auxiliary (alias: 7G″)

- **Status:** `deferred` — runs **after** Milestone **12** OOF
- **Plan:** [`plans/milestone-13-expression-auxiliary.md`](plans/milestone-13-expression-auxiliary.md)
  → detail [`plans/milestone-7g-double-prime-expression-auxiliary.md`](plans/milestone-7g-double-prime-expression-auxiliary.md)
- **Intent:** either **continue training** the OOF / pretrained methylation
  encoder or **finetune** it with gene-expression prediction as auxiliary (or
  primary) supervision — not a random-split TCGA copy of RSMethy-Net.
- **Data prerequisite (in plan):** download and catalog matched
  methylation–expression cohorts (RNA-seq / microarray) with study-grouped
  sample overlap; persist under `$MBS_DATA_ROOT` with manifests — **do not**
  start expression GPU work until download + overlap census land.
- **Depends on:** Milestone **12** OOF checkpoint(s) (or an explicit interim
  pretrained checkpoint from **10d** if OOF is postponed by ADR).
- **Not a gate** for **10** / **11** / **12**.

---

## 14. Optional layers (after core OOF is stable)

- **Status:** `deferred`
- **Rule:** Do not start until Milestone **12** OOF produces a real score
  pipeline (milestones **1–12**). Expression continue/finetune is **13**, not
  this section. Graph-layer cCRE for scoring is **7C/7F**. Tile **scores** are
  out (7F); leftover CpGs are direct. Full vision context:
  [`STRATEGIC_PLAN.md`](STRATEGIC_PLAN.md).

### Deferred candidates (only these)

| ID | Candidate | Intent / acceptance hint |
|----|-----------|--------------------------|
| **a** | EWAS Atlas enrichment | Compare significant gene–trait hits to Atlas curated associations / pathway enrichment |
| **b** | GEO-enriched training release | Phenotype `deepmat-data-geo-dev-v1` **built**; matrix convert + train launch still gated |
| **c** | Epivariants / episignatures | Explicit epivariant calling and clinical episignature work |
| **d** | Learned ProbeNormalizer (Level 2) | Bounded residual adapter; fold-fitted; after Level-1 |
| **e** | PROTRIDER-style / masked AE (Level 3) | Only if 7D/7E show Level-1 insufficient; Student-t or masked recon + phenotype; not default |
| **f** | ONT/PacBio long-read methylation | Needs new ingestion path — out of Milestone 10 scope |

Removed from this list (not Stage 0 optional gates here): ComBat-met, TileDB/Zarr
WGBS benchmark, ClickHouse, MethylGPT priors / richer FM fusion.



## Agent checklist

Before claiming a milestone `done`:

1. Acceptance criteria above are met with evidence (tests, report paths, artifact
   manifests).
2. Status in this file is updated in the same change set.
3. Required checks from `AGENTS.md` pass for code changes.
4. Do not mark complete because a stub module or download script exists.
5. Catalog CLI (`refresh-release`, `validate-release`, `phenotype-census`,
   `trait-eligibility`) **is implemented**. Do not advertise *other* placeholder
   commands as implemented. Default census dir is
   `reports/inspection/deepmat-data-v1` (hyphen); keep the committed snapshot
   under `reports/inspection/deepmat_data_v1/` (underscore).