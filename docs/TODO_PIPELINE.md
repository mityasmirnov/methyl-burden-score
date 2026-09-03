# Stage 0 scientific pipeline TODO

Authoritative milestone list for coding agents. Update status here when a
milestone is **truly done** (acceptance criteria met), not when scaffolding
exists. Cursor agents must read this file at session start and after finishing
work; see `.cursor/rules/pipeline-todo.mdc`.

Status values: `done` | `in_progress` | `pending` | `deferred`

True next milestone after bootstrap:

> … → max-N flat age/tissue/sex (5d — **done**) → hierarchical residual
> baseline (6 — **done**) → harmonized release + phenotype census (7A — **done**) →
> nine-pack matrices (7B — **done**) → architecture corrections (7C —
> **done**, fixture + Hub DeepRVAT A/B) → development CV (7E — **done**) →
> Hub multitask + hygiene (**7E′** — **done**) → **RBS→gene + direct topology
> (7F — done)** → **methylation-only full eval (7G — done)** →
> **Phase-2 gene-only cascade grid (7G′ Stage A — DeepRVAT screen
> in progress)** → **fold-selected panel + full model (7G′ Stage B — plumbing
> done, GPU pending)** → **final OOF cross-fitting (7)** → one score matrix.

**Current gate:** **7G′ Stage A DeepRVAT screen** (mixed pooling, RBS
diagnostic, vector cascade, one-hop). Required P2/P4/P5-max/`C-mvalue-*-G`
GPU arms already landed with test-only **`mbs_e2e`**; cascade does **not** beat
`C-mvalue-enet-G` by ≥0.03 F1. **P5 inactive.** Stage B GPU follows after the
screen selects (or rejects) gene aggregation. Plans:
[`plans/milestone-7g-prime-stage-a-deeprvat-screen.md`](plans/milestone-7g-prime-stage-a-deeprvat-screen.md),
[`plans/milestone-7g-prime-matched-probe-lightweight.md`](plans/milestone-7g-prime-matched-probe-lightweight.md).

| Layer | Status |
|-------|--------|
| Test-only `mbs_e2e` + lock refusal | **done** (code + unit tests) |
| `explicit_only` allocation + manifest | **done** (ADR 0010; YAML `*-explicit` run IDs) |
| Stage B selector / panel artifact / honest fusion names | **done** (code only) |
| 7G″ expression plan | **done** (deferred; no training code) |
| Stage A required GPU arms (P2/P4/P5-max/`C-*-G`) | **done** (P5 inactive thereafter) |
| Stage A DeepRVAT screen | **in_progress** ← **current ops gate** |
| Stage B GPU run | **pending** (after screen) |
| Milestone **7** 5×6 OOF | **blocked** |

**Trustworthy numbers (`explicit_only` panel, 51 375 gene-linked CpGs, test split):**

| Question | Best arm | Tissue macro-F1 | Caveat |
|----------|----------|----------------:|--------|
| Classical on gene panel | `C-mvalue-enet-G` | **0.388** (±0.018) | Still the tissue leader |
| Best MBS readout | `P2-G` **`mbs_enet`** | 0.385 (±0.053) | Same encoder as e2e; elastic-net heads |
| Best cascade product path | `P2-G` **`mbs_e2e`** | 0.373 (±0.038) | max/max, 15 epochs; P4 mean 0.370 |
| P4 mean/mean e2e | `P4-G` | 0.370 (±0.059) | Tied with P2 within noise |
| P5 longer max train | `P5-G-max` | 0.356 (±0.042) | Did not help |
| Invalid (do not cite) | pre-fix `mbs_e2e` on P*-G | ~0.67–0.70 | train+val+test leak |
| Invalid (do not cite) | pre-fix `N-light-gene-*` **`mbs_e2e`** | ~0.000–0.001 | orientation anchor + head/score mismatch; **linear/enet probes valid** |

**N-light repair (2026-03):** orientation contract **v2** (manifest-only flip;
eval applies `1-MBS` pairing for legacy negated heads). Re-eval without retrain:
`scripts/reeval_7g_light_stage_a.py`. One-hop is **not** rejected — frozen MBS
probes show representation signal. Diagnostic fold-0 mean run:
`configs/experiment/stage0_7g_gene_only_probe_light_mean_diag.yaml`.

**Static annotation channels (2026-09):** `cpg_context` (UCSC CGI, Milestone 3
artifact at `data/canonical/annotations/loci.parquet`) is now wired into
`build_gene_cols()` and `build_flat_region_gene_index()`. Regulatory channels
(cCRE/DHS/ChromHMM) remain zero — source files not on disk (Stage A non-goal).
Versioned annotation artifacts: `scripts/build_stage_a_locus_annotations.py`
→ `reports/inspection/stage0_7g_gene_only_probe/locus_annotations/`.
Ablation grid (fold 0, 8 ep, two seeds):
A0 `m_only` | A1 `m_role` | A2 `m_context` | A3 `m_role_context` | A4/A7 `full` |
N0 `obs_only` | N1 `anno_only` | N2 reg-permuted | N3 reg-zero.
Run all arms with `--fold 0 --device cuda` via `run_7g_gene_only_probe.py`.

**Provisional lock:** `P2-G` max/max, 15 epochs (`lock_recommendation.json`);
cascade is **not** clearly ahead of classical. Stage A screen may revise the
topology before Stage B. **7G″** expression pilot is **deferred**. Final
Milestone **7** OOF remains blocked.

Runners: Stage A `scripts/run_7g_gene_only_probe.py` · Stage B
`scripts/run_7g_prime_stage_b.py` (both `--device cuda` on GPU hosts).

**7G** methylation eval and **tissue probe P0–P3** are done (historical evidence
only — P2 ~0.38 used late fusion, not MBS-only). Plan:
[`plans/milestone-7g-cascade-tissue-investigation.md`](plans/milestone-7g-cascade-tissue-investigation.md).
Master brief:
[`plans/milestone-7g-prime-matched-probe-lightweight.md`](plans/milestone-7g-prime-matched-probe-lightweight.md).

**7F** (RBS→gene + direct leftover, no TBS) is **done**
(`reports/inspection/stage0_7f_rbs_gene_direct/`; [ADR 0009](adr/0009-drop-tbs-scores.md)).
Product scores remain gene-aggregated RBS (MBS), qualified per-region orphan
RBS, and indexed direct CpGs — **no TBS**. The present 7F
`direct_contrib.zarr` is a phenotype diagnostic (one task prediction per sample).
Association export requires `direct_cpg.zarr` (7G′ Stage B).

**ATS** = Age/Tissue/Sex frozen Hub GSM-union
`matrix-hub-age-tissue-sex-full-v1` (13 548). Graph-v2 is **on disk**
(`graph-grch38-gencode38-cgi-tile-v2`; inspection
`reports/inspection/annotation_graph_cgi_tile_v2/`; plan
[`plans/milestone-7c-graph-v2-topology.md`](plans/milestone-7c-graph-v2-topology.md)).
7E selected `N-multipath-l1a` under a 2-epoch / 8 192-locus budget and
linear fusion of **region means**, not neural scores — that is evidence, not
the shipped topology (7F replaces it). 7E′ hygiene is **done**. Readiness:
[`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md)
(§7G); 7F report
[`reports/inspection/stage0_7f_rbs_gene_direct/analysis.md`](../reports/inspection/stage0_7f_rbs_gene_direct/analysis.md).

**7A–7D** are `done` (7C = fixture + Hub smoke). Hub nine packs and 7B full
matrices are in the live DuckDB release (refresh 2026-09-02:
**149 244** samples, **1 584** studies, **216 476** phenotype rows, **21**
matrix artifacts including all nine Hub full packs). EWAS_db is incomplete
(**1 353**/1 989 studies, **132 289** GSM files) and **not** a gate. Authoritative
census: `reports/inspection/deepmat_data_v1/` (underscore; matches this
refresh). Ignore `reports/inspection/deepmat-data-v1/` if it shows ~5 GSM
(fixture leak into the CLI default hyphen path).
**Do not retrain v0.1** or start Milestone **7** until **7G′ Stage A and B** land
(7A–7E′ are already `done`; [ADR 0007](adr/0007-crossfit-prerequisites.md),
[ADR 0008](adr/0008-score-identifiability.md)). Programme brief:
[`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md)
(glossary: ADR / MBS / RBS / leftover-direct / Level-1 MAD / 3×2 / 5×6 OOF;
when training runs; DeepRVAT-style joint aggregation + linear heads).

**Training schedule:** 7A–7B data only; 7C–7D code/fixtures/smokes; **7E** =
architecture selection on frozen ATS (3×2, **done**); **7E′** = Hub multitask
+ hygiene (**done**); **7F** = RBS→gene cascade + direct leftover (no TBS
scores); **7G** = methylation-only re-eval that closes 7E evaluation gaps;
**7G′** = gene-only architecture selection + fold-selected panel / full model;
**7** = final 5×6 OOF after 7G′. Neural arms train shared score aggregation
**and** linear phenotype heads end-to-end (DeepRVAT pattern).
Architecture comparison tables use **methylation-input methods only**.
Metadata-only (study + platform IDs, no betas) is a leakage alarm from 7E′,
not a competitor.

Frozen references (do not overwrite): **deepMAT-flat-v0.1** /
**deepMAT-hierarchical-v0.1** / **deepmat-data-age-tissue-sex-v1**. Hierarchical
v0.1 is a valid baseline, not the preferred phenotype model. Residual-only
near-chance on an ordered 512-sample prefix is **not** evidence that noncoding
CpGs lack signal.

Hub **disease** profile zip is complete (2026-08-11). `EWAS_db` All-Data
mirror remains in progress (**1 353** local studies / 1 989 advertised;
**132 289** GSM files; failure manifest + post-download hook — see
[`EWAS_DATA.md`](EWAS_DATA.md)) and is **not** a gate for 7G′ or Milestone 7;
re-run `make catalog-refresh-release` as more `EWAS_db/{GSE}/` dirs arrive.

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
  **GEO sample backfill (pilot):** `make fetch-geo-sample-metadata` then
  `make catalog-refresh-release` (see [`plans/geo-metadata-backfill-ewas-db.md`](plans/geo-metadata-backfill-ewas-db.md)).
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
  Report selects architecture for Milestone 7. Eval-time branch masking and
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

---

## 7G. Methylation-only full evaluation (close 7E gaps)

- **Status:** `done`
- **Impl plan:** [`plans/milestone-7g-methylation-eval.md`](plans/milestone-7g-methylation-eval.md).
- **Why:** The 90-cell 7E bake-off **finished**. Gaps were evaluation quality,
  not a crashed trainer: 2 epochs and 8 192 / 482 379 loci; late fusion was
  not neural MBS; T-mean-region was not a named cell; HGB stood in for
  LightGBM; neural AUROC was a binary helper; PCA-SVA not Bioconductor sva;
  sex incomplete in the merged dump. Do **not** conclude “trees beat Deep
  Sets” from that table.
- **Comparators (methylation input only):** M-value ridge, penalised
  (elastic-net / SGD), histogram gradient boosting or LightGBM, optional
  train-fold PCA-SVA then linear. Same frozen 7E studies/folds. **Do not**
  put metadata-only (study + platform, no methylation) in architecture
  ranking tables. That control remains a 7E′ leakage alarm only.
- **Done when:**
  - 7F cascade retrained on the **same** frozen ATS folds with a documented
    budget above 7E’s 2-epoch / 8 192-locus ceiling (or an explicit remaining
    ceiling in the report)
  - Classical M-value models use methylation matrices only; ROC for **sex**
    and **tissue one-vs-rest** from the neural fusion scores (age stays MAE /
    R², no ROC)
  - Sex metrics present in the summary dump; region-mean transparent arm
    named if kept
  - Report under `reports/inspection/stage0_7g_methylation_eval/` names the
    Milestone 7 topology from methylation-input methods only
- **Depends on:** (7F).
- **Plan:** [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md)
  (§7G).
- **Report:** [`reports/inspection/stage0_7g_methylation_eval/analysis.md`](../reports/inspection/stage0_7g_methylation_eval/analysis.md)
  — ranking winner **`C-mvalue-enet`** (tissue macro-F1); **7F cascade** weak on
  tissue (~0.09 F1) at 65k/15ep. Follow-up:
  [`plans/milestone-7g-cascade-tissue-investigation.md`](plans/milestone-7g-cascade-tissue-investigation.md).
- **Tissue probe (P0–P3):** `done` (historical) — report
  [`reports/inspection/stage0_7g_cascade_tissue_probe/`](../reports/inspection/stage0_7g_cascade_tissue_probe/analysis.md).
  P0 ~0.09 F1 (late fusion); P2 ~0.38 F1 (**MBS train, late-fusion test** — not
  a clean MBS-only vs enet comparison). **Corrected gene-only grid (P2-G …
  P5-G, `C-mvalue-enet-G`)** scaffolding is **done**; honest GPU rerun on
  **`explicit_only`** is **7G′ Stage A** (**done** for required arms). Then
  **7G′ Stage B** before Milestone **7** OOF.

---

## 7G′. Gene-only architecture selection + matched-panel benchmark

- **Status:** scaffolding **`done`** · required Stage A GPU arms **`done`** ·
  DeepRVAT screen **`in_progress`** · Stage B code **`done`** / GPU **`pending`**
- **Plan:** [`plans/milestone-7g-prime-matched-probe-lightweight.md`](plans/milestone-7g-prime-matched-probe-lightweight.md)
- **Screen:** [`plans/milestone-7g-prime-stage-a-deeprvat-screen.md`](plans/milestone-7g-prime-stage-a-deeprvat-screen.md)
- **Runner (Stage A):** `scripts/run_7g_gene_only_probe.py` · background:
  `scripts/train_7g_gene_only_probe_background.sh`
- **Runner (Stage B):** `scripts/run_7g_prime_stage_b.py` · background:
  `scripts/train_7g_prime_stage_b_background.sh`

### Finished (code + docs — `4f5e022`)

- Test-only **`mbs_e2e`** (`eval_split=test`, `n_eval_samples`); regression tests
- Report writers refuse lock without valid e2e + completed **`C-mvalue-*-G`**
- [ADR 0010](adr/0010-gene-allocation-policy.md): **`explicit_only`** /
  `bounded_nearest` / `legacy_nearest`
- Stage B: study-grouped multitask panel selector, `fold_panels/fold_*_panel.json`,
  Stage A lock → training params, post-hoc fusion arm renames
- [7G″ plan](plans/milestone-7g-double-prime-expression-auxiliary.md) (deferred)

### Stage A — gene-only MBS architecture selection

- **Required arms done when:** Report under
  `reports/inspection/stage0_7g_gene_only_probe/` with per-fold tables,
  **test-only** **`mbs_e2e`** (`eval_split=test`), provisional lock, orphan
  fusion decision; arms `P2-G`, `P4-G`, `P5-G-max` (**inactive thereafter**),
  `C-mvalue-*-G`, optional `P2-orphan-ablation` on **`explicit_only`** panel.
- **Evidence (required arms):**
  - **`P2-G-explicit` / `P4-G-explicit` / `P5-G-max-explicit`:** 3/3 folds,
    `eval_split=test`. Best cascade **`mbs_e2e`:** P2-G **0.373**; P4-G 0.370;
    P5-G-max 0.356 (**P5 inactive** — do not extend).
  - **`mbs_enet`** on frozen P2/P4 MBS: P2 **0.385** (closest neural readout to classical).
  - **`C-mvalue-classical-G`** on **`explicit_only`** (51 375 cols): enet-G **0.388**
    still leads tissue F1; cascade not ≥0.03 ahead.
  - **`P2-orphan-ablation`:** orphan RBS negligible (Δ F1 ≈ 0)
  - **Provisional lock:** `P2-G` max/max, 15 epochs (`lock_recommendation.json`)
- **DeepRVAT screen (current):** mixed scalar pooling, RBS-only diagnostic,
  vector-region cascade, one-hop annotated DeepSet; always report tissue/age/sex;
  Tier 1 (5 ep) then Tier 2 (15 ep) for promoted arms. See screen plan.

### Stage B — fold-selected panel + full model

- **Done when:** Report under
  `reports/inspection/stage0_7g_prime_matched_probe/` with fold-safe
  `C-mvalue-enetS`, `N-cascade-S`, **`N-light-type`**, post-hoc
  **`N-mbs-posthoc-full-fusion`** / **`N-mbs-posthoc-mbs-direct`**, canonical
  `fold_panels/fold_*_panel.json`, and **`direct_cpg.zarr`** when direct loci exist.
- **Evidence:** selector + runner plumbing landed; **GPU run pending** (after
  Stage A screen). Seed-gene transfer is a **separate** design (screen plan §10).

### 7G″ — expression auxiliary (deferred)

- **Plan:** [`plans/milestone-7g-double-prime-expression-auxiliary.md`](plans/milestone-7g-double-prime-expression-auxiliary.md)
- **Not a gate** for 7G′ or Milestone **7**; run after Stage B if pursued.

- **Depends on:** (7G); Stage A lock for pooling/epochs.
- **Blocks:** Milestone **7** 5×6 OOF.

---

## 7. Run study-grouped cross-fitting (final OOF)

- **Status:** `pending` (**blocked until 7G′ Stage A and B**)
- **Done when:** Out-of-fold gene-aggregated RBS (and genome-wide RBS + direct
  contributions), age predictions, and tissue predictions are generated with
  leakage controls (no sample/donor/replicate/held-out study scored by a
  model that saw it). **No TBS scores.** Scores are **orientation-aligned**
  (ADR 0008) before averaging. Persist fold-specific normalization,
  presence/count/`norm_present` masks, complete model lineage. Protocol: 5
  outer folds × up to 6 restarts.
- **Depends on:** (7A)–(7E′), **7F**, **7G**, **7G′**; Hub multitask hygiene in 7E′.
- **Note:** A 3-fold / 1-restart smoke of *existing* machinery is allowed for
  plumbing; it does not complete this milestone and must not overwrite v0.1
  freezes ([ADR 0007](adr/0007-crossfit-prerequisites.md)).

---

## 8. Optional layers (after core pipeline is stable)

- **Status:** `deferred`
- **Rule:** Do not start these until milestones 1–7 produce a real OOF model
  pipeline (7A–7G then 7). Full vision: [`STRATEGIC_PLAN.md`](STRATEGIC_PLAN.md).
  Graph-layer cCRE for scoring is **7C/7F**, not this section. Tile **scores**
  are out (7F); leftover CpGs are direct.

### Stage 1+ roadmap (deferred candidates)

| Candidate | Intent / acceptance hint |
|-----------|--------------------------|
| PROTRIDER-style / masked AE (Level 3) | Only if 7D/7E show Level-1 insufficient; Student-t or masked recon + phenotype; not default |
| Learned ProbeNormalizer (Level 2) | Bounded residual adapter; fold-fitted; after Level-1 |
| ComBat-met (rpy2) | Beta-regression batch correction for user/custom IDAT or uncorrected cohorts; assert corrected betas stay in `[0, 1]`. Not required for Data Hub GMQN baselines |
| TileDB sparse / Zarr v3 WGBS benchmark | First representative WGBS cohort; catalog stays on DuckDB+Parquet ([ADR 0005](adr/0005-catalog-matrix-independence.md)) |
| ClickHouse | Only if multi-user OLAP portal needed; not training I/O |
| EWAS Atlas enrichment | Compare significant gene–trait hits to Atlas curated associations / pathway enrichment |
| MethylGPT priors / richer FM fusion | Ablation after multi-path scores stable |
| Epivariants / episignatures | Explicit epivariant calling and clinical episignature work |

---

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