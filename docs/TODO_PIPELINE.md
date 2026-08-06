# Stage 0 scientific pipeline TODO

Authoritative milestone list for coding agents. Update status here when a
milestone is **truly done** (acceptance criteria met), not when scaffolding
exists. Cursor agents must read this file at session start and after finishing
work; see `.cursor/rules/pipeline-todo.mdc`.

Status values: `done` | `in_progress` | `pending` | `deferred`

True next milestone after bootstrap:

> one source ingested cleanly → one graph built → one canonical matrix written →
> one baseline trained → one cross-fitted score matrix produced.

---

## Ops notes (plumbing)

### Done

- Project-local `/data` roots (`MBS_*`), DuckDB schema via `mbs catalog init`,
  shallow `mbs inspect source`, CpGCorpus download scripts with log tee +
  `scripts/download_cpgcorpus_background.sh` (nohup; no auto-download).
- Path / catalog / CLI / no-home tests green; vendor kept as gitlink submodules.
- ADR: `docs/adr/0001-project-local-workspace.md`.

### Blocked / caveats (do not treat as scientific blockers)

- **uv location:** Bootstrap installs/copies `uv` under `$MBS_ROOT/.tools/uv/bin`.
  Older runs may still have a copy in `$HOME/.local/bin`. Caches remain under
  `$MBS_CACHE_ROOT` (correct). Prefer
  `source scripts/activate_data_environment.sh` (puts `.tools` on `PATH` first).
- **CUDA:** Host driver may be older than the torch wheel; CPU path is fine for
  Stage 0 until training needs GPU.
- **Intentional non-goals until later milestones:** Parquet population / ingest
  into catalog tables; full annotation graph; foundation-model export/training;
  committing the shallow whole-corpus inventory under
  `reports/inspection/cpgcorpus/` (GSE/GPL report for milestone 1 is enough).

### Useful commands (already safe to re-run)

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
uv run mbs doctor --create-directories
uv run mbs catalog init
uv run mbs inspect source --source-id cpgcorpus
uv run mbs inspect cpgcorpus-gpl --gse GSE125367 --gpl GPL21145
```

### Downloads later (manual; requester-pays / large)

```bash
bash scripts/download_cpgcorpus_gse.sh
# or: bash scripts/download_cpgcorpus_background.sh gse
make download-ewas-atlas
make download-ewas-datahub
make download-manifests
```

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
- **Next action:** Milestone 2 — build the canonical annotation graph.

---

## 2. Build the canonical annotation graph

- **Status:** `pending`  ← **do this next**
- **Done when:** Stable locus registry and first
  `probe → locus → region → gene` mapping exist. Keep simple: promoter, body,
  UTR, and a few annotation flags. Do **not** expand to full MethylGPT /
  MethylCapsNet topology yet.
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

---

## 3. Export static locus features

- **Status:** `pending`
- **Done when:** Offline CpGPT sequence-adapter embedding artifact is the default
  static feature, with a complete static-feature manifest (commit, checkpoint
  hash, vocabulary/locus-table hash, dims, dtype, genome build, export command).
  MethylGPT token priors remain ablation-only.
- **Depends on:** (2) locus registry.

---

## 4. Convert one pilot matrix into canonical storage

- **Status:** `pending`
- **Done when:** One pilot source is written in project-local canonical matrix
  format; slices round-trip correctly from raw file to matrix store (checksum /
  equality checks in tests or inspection report).
- **Depends on:** (1), (2).

---

## 5. Train the flat DeepRVAT-style baseline

- **Status:** `pending`
- **Done when:** Exact CpG-to-gene max-pooling baseline runs end to end: overfits
  a tiny fixture, then trains on the pilot source. Checkpoints + resolved config
  under `$MBS_ARTIFACT_ROOT`.
- **Depends on:** (4). Model module scaffolding alone is not sufficient.

---

## 6. Add the hierarchical model

- **Status:** `pending`
- **Done when:** Region layer is trained after the flat baseline is stable;
  promoter/body (and related roles) can be compared to the flat model on the
  same pilot folds.
- **Depends on:** (5).

---

## 7. Run study-grouped cross-fitting

- **Status:** `pending`
- **Done when:** Out-of-fold MBS scores, age predictions, and tissue predictions
  are generated with leakage controls (no sample/donor/replicate/held-out study
  scored by a model that saw it). Score matrix + fold-assignment hash stored.
- **Depends on:** (5) at minimum; (6) preferred for hierarchical OOF scores.

---

## 8. Optional layers (after core pipeline is stable)

- **Status:** `deferred`
- **Candidates:** MethylGPT priors, richer MethylCapsNet-inspired annotations,
  intergenic tiles, epivariants, episignature work.
- **Rule:** Do not start these until milestones 1–7 produce a real model pipeline.

---

## Agent checklist

Before claiming a milestone `done`:

1. Acceptance criteria above are met with evidence (tests, report paths, artifact
   manifests).
2. Status in this file is updated in the same change set.
3. Required checks from `AGENTS.md` pass for code changes.
4. Do not mark complete because a stub module or download script exists.
