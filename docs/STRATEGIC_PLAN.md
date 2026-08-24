# Strategic plan: epigenome-wide data integration and deep set modeling

Project-owned distillation of the methyl-burden-score strategic implementation
vision. Normative Stage 0 engineering rules remain in [`AGENTS.md`](../AGENTS.md),
[`ARCHITECTURE.md`](ARCHITECTURE.md), and [`TODO_PIPELINE.md`](TODO_PIPELINE.md).
Primary open data source decision: [`adr/0002-ewas-datahub-primary-source.md`](adr/0002-ewas-datahub-primary-source.md).
Post-v0 sequencing: [`adr/0007-crossfit-prerequisites.md`](adr/0007-crossfit-prerequisites.md);
programme brief: [`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).

## Motivation

Rare-variant burden methods (especially DeepRVAT) aggregate sparse, variable
variant sets into a gene-level impairment score, reducing multiple-testing
burden and capturing combined locus effects. Methyl burden scores apply the same
idea to DNA methylation: map variable observed CpG sets into typed regulatory
regions, then to one scalar score per sample and gene for downstream association
and prediction. deepMAT extends this with optional non-gene regulatory (RBS),
intergenic tile (TBS), and direct CpG paths ([ADR 0006](adr/0006-multipath-noncoding-scores.md)).

## Data strategy

Prefer open, curated CNCB EWAS Open Platform resources over paywalled corpora:

| Resource | Role |
|----------|------|
| **EWAS Data Hub** | Primary Stage 0 pilot matrices and open-scale GMQN-normalized baselines (`EWAS_db/` per-study betas; `download/*_v1.zip` tissue/blood/brain/disease/age packs). |
| **EWAS Atlas** | Curated association knowledge for post–Stage 0 enrichment / calibration checks. |
| **CpGCorpus** | Optional alternate (requester-pays S3). Historical milestone-1 inspection used `GSE125367`; not the default pilot path going forward. |

Inventory and download commands: [`EWAS_DATA.md`](EWAS_DATA.md). Labeling GSE
notes and the CpGCorpus alternate path: [`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md).

Manuscript FTP names such as `tissue_methylation.zip` correspond to Data Hub
archives under `download/` with `_v1` suffixes (for example
`tissue_methylation_v1.zip`). Prefer HTTP mirrors via
`make download-ewas-datahub` / `scripts/download_ewas_datahub.sh` on this host.

### Storage layers ([ADR 0005](adr/0005-catalog-matrix-independence.md))

```text
DuckDB + Parquet   metadata, phenotypes, provenance, splits, eligibility
Zarr               dense sample × locus / score matrices
```

Keep the catalog independent of the matrix-store implementation. Do **not** add
ClickHouse for Stage 0. Benchmark TileDB sparse (or sharded Zarr v3) only when
the first representative WGBS cohort arrives. Never materialize sample×CpG longs
in DuckDB.

## Stage 0 architecture (current order)

Stage 0 does **not** require a PROTRIDER-style autoencoder, ComBat-met, or
REGENIE as defaults. The critical path is:

1. Inspect one small real source (done; historical CpGCorpus evidence).
2. Canonical annotation graph (simple regions).
3. Static locus features (CpGPT default).
4. One pilot canonical matrix from **EWAS Data Hub**.
5. Flat DeepRVAT-style max-pooling baseline.
5b–5d. Phenotype registry, Hub pack matrices, multitask, max-N age/tissue/sex.
6. Hierarchical region model (frozen v0.1 residual baseline).
7A. Harmonized data release + phenotype census (**done**).
7B. Complete nine-pack canonical matrices (**current gate**).
7C. Architecture corrections (multi-path scores, splits, heads).
7D. Fold-fitted Level-1 normalization.
7E. Development CV (architecture selection).
7. Final study-grouped OOF cross-fitting (blocked until 7A–7E).

Model contracts: [`ARCHITECTURE.md`](ARCHITECTURE.md). Milestone checklist:
[`TODO_PIPELINE.md`](TODO_PIPELINE.md).

Frozen references: **deepMAT-flat-v0.1**, **deepMAT-hierarchical-v0.1**,
**deepmat-data-age-tissue-sex-v1**. Hierarchical v0.1 underperformed flat on
tissue accuracy and age MAE; it remains a baseline, not the preferred phenotype
model.

## Post–Stage 0 multimodal stack

After Milestone 7 produces OOF scores, the longer-term vision adds further
layers. Do not start these while 7A–7E / 7 are open (see
`.cursor/rules/pipeline-todo.mdc`).

### Module A — Normalization and epimutation features

**Stage 0 next (7D):** fold-fitted robust per-CpG M-deviation channels (Level
1), trained on training studies only. Hub GMQN betas remain the canonical raw
matrix.

**Later ablation (not default):** learned ProbeNormalizer (Level 2); masked /
PROTRIDER-style autoencoder (Level 3) only if Level 1 is insufficient on held-
out phenotype and stability metrics. Do not select on reconstruction loss alone.

### Module B — Contextual epigenetic embedding

Enrich scalar betas (and optional epimutation scores) with static sequence /
position / regulatory context drawn from foundation-model practice (CpGPT
sequence-adapter default in Stage 0; MethylGPT priors as ablation). Full DNA-LM
fusion and ChromHMM-rich embeddings remain optional later.

### Module C — Deep Set aggregation

Already the Stage 0 core: shared φ per CpG, permutation-invariant pooling
(max), shared ρ → sigmoid MBS in `[0, 1]`. Hierarchical CpG→region→gene is the
Stage 0 upgrade over flat DeepRVAT-style pooling. Milestone **7C** adds RBS /
TBS / direct CpG paths so ~30% unassigned loci are not compressed to one scalar.

### Association testing — REGENIE

Export gene-level MBS as genetic pseudodosages (scale by 2.0 for diploid-like
dosage in `[0, 2]`), write BGEN/VCF, then REGENIE Step 1 (polygenic background
on common SNPs) and Step 2 (test burden scores; Firth for imbalanced binary
traits). Validation against EWAS Atlas enrichments is deferred with this layer.

## Engineering concerns (deferred ingest / custom cohorts)

- **GMQN:** Data Hub baselines are already GMQN-normalized; do not re-invent
  reference fitting for Hub-ingested matrices.
- **ComBat-met:** Beta-regression batch correction for user-supplied IDAT /
  uncorrected cohorts via an R bridge later; not required for Hub GMQN data.
- **EPICv2 IlmnID duplicates:** Explicit mean/robust-mean collapse with recorded
  probe IDs is Milestone **7B** (not a forever-deferred §8 item).
- **Array harmonization:** Manifest parsing for 450K / EPIC v1 / EPIC v2 aligns
  with milestone 2 annotation work using vendor references, not a greenfield
  Bioconductor-only stack.

## Compatibility with this repository

| Manuscript / greenfield prompt assumption | This repository |
|-------------------------------------------|-----------------|
| Conda, Python 3.10, `environment.yml` | `uv`, Python ≥3.11 (`AGENTS.md`) |
| `src/data_engineering/`, `src/models/` | Package under `src/mbs/` |
| New FTP-only downloader from scratch | Existing HTTP downloaders + Makefile targets |
| Implement AE → Deep Set → REGENIE before a pilot matrix | Stage 0 order above; AE/REGENIE deferred |

## Appendix: manuscript agent-prompt phases

An external strategic manuscript listed Cursor-style prompt phases (scaffold,
FTP, ComBat-met, autoencoder, Deep Set Lightning module, BGEN/REGENIE). Those
phases are **historical inspiration only**. They must not be executed as a
greenfield rewrite: they conflict with existing ADRs, package layout, and the
Stage 0 milestone order. Map useful ideas into deferred todos in
[`TODO_PIPELINE.md`](TODO_PIPELINE.md) section 8 and implement inside `src/mbs/`
when Stage 0 is complete.
