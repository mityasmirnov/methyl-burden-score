# GSE35069: EWAS Data Hub vs CpGCorpus

Generated: `2026-08-06T12:53:55Z`

## Sources

- Hub: `/data/projects/methyl-burden-score/data/raw/ewas_datahub/EWAS_db/GSE35069` (60 GSM*.txt)
- CpGCorpus: `/data/projects/methyl-burden-score/data/raw/cpgcorpus/GSE35069/GPL13534/betas/gse_betas.arrow` (GPL13534, 60 samples × 485577 probes)

## Sample overlap

- Shared: **60** / Hub 60 / CpGCorpus 60
- Hub-only: none
- CpGCorpus-only: none

## Probe overlap

- Shared: **485512**
- Hub-only: 0 (e.g. [])
- CpGCorpus-only: 65 (e.g. ['rs10033147', 'rs1019916', 'rs1040870', 'rs10457834', 'rs10774834'])

## Beta values (shared sample × shared probe, both finite)

- Cells compared (both finite): **28,936,798**
- Pearson r: **0.976568**
- MAE: **0.0503461**
- RMSE: **0.0861395**
- Exact match fraction: 0.0013%
- |Δ| ≤ 1e-4: 0.2036%
- |Δ| ≤ 1e-3: 2.0309%
- |Δ| ≤ 1e-2: 20.9942%
- Max |Δ|: 1 at sample `GSM861646` probe `cg09644806`
- Missingness asymmetric cells: 193,480

## Interpretation

Hub All Data is GMQN-normalized; CpGCorpus stores GEO-derived betas.
High correlation with non-zero MAE means the same study/samples/probes
are aligned but values are not bit-identical — expected across pipelines.
Do not treat CpGCorpus as the Stage 0 pilot default (ADR 0002).
