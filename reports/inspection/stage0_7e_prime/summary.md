# Milestone 7E′ inspection summary

Generated after Hub virtual multitask + hygiene landed and Hub-wide train
`stage0-flat-hub-multitask-v1` wrote `metrics.json`.

## Artifacts

| Artifact | Status |
|----------|--------|
| `matrix-hub-nine-pack-virtual-v1` | **34 234** GSM × **482 379** loci (route + indices; no dense Zarr) |
| `sample_phenotype_table_hub_nine_pack_v1.parquet` | age 10 002 / tissue 13 457 / sex 12 445 / disease pack 12 218 / cancer pack 10 101 |
| Split `hub-nine-pack-full-auto-v1` | 24 010 / 5 576 / 4 648 (regenerate full GSM lists via `scripts/report_7e_prime_metadata_only.py`) |
| Metadata-only ceiling | [`metadata_only.md`](metadata_only.md) (ATS freeze + Hub split) |
| Smoke run | `artifacts/runs/stage0-flat-hub-multitask-smoke-v1/` (5 heads; disease/cancer AUROC present) |
| Hub-wide run | `artifacts/runs/stage0-flat-hub-multitask-v1/` — 2 epochs, `--max-loci 8192`, full cohort; disease/cancer label lists + `val_auroc`; metadata-only sidecar in `metrics.json` |

## Hub-wide run notes

- Matrix: `matrix-hub-nine-pack-virtual-v1` (virtual multi-store)
- Split: `hub-nine-pack-full-auto-v1` (does **not** overwrite ATS freeze)
- Heads: age / tissue / sex / disease / cancer; `lambda_blood: 0`; `empty_as_control: false`
- Disease labels scored in training: 30; cancer: 43
- Best epoch: 2; device `cuda:0`
- Uncapped (full loci, 40 epochs): `scripts/train_hub_multitask_7e_prime.sh`

## Locked layout

Virtual multi-store via pack priority
`age > tissue > sex > disease > cancer > blood > brain > bmi > ancestry`.
ATS freeze and v0.1 untouched. Blood `cell_component` is **not** a pack-wide head.

## Hygiene

- Catalog: `450K` → `HM450` via `normalize_platform` (catalog refresh + future Hub pack manifests)
- Census: donor/replicate counts when present; age/BMI-by-study views
- CLI tests use temp `--report-dir` (no hyphen-default clobber)
- `*.RData` gitignored (including `reports/inspection/ewas_datahub_samples/**/*.RData`)

## Commands

```bash
mbs matrix build-hub-virtual
mbs phenotypes build-hub-union-table
uv run python scripts/report_7e_prime_metadata_only.py
scripts/train_hub_multitask_7e_prime.sh   # or: mbs train flat --config configs/experiment/stage0_flat_hub_multitask_v1.yaml --run-id stage0-flat-hub-multitask-v1
```
