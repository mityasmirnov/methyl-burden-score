# deepmat-data-geo-dev-v1

- Generated: `2026-09-07T18:11:06Z`
- Mutates ATS / Hub packs: **false** / **false**

## Arm counts

| arm | n_samples |
| --- | ---: |
| `hub_only` | 34090 |
| `geo_only` | 68416 |
| `hub_geo` | 102506 |
| `metadata_only` | 102506 |

## Artifacts

- `/data/projects/methyl-burden-score/data/canonical/releases/deepmat-data-geo-dev-v1/manifest.json`
- `/data/projects/methyl-burden-score/data/canonical/releases/deepmat-data-geo-dev-v1/phenotypes/sample_info_*.parquet`
- `configs/experiment/geo_dev/*.yaml` (stubs; no training launch)

## Notes

- Phenotype-only development release (no matrix convert).
- GEO-only requires EWAS_db assay .txt on disk + acceptable age/sex/tissue.
- Disease/cancer training_ok flags are false in this builder.
- Do not enlarge matrix-hub-age-tissue-sex-full-v1.
- Larger GEO crawl and training launch remain gated.
