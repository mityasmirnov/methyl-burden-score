# EWAS DataHub sample-info packs (unzipped)

Unpacked from `data/raw/ewas_datahub/download/sample_*.zip`.
Zips are **deleted after extract**; re-download from
https://download.cncb.ac.cn/ewas/datahub/download/ if archives are needed again.

Each pack contains:
- `sample_*.txt` — tabular sample metadata (Cursor-indexed)
- `sample_*.RData` — same info as serialized R (ignored by Cursor)

Present families: age, blood, brain, disease, sex, tissue, ancestry (`sample_race.txt`), bmi, cancer.

Note: the ancestry pack zip is named `sample_ancestry_category_*` but the
member files are `sample_race.txt` / `sample_race.RData`.

Regenerate structure report: `uv run mbs inspect ewas-metadata`
