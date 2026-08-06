# Reference implementations

This directory contains read-only Git submodules after running:

```bash
source scripts/activate_data_environment.sh
make references
```

Expected submodules:

```text
vendor/deeprvat
vendor/cpgpt
vendor/methylgpt
vendor/methylcapsnet
vendor/methylartist
vendor/epicv2_manifest
vendor/infinium_annotation
```

## Policy

- Reference repositories are recorded only as Git submodule **gitlinks** (index mode `160000`) plus `.gitmodules` URLs and the pins in `SOURCES.lock.yaml`.
- After `git submodule update --init`, full working trees appear on disk under `vendor/<name>/`. That checkout is expected. Do **not** copy upstream trees into the repo as ordinary tracked directories (`cp -r`, vendored tarballs, or un-submoduled clones).
- Agents must not commit files inside submodule working trees from this repository; patches belong in the upstream fork and are pinned here by commit SHA.
- Submodule commits are pinned and reviewed before every reported experiment.
- Project code must not import from `vendor/` at runtime.
- Foundation-model export tools use separate environments and consume a pinned submodule plus a pinned checkpoint.
- Changes to a reference project belong in its own fork, not inside this repository.
- Dataset and checkpoint licenses are reviewed separately from source-code licenses.
- Downloaded research data and Zenodo artifacts belong under `$MBS_DATA_ROOT/raw/...`, never under `vendor/`.
- Cursor indexing of bulky submodule paths is limited via `.cursorignore`; prefer targeted reads over recursive vendor scans.

## Canonical upstreams

```text
DeepRVAT              https://github.com/PMBio/deeprvat
CpGPT                 https://github.com/lucascamillomd/CpGPT
MethylGPT             https://github.com/albert-ying/MethylGPT
MethylCapsNet         https://github.com/Christensen-Lab-Dartmouth/MethylCapsNet
methylartist          https://github.com/adamewing/methylartist
EPICv2_manifest       https://github.com/bethan-mallabar-rimmer/EPICv2_manifest
InfiniumAnnotation    https://github.com/zhou-lab/InfiniumAnnotation
```

Reannotated EPICv2 manifest tables are published on Zenodo ([doi:10.5281/zenodo.14933468](https://doi.org/10.5281/zenodo.14933468)); store those files under `data/raw/manifests/`. The Git submodule is the reannotation code only.

`vendor/infinium_annotation` ([zhou-lab/InfiniumAnnotation](https://github.com/zhou-lab/InfiniumAnnotation)) provides hg38 probe ordering, coordinates, design/quality masks, and KnowYourCG feature sets for HM27/HM450/EPIC/EPICv2 (and other arrays). Together with `vendor/epicv2_manifest` and the capsule/region taxonomy approach in `vendor/methylcapsnet`, it is a primary reference for Stage 0 Milestone 2 (canonical annotation graph). Do not treat the vendor trees as a runtime import path; export or convert needed tables into `$MBS_DATA_ROOT`.

The user-owned forks remain useful for experimental patches and historical work, but they should replace an upstream submodule only through an explicit architecture decision and a pinned commit.

## Installation on the server

Because GitHub's file API cannot create Git submodule objects, this step must be performed in the local clone:

```bash
cd /data/projects/methyl-burden-score
git switch main
git pull --ff-only
bash scripts/add_reference_submodules.sh
git push origin main
```

Review exact commits:

```bash
git submodule status
git diff --submodule=log --cached
```

Do not recursively index submodule datasets, checkpoints, results, or dependency caches in Cursor.