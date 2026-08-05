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
```

## Policy

- Reference repositories are never copied into this project as ordinary source directories.
- Submodule commits are pinned and reviewed before every reported experiment.
- Project code must not import from `vendor/` at runtime.
- Foundation-model export tools use separate environments and consume a pinned submodule plus a pinned checkpoint.
- Changes to a reference project belong in its own fork, not inside this repository.
- Dataset and checkpoint licenses are reviewed separately from source-code licenses.

## Canonical upstreams

```text
DeepRVAT       https://github.com/PMBio/deeprvat
CpGPT          https://github.com/lucascamillomd/CpGPT
MethylGPT      https://github.com/albert-ying/MethylGPT
MethylCapsNet  https://github.com/Christensen-Lab-Dartmouth/MethylCapsNet
methylartist   https://github.com/adamewing/methylartist
```

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