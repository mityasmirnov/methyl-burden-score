#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${MBS_ROOT:-${MBS_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}}"
cd "$PROJECT_ROOT"

if [[ -n "$(git status --porcelain)" ]]; then
  printf 'ERROR: working tree must be clean before adding submodules.\n' >&2
  exit 1
fi

mkdir -p vendor

add_submodule() {
  local url="$1"
  local path="$2"

  if [[ -e "$path/.git" || -f "$path/.git" ]]; then
    printf 'Already present: %s\n' "$path"
    return
  fi

  git submodule add "$url" "$path"
}

# Canonical upstream repositories. Use forks only when a project-specific patch
# is required and record that decision in vendor/SOURCES.lock.yaml.
add_submodule https://github.com/PMBio/deeprvat.git vendor/deeprvat
add_submodule https://github.com/lucascamillomd/CpGPT.git vendor/cpgpt
add_submodule https://github.com/albert-ying/MethylGPT.git vendor/methylgpt
add_submodule https://github.com/Christensen-Lab-Dartmouth/MethylCapsNet.git vendor/methylcapsnet
add_submodule https://github.com/adamewing/methylartist.git vendor/methylartist
add_submodule https://github.com/bethan-mallabar-rimmer/EPICv2_manifest.git vendor/epicv2_manifest
add_submodule https://github.com/zhou-lab/InfiniumAnnotation.git vendor/infinium_annotation

git submodule update --init --recursive

git add .gitmodules vendor/SOURCES.lock.yaml vendor/README.md
git submodule status

printf '\nReview exact SHAs above, then commit and push if needed.\n'
