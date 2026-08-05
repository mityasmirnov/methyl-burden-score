#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MBS_PROJECT_ROOT:-/data/projects/methyl-burden-score}"
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

git submodule update --init --recursive

git add .gitmodules vendor

git commit -m "Add pinned reference implementations as submodules"

printf '\nReference repositories added. Review exact SHAs with:\n'
git submodule status
printf '\nThen push with:\n  git push origin main\n'
