#!/usr/bin/env bash
# Mirror all public EWAS DataHub trees from
#   https://download.cncb.ac.cn/ewas/datahub/
# into $MBS_DATA_ROOT/raw/ewas_datahub/{EWAS_db,add_ewas_db,download}/
#
# Remote index children: EWAS_db/  add_ewas_db/  download/
#
# Note: nginx listings are HTML+JS. Plain recursive wget often fails to pull
# GSM*.txt under EWAS_db; use the HTML-index parser below instead.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/ewas_datahub"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET" "$LOGDIR" "$MBS_SCRATCH_ROOT/downloads"

HTTP_ROOT="https://download.cncb.ac.cn/ewas/datahub"

BASELINE_FILES=(
  tissue_methylation_v1.zip
  sample_tissue_methylation_v1.zip
  brain_methylation_v1.zip
  sample_brain_methylation_v1.zip
  blood_methylation_v1.zip
  sample_blood_methylation_v1.zip
  sex_methylation_v1.zip
  sample_sex_methylation_v1.zip
  age_methylation_v1.zip
  sample_age_methylation_v1.zip
  ancestry_category_methylation_v1.zip
  sample_ancestry_category_methylation_v1.zip
  bmi_methylation_v1.zip
  sample_bmi_methylation_v1.zip
  cancer_methylation_v1.zip
  sample_cancer_methylation_v1.zip
  disease_methylation_v1.zip
  sample_disease_methylation_v1.zip
  GMQN.zip
)

cat > "$TARGET/SOURCE.txt" <<EOF
source: EWAS DataHub @ EWAS Open Platform
http_root: ${HTTP_ROOT}/
trees: EWAS_db add_ewas_db download
policy: mirror all public DataHub trees via HTTP
ftp_note: ftp://download.big.ac.cn/ewas/datahub/ often stalls on this host; prefer HTTP
gmqn_pmid: 35069703
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

list_hrefs() {
  local url="$1"
  curl -fsSL --max-time 120 "$url" \
    | grep -oE 'href="[^"]+"' \
    | sed 's/^href="//;s/"$//' \
    | grep -Ev '^\.\./|^/|^javascript:|^#|^https?:'
}

mirror_download_zips() {
  mkdir -p "$TARGET/download"
  printf '=== Baseline/download packs -> %s/download ===\n' "$TARGET"
  for name in "${BASELINE_FILES[@]}"; do
    printf '  %s\n' "$name"
    if [[ -f "$TARGET/$name" && ! -f "$TARGET/download/$name" ]]; then
      mv "$TARGET/$name" "$TARGET/download/$name"
    fi
    wget -c -O "$TARGET/download/$name" "${HTTP_ROOT}/download/${name}"
  done
}

mirror_ewas_db() {
  local root_url="${HTTP_ROOT}/EWAS_db/"
  local dest="$TARGET/EWAS_db"
  mkdir -p "$dest"
  printf '=== All Data EWAS_db -> %s ===\n' "$dest"
  printf 'Listing studies from %s\n' "$root_url"
  mapfile -t studies < <(list_hrefs "$root_url" | grep -E '/$' | sed 's|/$||' | grep -Ev '^(index\.html)?$')
  printf 'Found %s study directories\n' "${#studies[@]}"
  local i=0
  for study in "${studies[@]}"; do
    i=$((i + 1))
    mkdir -p "$dest/$study"
    printf '[%s/%s] %s\n' "$i" "${#studies[@]}" "$study"
    mapfile -t files < <(list_hrefs "${root_url}${study}/" | grep -Ev '/$')
    for f in "${files[@]}"; do
      wget -c -q -O "$dest/$study/$f" "${root_url}${study}/${f}" || {
        printf 'WARN: failed %s/%s\n' "$study" "$f" >&2
      }
    done
  done
}

mirror_add_ewas_db() {
  local root_url="${HTTP_ROOT}/add_ewas_db/"
  local dest="$TARGET/add_ewas_db"
  mkdir -p "$dest"
  printf '=== add_ewas_db -> %s ===\n' "$dest"
  mapfile -t entries < <(list_hrefs "$root_url")
  for entry in "${entries[@]}"; do
    if [[ "$entry" == */ ]]; then
      local sub="${entry%/}"
      mkdir -p "$dest/$sub"
      mapfile -t files < <(list_hrefs "${root_url}${sub}/" | grep -Ev '/$')
      printf '  %s (%s files)\n' "$sub" "${#files[@]}"
      for f in "${files[@]}"; do
        wget -c -q -O "$dest/$sub/$f" "${root_url}${sub}/${f}" || true
      done
    else
      wget -c -q -O "$dest/$entry" "${root_url}${entry}" || true
    fi
  done
}

MODE="${1:-all}"
case "$MODE" in
  EWAS_db) mirror_ewas_db ;;
  add_ewas_db) mirror_add_ewas_db ;;
  download|baseline) mirror_download_zips ;;
  all)
    mirror_download_zips
    mirror_add_ewas_db
    mirror_ewas_db
    ;;
  *)
    printf 'Usage: %s [all|EWAS_db|add_ewas_db|download|baseline]\n' "$(basename "$0")" >&2
    exit 1
    ;;
esac

printf 'EWAS DataHub download finished: %s\n' "$TARGET"
