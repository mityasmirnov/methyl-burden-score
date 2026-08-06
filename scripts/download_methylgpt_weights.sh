#!/usr/bin/env bash
# Download MethylGPT pretrained checkpoints + probe vocabulary into
# $MBS_DATA_ROOT/raw/methylgpt (never under vendor/).
#
# Sources: vendor/methylgpt/README.md (Google Drive folders + Dropbox probe IDs).
# Default: methylgpt-medium (recommended). Pass --base and/or --large for others.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

PYTHON="${MBS_METHYLGPT_PYTHON:-$MBS_ROOT/.venv-methylgpt/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  PYTHON="$MBS_ROOT/.venv/bin/python"
fi

TARGET="$MBS_DATA_ROOT/raw/methylgpt"
MODELS_DIR="$TARGET/pretrained_models"
VOCAB_DIR="$TARGET/vocab"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
SCRATCH="$MBS_SCRATCH_ROOT/downloads/methylgpt"
mkdir -p "$MODELS_DIR" "$VOCAB_DIR" "$LOGDIR" "$SCRATCH"

WANT_BASE=0
WANT_MEDIUM=1
WANT_LARGE=0
for arg in "$@"; do
  case "$arg" in
    --base) WANT_BASE=1 ;;
    --medium) WANT_MEDIUM=1 ;;
    --large) WANT_LARGE=1 ;;
    --all) WANT_BASE=1; WANT_MEDIUM=1; WANT_LARGE=1 ;;
    --help|-h)
      printf 'Usage: %s [--base] [--medium] [--large] [--all]\n' "$0"
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$arg" >&2
      exit 1
      ;;
  esac
done

# Folder IDs from vendor/methylgpt/README.md
BASE_FOLDER_ID="1kWdmkkVQpU17uzUC6-wpNR_4UEdxGx6k"
MEDIUM_FOLDER_ID="14M4wdS83el9PAgh9TdfjSCeEcDPbz34f"
LARGE_FOLDER_ID="1lt8SF9MvoytPN3DeaxIss_ED9zNpf_Le"

PROBE_URL="https://www.dropbox.com/scl/fi/2n6bx7j8v0aon0kwfsghp/probe_ids_type3.csv?rlkey=ly133xlce1xxjiku6tiski6qq&st=pig4e41h&dl=1"

download_folder() {
  local folder_id="$1"
  local name="$2"
  local out="$MODELS_DIR/$name"
  mkdir -p "$out"
  if compgen -G "$out/*.pt" > /dev/null; then
    printf 'Skip %s (checkpoint already present in %s)\n' "$name" "$out"
    return 0
  fi
  printf 'Downloading %s (Google Drive folder %s) -> %s\n' "$name" "$folder_id" "$out"
  local tmp="$SCRATCH/$name"
  rm -rf "$tmp"
  mkdir -p "$tmp"
  "$PYTHON" - <<PY
import gdown
from pathlib import Path
url = "https://drive.google.com/drive/folders/${folder_id}"
out = Path(r"${tmp}")
gdown.download_folder(url=url, output=str(out), quiet=False, use_cookies=False)
PY
  # gdown may nest another directory; flatten one level if needed
  if [[ -z "$(find "$tmp" -maxdepth 1 -type f | head -n 1)" ]]; then
    local nested
    nested="$(find "$tmp" -mindepth 1 -maxdepth 1 -type d | head -n 1 || true)"
    if [[ -n "$nested" ]]; then
      mv "$nested"/* "$out"/
    fi
  else
    mv "$tmp"/* "$out"/
  fi
  rm -rf "$tmp"
  if ! compgen -G "$out/*.pt" > /dev/null; then
    printf 'ERROR: no .pt checkpoint found after downloading %s\n' "$name" >&2
    return 1
  fi
}

if [[ "$WANT_BASE" -eq 1 ]]; then
  download_folder "$BASE_FOLDER_ID" "methylgpt-base"
fi
if [[ "$WANT_MEDIUM" -eq 1 ]]; then
  download_folder "$MEDIUM_FOLDER_ID" "methylgpt-medium"
fi
if [[ "$WANT_LARGE" -eq 1 ]]; then
  download_folder "$LARGE_FOLDER_ID" "methylgpt-large"
fi

PROBE_PATH="$VOCAB_DIR/probe_ids_type3.csv"
if [[ -f "$PROBE_PATH" ]]; then
  printf 'Skip probe IDs (already at %s)\n' "$PROBE_PATH"
else
  printf 'Downloading probe_ids_type3.csv -> %s\n' "$PROBE_PATH"
  wget -c -O "$PROBE_PATH" "$PROBE_URL"
fi

cat > "$TARGET/SOURCE.txt" <<EOF
source: MethylGPT pretrained models + type3 probe vocabulary
repository: https://github.com/albert-ying/MethylGPT
vendor_path: vendor/methylgpt
medium_drive: https://drive.google.com/drive/folders/${MEDIUM_FOLDER_ID}
base_drive: https://drive.google.com/drive/folders/${BASE_FOLDER_ID}
large_drive: https://drive.google.com/drive/folders/${LARGE_FOLDER_ID}
probe_ids: Dropbox probe_ids_type3.csv (type3 default)
layout:
  pretrained_models/methylgpt-{base,medium,large}/
  vocab/probe_ids_type3.csv
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'MethylGPT weight download finished: %s\n' "$TARGET"
ls -lah "$MODELS_DIR"/*/ 2>/dev/null || true
ls -lah "$VOCAB_DIR"
