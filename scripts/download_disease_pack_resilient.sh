#!/usr/bin/env bash
# Resilient resume for Hub disease_methylation_v1.zip.
# wget alone can exit after Connection refused / bogus 416 "fully retrieved".
# Usage: bash scripts/download_disease_pack_resilient.sh
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

URL="https://download.cncb.ac.cn/ewas/datahub/download/disease_methylation_v1.zip"
EXPECTED=21589344448
TARGET="$MBS_DATA_ROOT/raw/ewas_datahub/download"
OUT="$TARGET/disease_methylation_v1.zip"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET" "$LOGDIR"
LOG="$LOGDIR/ewas_disease_resilient_$(date -u +%Y%m%dT%H%M%SZ).log"

log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

eocd_ok() {
  python3 - "$1" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
with p.open("rb") as f:
    f.seek(max(0, p.stat().st_size - 65536))
    sys.exit(0 if f.read().rfind(b"PK\x05\x06") >= 0 else 1)
PY
}

log "start expected=$EXPECTED out=$OUT"
while true; do
  local_size="$(stat -c%s "$OUT" 2>/dev/null || echo 0)"
  if [[ "$local_size" -eq "$EXPECTED" ]]; then
    if eocd_ok "$OUT"; then
      log "COMPLETE size=$local_size EOCD_ok"
      exit 0
    fi
    log "size matches Content-Length but EOCD missing — removing and refetching"
    rm -f "$OUT"
    local_size=0
  fi

  log "wget -c resume local=$local_size / $EXPECTED"
  # --tries=0 infinite; retry connection refused; tolerate flaky CNCB
  wget -c --tries=0 --retry-connrefused --waitretry=30 --timeout=60 --read-timeout=120 \
    -O "$OUT" "$URL" >>"$LOG" 2>&1 || true
  local_size="$(stat -c%s "$OUT" 2>/dev/null || echo 0)"
  log "wget returned local=$local_size"

  if [[ "$local_size" -eq "$EXPECTED" ]]; then
    continue
  fi
  # Aug-6 failure mode: 416 + "fully retrieved" while still short
  if tail -n 40 "$LOG" | grep -q '416 Requested Range Not Satisfiable' \
    && [[ "$local_size" -gt 0 && "$local_size" -lt "$EXPECTED" ]]; then
    log "416 with short file — truncate 1 byte to force Range resume"
    truncate -s "$((local_size - 1))" "$OUT" || true
  fi
  sleep 30
done
