#!/usr/bin/env bash
# Milestone 10 GPU-0 keeper queue.
#
# Keeps CUDA device 0 busy with chained Milestone-10 work. Does NOT interrupt
# the live nine-pack full refs. After they exit:
#   1) repair m-only fold-0 if smoke poisoned skip-if-done
#   2) refresh nine-pack report
#   3) Track B.4 ATS seed-43 2×2 pooling (4 arms)
#   4) ATS one-hop light-mean seed-43 (extra keep-busy)
# Soft-stop before Milestone 11 Stage B / 13 OOF / disease GPU.
#
# Policy: prefer continuous GPU-0 occupancy; poll frequently on handoff.
#
# Usage:
#   nohup bash scripts/run_7h_next_queue.sh >> scratch/logs/7h_next_queue.log 2>&1 &
#   echo $! > scratch/logs/7h_next_queue.pid
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

LOG_DIR=scratch/logs
mkdir -p "$LOG_DIR"
FULL_PID_FILE="$LOG_DIR/7h_nine_pack_full.pid"
POLL_SEC="${MBS_GPU0_POLL_SEC:-30}"

log() { echo "[$(date -Is)] $*"; }

wait_for_pidfile() {
  local pidf="$1"
  local label="$2"
  if [[ ! -f "$pidf" ]]; then
    log "no $pidf — assuming $label already finished"
    return 0
  fi
  local pid
  pid="$(tr -d '[:space:]' <"$pidf" || true)"
  if [[ -z "$pid" ]]; then
    log "empty $pidf — assuming $label finished"
    return 0
  fi
  if ! kill -0 "$pid" 2>/dev/null; then
    log "$label pid $pid not alive — continuing"
    return 0
  fi
  log "waiting for $label pid=$pid (poll ${POLL_SEC}s)"
  while kill -0 "$pid" 2>/dev/null; do
    sleep "$POLL_SEC"
  done
  log "$label pid=$pid exited"
}

wait_for_ninepack_gpu_owner() {
  while pgrep -af 'run_7h_nine_pack_smoke.py --phase full|stage0-7h-nine-pack-P2-G[^-]|stage0-7h-nine-pack-m-only-f' \
      | rg -v 'rg |pgrep|next_queue|m-only-smoke|m-only-f0-smoke|smoke-3ep' >/dev/null; do
    log "nine-pack GPU owner still present — sleep ${POLL_SEC}s"
    sleep "$POLL_SEC"
  done
}

repair_m_only_fold0_if_smoke() {
  local run=artifacts/runs/stage0-7h-nine-pack-m-only-f0
  local metrics="$run/metrics.json"
  local need_retrain=0
  if [[ ! -f "$metrics" ]]; then
    log "m-only-f0 missing — will train full budget"
    need_retrain=1
  elif python3 - <<'PY'
import json
from pathlib import Path
b = json.loads(Path("artifacts/runs/stage0-7h-nine-pack-m-only-f0/metrics.json").read_text())
hist = b.get("history") or []
# Smoke used max_epochs=3 → history length 3. Full budget is 16.
raise SystemExit(0 if len(hist) <= 3 else 1)
PY
  then
    local dest="${run}-smoke-3ep-legacy"
    log "smoke-poisoned m-only-f0 (history<=3) → move to $dest"
    rm -rf "$dest"
    mv "$run" "$dest"
    need_retrain=1
  else
    log "m-only-f0 already full-budget — skip repair"
  fi
  if [[ "$need_retrain" -ne 1 ]]; then
    return 0
  fi
  log "=== GPU-0: retrain m-only fold-0 @ full budget ==="
  uv run python - <<'PY'
import sys
from pathlib import Path
ROOT = Path(".").resolve()
sys.path.insert(0, str(ROOT / "scripts"))
from run_7h_nine_pack_smoke import flat_train

flat_train(fold_filter=0, max_epochs=None, run_prefix="stage0-7h-nine-pack-m-only")
print("[queue] m-only-f0 retrain done", flush=True)
PY
}

train_ats_light_mean_s2() {
  log "=== GPU-0: ATS one-hop light-mean seed-43 ==="
  mkdir -p configs/experiment/_generated
  local cfg=configs/experiment/_generated/stage0_7h_ats_light_mean_s2.yaml
  local src=configs/experiment/stage0_7g_gene_only_probe_light_mean.yaml
  if [[ ! -f "$src" ]]; then
    log "missing $src — skip light-mean s2"
    return 0
  fi
  python3 - <<PY
from pathlib import Path
import yaml
src = Path("$src")
cfg = yaml.safe_load(src.read_text())
cfg.setdefault("experiment", {})["seed"] = 43
cfg["experiment"]["name"] = str(cfg["experiment"].get("name", "light_mean")) + "_s2"
cfg.setdefault("pilot", {}).setdefault("matrix_id", "matrix-hub-age-tissue-sex-full-v1")
cfg.setdefault("split_id", "hub-ats-7e-3fold-v1")
Path("$cfg").write_text(yaml.safe_dump(cfg, sort_keys=False))
print("wrote", "$cfg")
PY
  uv run python - <<'PY'
import sys
from pathlib import Path
from mbs.paths import DataPaths

sys.path.insert(0, "scripts")
from run_7g_gene_only_probe import train_flat_region_arm

paths = DataPaths.from_environment()
cfg = Path("configs/experiment/_generated/stage0_7h_ats_light_mean_s2.yaml")
train_flat_region_arm(
    paths=paths,
    config_path=cfg,
    run_prefix="stage0-7h-ats-light-mean-s2",
    device="cuda",
    report_dir=Path("reports/inspection/stage0_7h_ats_pooling_s2"),
    fold_filter=None,
)
print("[queue] light-mean s2 done", flush=True)
PY
}

log "=== 7H GPU-0 keeper queue start (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES) ==="
wait_for_pidfile "$FULL_PID_FILE" "nine-pack-full"
wait_for_ninepack_gpu_owner

repair_m_only_fold0_if_smoke

log "=== refresh nine-pack full report ==="
uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || true

log "=== Track B.4: ATS seed-43 2×2 pooling (occupy GPU 0) ==="
bash scripts/run_7h_ats_pooling_s2.sh

train_ats_light_mean_s2

log "=== SOFT STOP (GPU 0 may go idle) ==="
log "Auto-queue will not launch Milestone 11 Stage B / 13 OOF / disease GPU."
log "To keep GPU 0 busy next: review 10a/10b reports, then manually start Stage B"
log "  or 10c-gated trait arms. Prefer CUDA_VISIBLE_DEVICES=0."
log "Reports:"
log "  reports/inspection/stage0_7h_nine_pack_smoke/analysis.md"
log "  reports/inspection/stage0_7h_ats_pooling_s2/analysis.md"
log "=== queue done ==="
