#!/usr/bin/env bash
# Milestone 10 GPU-0 keeper queue (priority order).
#
# After live nine-pack full / m-only-f0 repair:
#   1) refresh P2-G / m-only report
#   2) HIGH: nine-pack vector RBS @ 15 ep (mean→max, max→max)
#   3) HIGH: one-hop correctness smokes (seed-mask G0/G1 + multi-seed) — cheap
#   4) Track B.4 ATS seed-43 2×2 pooling (valuable; already planned)
# Soft-stop before Stage B / Milestone 12 OOF / disease GPU.
# Demoted: ATS light-mean s2 keep-busy (lower value than 2–3).
#
# Architecture policy: P2-G provisional only. Milestone 12 = finalists only
# (winning cascade + N-light). Trait expansion = age/tissue/sex until 10c
# label_status write; blood/brain/bmi/ancestry deferred.
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
REPORT=reports/inspection/stage0_7h_nine_pack_smoke

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
  [[ -z "$pid" ]] && { log "empty $pidf — assuming $label finished"; return 0; }
  if ! kill -0 "$pid" 2>/dev/null; then
    log "$label pid $pid not alive — continuing"
    return 0
  fi
  log "waiting for $label pid=$pid (poll ${POLL_SEC}s)"
  while kill -0 "$pid" 2>/dev/null; do sleep "$POLL_SEC"; done
  log "$label pid=$pid exited"
}

wait_for_ninepack_gpu_owner() {
  while pgrep -af 'run_7h_nine_pack_smoke.py --phase full|stage0-7h-nine-pack-P2-G[^-]|stage0-7h-nine-pack-m-only-f|stage0-7h-nine-pack-vector' \
      | rg -v 'rg |pgrep|next_queue|m-only-smoke|smoke-3ep|handoff' >/dev/null; do
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
raise SystemExit(0 if len(b.get("history") or []) <= 3 else 1)
PY
  then
    local dest="${run}-smoke-3ep-legacy"
    log "smoke-poisoned m-only-f0 → $dest"
    rm -rf "$dest"
    mv "$run" "$dest"
    need_retrain=1
  else
    log "m-only-f0 already full-budget — skip repair"
  fi
  [[ "$need_retrain" -eq 1 ]] || return 0
  log "=== GPU-0: retrain m-only fold-0 @ full budget ==="
  uv run python - <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve() / 'scripts'))
from run_7h_nine_pack_smoke import flat_train
flat_train(fold_filter=0, max_epochs=None, run_prefix='stage0-7h-nine-pack-m-only')
print('[queue] m-only-f0 retrain done', flush=True)
PY
}

train_ninepack_vector_rbs() {
  log "=== HIGH: nine-pack vector RBS scale screen vs provisional scalar P2-G ==="
  local arms=(
    "configs/experiment/stage0_7h_nine_pack_vector_mean_max.yaml|stage0-7h-nine-pack-vector-mean-max|_staging_vector_mean_max"
    "configs/experiment/stage0_7h_nine_pack_vector_max_max.yaml|stage0-7h-nine-pack-vector-max-max|_staging_vector_max_max"
  )
  local entry cfg run_id staging
  for entry in "${arms[@]}"; do
    IFS='|' read -r cfg run_id staging <<<"$entry"
    log "=== train $run_id ==="
    uv run mbs train cascade \
      --config "$cfg" \
      --run-id "$run_id" \
      --device cuda \
      --skip-if-done \
      --report-dir "$REPORT/$staging"
  done
  uv run python - <<'PY'
import json, statistics
from datetime import datetime, timezone
from pathlib import Path
ART = Path('artifacts/runs')
arms = {
    'P2-G scalar max/max': 'stage0-7h-nine-pack-P2-G',
    'vector mean→max': 'stage0-7h-nine-pack-vector-mean-max',
    'vector max→max': 'stage0-7h-nine-pack-vector-max-max',
}
def summarize(run_id=None, m_only=False):
    rows=[]
    if m_only:
        for i in range(3):
            mp = ART / f'stage0-7h-nine-pack-m-only-f{i}' / 'metrics.json'
            if not mp.is_file():
                continue
            b=json.loads(mp.read_text())
            if i==0 and len(b.get('history') or [])<=3:
                continue
            e=(b.get('evaluations')or{}).get('mbs_e2e')or{}
            m=e.get('metrics')or{}
            rows.append({'tissue_f1':(m.get('tissue')or{}).get('macro_f1'),'age_mae':(m.get('age')or{}).get('mae'),'sex_auroc':(m.get('sex')or{}).get('auroc')})
    else:
        for i in range(3):
            mp = ART / run_id / f'fold_{i}' / 'metrics.json'
            if not mp.is_file():
                continue
            b=json.loads(mp.read_text())
            e=(b.get('evaluations')or{}).get('mbs_e2e')or{}
            m=e.get('metrics')or{}
            rows.append({'tissue_f1':(m.get('tissue')or{}).get('macro_f1'),'age_mae':(m.get('age')or{}).get('mae'),'sex_auroc':(m.get('sex')or{}).get('auroc')})
    def mean(k):
        vals=[float(r[k]) for r in rows if isinstance(r.get(k),(int,float))]
        return float(statistics.mean(vals)) if vals else None
    return {'n_folds':len(rows),'tissue_f1_mean':mean('tissue_f1'),'age_mae_mean':mean('age_mae'),'sex_auroc_mean':mean('sex_auroc'),'folds':rows}
out={'generated_at':datetime.now(timezone.utc).isoformat(),'note':'P2-G provisional. Pick cascade finalist only after this screen. N-light is co-equal OOF finalist.','arms':{}}
for label, rid in arms.items():
    out['arms'][label]=summarize(rid)
out['arms']['N-light m-only']=summarize(m_only=True)
rep=Path('reports/inspection/stage0_7h_nine_pack_smoke')
rep.mkdir(parents=True, exist_ok=True)
(rep/'vector_vs_scalar.json').write_text(json.dumps(out, indent=2)+'\n')
lines=['# Nine-pack scalar vs vector RBS','',f"Generated: `{out['generated_at']}`",'',out['note'],'','| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |','|-----|------:|----------:|--------:|----------:|']
for label, arm in out['arms'].items():
    def fmt(x):
        return '—' if x is None else f'{x:.3f}'
    lines.append(f"| {label} | {arm['n_folds']} | {fmt(arm['tissue_f1_mean'])} | {fmt(arm['age_mae_mean'])} | {fmt(arm['sex_auroc_mean'])} |")
lines += ['','Milestone **12** OOF = finalists only (winning cascade + N-light).','']
(rep/'vector_vs_scalar.md').write_text('\n'.join(lines)+'\n')
print('wrote', rep/'vector_vs_scalar.md')
PY
}

run_onehop_correctness() {
  log "=== HIGH: one-hop correctness smokes (seed-mask + multi-seed, fold0, 5 ep) ==="
  uv run python -u scripts/run_7h_onehop_correctness_smokes.py --device cuda --epochs 5
}

log "=== 7H GPU-0 keeper queue start (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES) ==="
wait_for_pidfile "$FULL_PID_FILE" "nine-pack-full"
wait_for_ninepack_gpu_owner
repair_m_only_fold0_if_smoke

log "=== refresh nine-pack P2-G / m-only report ==="
uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || true

train_ninepack_vector_rbs
run_onehop_correctness

log "=== Track B.4: ATS seed-43 2×2 pooling ==="
bash scripts/run_7h_ats_pooling_s2.sh

log "=== SOFT STOP ==="
log "Skipped low-value ATS light-mean s2 keep-busy."
log "Do NOT auto-launch Milestone 11 Stage B / 12 OOF / disease GPU."
log "Next manual: review vector_vs_scalar.md + onehop correctness; 10c uses label_status table."
log "Usable traits today: age/tissue/sex only. Disease/cancer need label_status write; blood/brain/bmi/ancestry deferred."
log "=== queue done ==="
