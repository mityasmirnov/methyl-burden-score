#!/usr/bin/env bash
# Milestone 10 GPU-0 keeper queue.
#
# Does NOT interrupt live nine-pack full refs. After they exit:
#   1) repair m-only fold-0 if smoke poisoned skip-if-done
#   2) refresh nine-pack P2-G / m-only report
#   3) nine-pack vector RBS @ 15 ep (mean→max, max→max) — scale test vs scalar P2-G
#   4) Track B.4 ATS seed-43 2×2 pooling
#   5) ATS one-hop light-mean seed-43
# Soft-stop before Milestone 11 Stage B / 12 OOF / disease GPU.
#
# Architecture policy: P2-G is provisional only. Do not declare a cascade
# primary until nine-pack scalar-vs-vector (and N-light) are compared.
# Milestone 12 OOF is for finalists only (winning cascade + N-light), not all arms.
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
  while pgrep -af 'run_7h_nine_pack_smoke.py --phase full|stage0-7h-nine-pack-P2-G[^-]|stage0-7h-nine-pack-m-only-f|stage0-7h-nine-pack-vector' \
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

train_ninepack_vector_rbs() {
  log "=== GPU-0: nine-pack vector RBS scale screen (vs provisional scalar P2-G) ==="
  log "P2-G is NOT locked — vector may win at 34k samples."
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
  log "=== write vector-vs-scalar comparison snippet ==="
  uv run python - <<'PY'
import json, statistics
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(".")
ART = ROOT / "artifacts" / "runs"
arms = {
    "P2-G scalar max/max": "stage0-7h-nine-pack-P2-G",
    "vector mean→max": "stage0-7h-nine-pack-vector-mean-max",
    "vector max→max": "stage0-7h-nine-pack-vector-max-max",
    "N-light m-only": None,  # per-fold run ids
}

def fold_metrics(run_root: Path):
    rows = []
    for i in range(3):
        mp = run_root / f"fold_{i}" / "metrics.json"
        if not mp.is_file():
            continue
        b = json.loads(mp.read_text())
        e = (b.get("evaluations") or {}).get("mbs_e2e") or {}
        m = e.get("metrics") or {}
        rows.append({
            "fold": i,
            "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
            "age_mae": (m.get("age") or {}).get("mae"),
            "sex_auroc": (m.get("sex") or {}).get("auroc"),
            "best_epoch": b.get("best_epoch") or (b.get("checkpoint_selection") or {}).get("best_epoch"),
        })
    return rows

def summarize(rows):
    def mean(k):
        vals = [float(r[k]) for r in rows if isinstance(r.get(k), (int, float))]
        return float(statistics.mean(vals)) if vals else None
    return {
        "n_folds": len(rows),
        "folds": rows,
        "tissue_f1_mean": mean("tissue_f1"),
        "age_mae_mean": mean("age_mae"),
        "sex_auroc_mean": mean("sex_auroc"),
    }

out = {
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "note": (
        "Nine-pack scalar P2-G vs vector RBS at matched 15-ep / 65k loci. "
        "P2-G remains provisional until this comparison is complete. "
        "N-light m-only is a co-equal light-model finalist for Milestone 12 OOF."
    ),
    "arms": {},
}
for label, rid in arms.items():
    if rid is None:
        rows = []
        for i in range(3):
            mp = ART / f"stage0-7h-nine-pack-m-only-f{i}" / "metrics.json"
            if not mp.is_file():
                continue
            b = json.loads(mp.read_text())
            # skip smoke-poisoned short history
            if len(b.get("history") or []) <= 3 and i == 0:
                continue
            e = (b.get("evaluations") or {}).get("mbs_e2e") or {}
            m = e.get("metrics") or {}
            rows.append({
                "fold": i,
                "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
                "age_mae": (m.get("age") or {}).get("mae"),
                "sex_auroc": (m.get("sex") or {}).get("auroc"),
                "best_epoch": b.get("best_epoch"),
            })
        out["arms"][label] = summarize(rows)
    else:
        out["arms"][label] = summarize(fold_metrics(ART / rid))

rep = ROOT / "reports" / "inspection" / "stage0_7h_nine_pack_smoke"
rep.mkdir(parents=True, exist_ok=True)
(rep / "vector_vs_scalar.json").write_text(json.dumps(out, indent=2) + "\n")
lines = [
    "# Nine-pack scalar vs vector RBS (Milestone 10 scale screen)",
    "",
    f"Generated: `{out['generated_at']}`",
    "",
    out["note"],
    "",
    "| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |",
    "|-----|------:|----------:|--------:|----------:|",
]
for label, arm in out["arms"].items():
    t, a, s = arm.get("tissue_f1_mean"), arm.get("age_mae_mean"), arm.get("sex_auroc_mean")
    def fmt(x):
        return "—" if x is None else f"{x:.3f}"
    lines.append(f"| {label} | {arm['n_folds']} | {fmt(t)} | {fmt(a)} | {fmt(s)} |")
lines += [
    "",
    "**Do not** call a cascade primary from ATS alone. Narrow Milestone **12** OOF",
    "to finalists: winning cascade (scalar or vector) + **N-light** light model.",
    "",
]
(rep / "vector_vs_scalar.md").write_text("\n".join(lines) + "\n")
print("wrote", rep / "vector_vs_scalar.md")
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

log "=== refresh nine-pack P2-G / m-only report ==="
uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || true

train_ninepack_vector_rbs

log "=== Track B.4: ATS seed-43 2×2 pooling ==="
bash scripts/run_7h_ats_pooling_s2.sh

train_ats_light_mean_s2

log "=== SOFT STOP (GPU 0 may go idle) ==="
log "Do NOT auto-launch Milestone 11 Stage B or Milestone 12 OOF."
log "Review vector_vs_scalar.md — pick cascade finalist; keep N-light as co-equal light model."
log "Milestone 12 OOF = finalists only (not all ~9 arms)."
log "Reports:"
log "  $REPORT/analysis.md"
log "  $REPORT/vector_vs_scalar.md"
log "  reports/inspection/stage0_7h_ats_pooling_s2/analysis.md"
log "=== queue done ==="
