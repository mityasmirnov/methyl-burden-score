#!/usr/bin/env bash
# Resume N-light 30-ep retrain after OOM: keep completed folds, redo incomplete.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

log() { echo "[$(date '+%H:%M:%S')] $*"; }
export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

MAIN_CFG=configs/experiment/stage0_7g_gene_only_probe.yaml
RUN=scripts/run_7g_gene_only_probe.py
REPORT=scripts/write_7g_gene_only_probe_report.py
ART=artifacts/runs

sync_per_arm_and_report() {
  local arm_id="$1"
  local prefix="$2"
  python3 - <<PY
import json
from pathlib import Path
art = Path("$ART")
report = Path("reports/inspection/stage0_7g_gene_only_probe/per_arm")
report.mkdir(parents=True, exist_ok=True)
folds = []
for i in range(3):
    p = art / f"$prefix-f{i}" / "metrics.json"
    if p.is_file():
        folds.append(json.loads(p.read_text()))
out = {"arm_id": "$arm_id", "kind": "flat_region_train", "folds": folds, "n_folds": len(folds)}
(report / f"$arm_id.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
print(f"synced $arm_id n_folds={len(folds)}")
PY
  uv run python "$REPORT"
  ARM_ID="$arm_id" python3 - <<'PY'
import json, os
from pathlib import Path
arm = os.environ["ARM_ID"]
d = json.loads(Path(f"reports/inspection/stage0_7g_gene_only_probe/per_arm/{arm}.json").read_text())
folds = d.get("folds") or []
e2es, bests, rans = [], [], []
for i, f in enumerate(folds):
    ev = ((f.get("evaluations") or {}).get("mbs_e2e") or {}).get("metrics") or {}
    tf = (ev.get("tissue") or {}).get("macro_f1")
    be = (f.get("checkpoint_selection") or {}).get("best_epoch", f.get("best_epoch"))
    hist = f.get("history") or []
    ran = hist[-1].get("epoch") if hist and isinstance(hist[-1], dict) else len(hist)
    print(f"  fold{i}: e2e_f1={tf} best_epoch={be} epochs_ran={ran}")
    if tf is not None:
        e2es.append(float(tf))
    if be is not None:
        bests.append(int(be))
    if ran:
        rans.append(int(ran))
if e2es:
    print(f"[{arm}] mean_e2e_f1={sum(e2es)/len(e2es):.3f} n={len(e2es)} best_eps={bests} ran={rans}")
PY
}

# Drop incomplete f1 (metrics.jsonl only, no metrics.json)
for d in \
  "$ART/stage0-7g-gene-probe-light-max-f1" \
  "$ART/stage0-7g-gene-probe-light-max-f2" \
  "$ART/stage0-7g-gene-probe-light-mean-f0" \
  "$ART/stage0-7g-gene-probe-light-mean-f1" \
  "$ART/stage0-7g-gene-probe-light-mean-f2"
do
  if [[ -d "$d" && ! -f "$d/metrics.json" ]]; then
    rm -rf "$d"
    log "removed incomplete $d"
  fi
done

log "=== refresh analysis.md from completed folds ==="
sync_per_arm_and_report N-light-gene-max stage0-7g-gene-probe-light-max

log "=== resume N-light-gene-max (skip-if-done keeps f0) ==="
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm N-light-gene-max
log "N-light-gene-max done"
sync_per_arm_and_report N-light-gene-max stage0-7g-gene-probe-light-max

log "=== N-light-gene-mean all folds ==="
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm N-light-gene-mean
log "N-light-gene-mean done"
sync_per_arm_and_report N-light-gene-mean stage0-7g-gene-probe-light-mean

log "=== All done ==="
