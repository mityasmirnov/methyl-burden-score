#!/usr/bin/env bash
# P5-style 30-epoch early-stop retrain for N-light-gene-max then N-light-gene-mean.
# Sequential on GPU 0; skip enet; refresh analysis.md after each arm.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

log() { echo "[$(date '+%H:%M:%S')] $*"; }
export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1

MAIN_CFG=configs/experiment/stage0_7g_gene_only_probe.yaml
RUN=scripts/run_7g_gene_only_probe.py
REPORT=scripts/write_7g_gene_only_probe_report.py
ART=artifacts/runs

summarize_arm() {
  local arm_id="$1"
  local json="reports/inspection/stage0_7g_gene_only_probe/per_arm/${arm_id}.json"
  if [[ ! -f "$json" ]]; then
    log "WARNING: missing $json"
    return
  fi
  ARM_ID="$arm_id" JSON_PATH="$json" python - <<'PY'
import json, os
from pathlib import Path
arm_id = os.environ["ARM_ID"]
d = json.loads(Path(os.environ["JSON_PATH"]).read_text())
folds = d.get("folds") or []
e2es, lins, ages, sexes = [], [], [], []
for f in folds:
    ev = f.get("evaluations") or {}
    e2e = (ev.get("mbs_e2e") or {}).get("metrics") or {}
    lin = (ev.get("mbs_linear_probe") or {}).get("metrics") or {}
    if "tissue_macro_f1" in e2e:
        e2es.append(float(e2e["tissue_macro_f1"]))
    if "tissue_macro_f1" in lin:
        lins.append(float(lin["tissue_macro_f1"]))
    if "age_mae" in e2e:
        ages.append(float(e2e["age_mae"]))
    sex = e2e.get("sex_auroc") or lin.get("sex_auroc")
    if sex is not None:
        sexes.append(float(sex))
    be = (f.get("checkpoint_selection") or {}).get("best_epoch")
    print(f"  fold best_epoch={be} e2e_f1={e2e.get('tissue_macro_f1')} linear_f1={lin.get('tissue_macro_f1')}")

def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")
print(
    f"[{arm_id}] n_folds={len(folds)} mean_e2e_f1={mean(e2es):.3f} "
    f"mean_linear_f1={mean(lins):.3f} mean_age_mae={mean(ages):.2f} "
    f"mean_sex_auroc={mean(sexes):.3f}"
)
PY
}

# Move 5-epoch runs aside so train does not skip-if-done
for prefix in stage0-7g-gene-probe-light-max stage0-7g-gene-probe-light-mean; do
  for f in 0 1 2; do
    src="$ART/${prefix}-f${f}"
    dst="$ART/${prefix}-f${f}.stale-5ep"
    if [[ -d "$src" && ! -e "$dst" ]]; then
      mv "$src" "$dst"
      log "moved $src → $dst"
    elif [[ -d "$src" ]]; then
      rm -rf "$src"
      log "removed $src (stale-5ep already present)"
    fi
  done
done

log "=== 1/2 N-light-gene-max (30 ep, patience 5) ==="
# Omit --fold so all 3 folds train (--fold is single-value, last wins).
uv run python "$RUN" --config "$MAIN_CFG" --device cuda \
  --arm N-light-gene-max
log "N-light-gene-max done — refreshing analysis.md"
uv run python "$REPORT"
summarize_arm N-light-gene-max

log "=== 2/2 N-light-gene-mean (30 ep, patience 5) ==="
uv run python "$RUN" --config "$MAIN_CFG" --device cuda \
  --arm N-light-gene-mean
log "N-light-gene-mean done — refreshing analysis.md"
uv run python "$REPORT"
summarize_arm N-light-gene-mean

log "=== All done ==="
summarize_arm N-light-gene-max
summarize_arm N-light-gene-mean
