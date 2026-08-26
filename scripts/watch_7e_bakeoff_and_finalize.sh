#!/usr/bin/env bash
# Wait for parallel 7E bake-off slices, merge report, mark Milestone 7E done.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source scripts/activate_data_environment.sh >/dev/null

FLAT="$ROOT/reports/inspection/stage0_7e_dev_cv_flat/summary.json"
HIER="$ROOT/reports/inspection/stage0_7e_dev_cv_hier/summary.json"
MP="$ROOT/reports/inspection/stage0_7e_dev_cv_mp/summary.json"
OUT="$ROOT/reports/inspection/stage0_7e_dev_cv"
LOG="$ROOT/scratch/stage0_7e_watch.log"

log() { echo "[$(date -Is)] $*" | tee -a "$LOG"; }

log "watching for flat/hier/mp summaries…"
while true; do
  flat_n=0; hier_n=0; mp_n=0
  for d in artifacts/runs/stage0-7e-*/; do
    [[ -f "$d/metrics.json" ]] || continue
    b=$(basename "$d")
    case "$b" in
      *flat*|*T-mean*|*T-enet*) flat_n=$((flat_n+1)) ;;
      *hier*) hier_n=$((hier_n+1)) ;;
      *) mp_n=$((mp_n+1)) ;;
    esac
  done
  log "metrics: flatish=$flat_n hier=$hier_n mpish=$mp_n | slices: flat=$([[ -f $FLAT ]] && echo y || echo n) hier=$([[ -f $HIER ]] && echo y || echo n) mp=$([[ -f $MP ]] && echo y || echo n)"
  if [[ -f "$FLAT" && -f "$HIER" && -f "$MP" ]]; then
    break
  fi
  # keep jobs alive note
  if ! pgrep -f "mbs train dev-cv" >/dev/null; then
    log "WARNING: no mbs train dev-cv processes; waiting for existing summaries only"
  fi
  sleep 300
done

log "all slices present; merging"
uv run python scripts/merge_7e_dev_cv_reports.py | tee -a "$LOG"

if [[ ! -f "$OUT/summary.md" ]]; then
  log "merge failed: missing $OUT/summary.md"
  exit 1
fi

WINNER=$(python3 -c "import json; print(json.load(open('$OUT/summary.json'))['winner'].get('arm'))")
log "winner=$WINNER"

# Update pipeline + plan docs only when report names a winner.
python3 - <<'PY' | tee -a "$ROOT/scratch/stage0_7e_watch.log"
from pathlib import Path
import json
import re

root = Path("/data/projects/methyl-burden-score")
summary = json.loads((root / "reports/inspection/stage0_7e_dev_cv/summary.json").read_text())
winner = summary.get("winner") or {}
arm = winner.get("arm")
if not arm:
    raise SystemExit("no winner in summary")

todo = root / "docs/TODO_PIPELINE.md"
text = todo.read_text(encoding="utf-8")
old = "- **Status:** `pending` (**current coding gate**)"
# Only the 7E section uses this exact phrase for pending current gate.
if old not in text:
    raise SystemExit("could not find 7E status line")
# Replace first occurrence after '## 7E. Development'
idx = text.find("## 7E. Development cross-validation")
if idx < 0:
    raise SystemExit("7E heading missing")
idx2 = text.find(old, idx)
if idx2 < 0 or idx2 > idx + 800:
    raise SystemExit("7E status line not near heading")
text = (
    text[:idx2]
    + f"- **Status:** `done` (winner `{arm}`; report `reports/inspection/stage0_7e_dev_cv/`)"
    + text[idx2 + len(old) :]
)
# Next action under 7E
text = text.replace(
    "- **Next action:** Start 7E ATS bake-off; run **7E′** Hub multitask + hygiene\n"
    "  in parallel (required before Milestone 7).",
    "- **Next action:** Run **7E′** Hub multitask + hygiene (required before Milestone 7);\n"
    f"  architecture for Milestone 7 is `{arm}`.",
    1,
)
todo.write_text(text, encoding="utf-8")
print("updated TODO_PIPELINE.md")

plan = root / "docs/plans/milestone-7e-development-cv.md"
p = plan.read_text(encoding="utf-8")
p = p.replace(
    "Status: **pending** (current coding gate).",
    f"Status: **done**. Winner: `{arm}`. Report: `reports/inspection/stage0_7e_dev_cv/`.",
    1,
)
if "## Evidence" not in p:
    p += (
        "\n## Evidence\n\n"
        "- Report: `reports/inspection/stage0_7e_dev_cv/{summary.md,summary.json}`\n"
        f"- Winner: `{arm}` "
        f"(tissue macro-F1={winner.get('tissue_macro_f1')}, "
        f"age MAE={winner.get('age_mae')})\n"
        "- CV budget: max_loci=8192, max_epochs=2, 3 outer folds × 2 restarts\n"
        "- Protocol: shared `hub-ats-7e-3fold-v1` folds; Level-1 A/B; CpGPT-off; "
        "late-fusion multipath/gene-direct on graph-v2\n"
    )
plan.write_text(p, encoding="utf-8")
print("updated milestone-7e-development-cv.md")

post = root / "docs/plans/post-v0-scientific-programme.md"
if post.is_file():
    t = post.read_text(encoding="utf-8")
    t2 = t
    t2 = t2.replace("`graph-v2` not on disk", "`graph-v2` on disk (closed in 7E prep)")
    t2 = t2.replace("branch masks not wired", "branch masks train-time (closed in 7E prep)")
    if t2 != t:
        post.write_text(t2, encoding="utf-8")
        print("updated post-v0-scientific-programme.md")

rule = root / ".cursor/rules/pipeline-todo.mdc"
if rule.is_file():
    r = rule.read_text(encoding="utf-8")
    r2 = r.replace("(**current gate**; graph-v2 on disk; multi-path unblocked)", "(**done**)")
    r2 = r2.replace("coding gate is **7E**", "coding gate is **7E′** (7E done)")
    if r2 != r:
        rule.write_text(r2, encoding="utf-8")
        print("updated pipeline-todo.mdc")
print("finalize ok")
PY

log "7E watch complete"
touch "$ROOT/scratch/stage0_7e_watch.done"
