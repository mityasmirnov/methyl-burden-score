#!/usr/bin/env bash
# run_nlight_retrain_all.sh
#
# Orchestrate the full N-light-gene retrain + annotation ablation grid on GPU 0.
#
# Order:
#   1. N-light-gene-max  folds 1 & 2  (both bugs fixed: orientation-v2 + ckpt-selection)
#   2. N-light-gene-mean folds 0, 1, 2
#   3. post-hoc mbs_enet on light-max and light-mean
#   4. Annotation ablation grid — fold 0, 18 arms (9 feature modes × 2 seeds)
#   5. Final report regeneration
#   6. git commit + push
#
# Usage (from repo root, after activate):
#   nohup bash scripts/run_nlight_retrain_all.sh \
#     > logs/nlight_retrain_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#
# GPU policy: CUDA_VISIBLE_DEVICES=0 (RTX 6000 Ada, ~48 GB free).
# Batch-token-budget 16M keeps peak VRAM ~6-8 GB per fold for 51k CpG panel.
# bf16-mixed is on in all training configs.

set -euo pipefail
cd "$(dirname "$0")/.."

# ── helpers ────────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

DEVICE=cuda
export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1

MAIN_CFG=configs/experiment/stage0_7g_gene_only_probe.yaml
RUN_SCRIPT=scripts/run_7g_gene_only_probe.py
ENET_SCRIPT=scripts/eval_mbs_enet_from_scores.py
REPORT_SCRIPT=scripts/write_7g_gene_only_probe_report.py

# ── sanity: both fixes must be on disk ────────────────────────────────────────
log "Sanity: confirming both code fixes..."
grep -q "stage_a_per_epoch_eval = bool.*or use_tissue_rank" src/mbs/training/loop.py \
    || die "checkpoint-selection fix missing in loop.py"
grep -q "mbs_for_heads = mbs_arr" src/mbs/training/flat_stage_a_eval.py \
    || die "orientation v2 fix missing in flat_stage_a_eval.py"
FIXES_COMMIT=$(git log -1 --format='%H %cI' -- src/mbs/training/loop.py src/mbs/training/flat_stage_a_eval.py)
log "Both fixes confirmed: $FIXES_COMMIT"

refresh_report() {
    log "Refreshing analysis.md..."
    uv run python "$REPORT_SCRIPT" 2>&1 || log "WARNING: report refresh failed (non-fatal)"
}

# ── 1. N-light-gene-max folds 1 & 2 ──────────────────────────────────────────
log "=== 1/6  N-light-gene-max folds 1 + 2 ==="
uv run python "$RUN_SCRIPT" \
    --config "$MAIN_CFG" --device "$DEVICE" \
    --arm N-light-gene-max --fold 1 --fold 2
log "N-light-gene-max f1+f2 done"
refresh_report

# ── 2. N-light-gene-mean folds 0, 1, 2 ───────────────────────────────────────
log "=== 2/6  N-light-gene-mean all 3 folds ==="
uv run python "$RUN_SCRIPT" \
    --config "$MAIN_CFG" --device "$DEVICE" \
    --arm N-light-gene-mean --fold 0 --fold 1 --fold 2
log "N-light-gene-mean all folds done"
refresh_report

# ── 3. Post-hoc mbs_enet ──────────────────────────────────────────────────────
log "=== 3/6  post-hoc mbs_enet on light-max and light-mean ==="
uv run python "$ENET_SCRIPT" \
    --run-prefix stage0-7g-gene-probe-light-max \
    --config "$MAIN_CFG" || log "WARNING: enet light-max failed (non-fatal)"
uv run python "$ENET_SCRIPT" \
    --run-prefix stage0-7g-gene-probe-light-mean \
    --config "$MAIN_CFG" || log "WARNING: enet light-mean failed (non-fatal)"
log "post-hoc enet done"
refresh_report

# ── 4. Annotation ablation grid (fold 0 only, 18 arms) ───────────────────────
log "=== 4/6  annotation ablation grid (18 arms, fold 0) ==="

ABLATION_ARM_IDS=(
    N-light-gene-ablation-m-only
    N-light-gene-ablation-m-only-s2
    N-light-gene-ablation-m-role
    N-light-gene-ablation-m-role-s2
    N-light-gene-ablation-m-context
    N-light-gene-ablation-m-context-s2
    N-light-gene-ablation-m-role-context
    N-light-gene-ablation-m-role-context-s2
    N-light-gene-ablation-full
    N-light-gene-ablation-full-s2
    N-light-gene-ablation-n0-obs-only
    N-light-gene-ablation-n0-obs-only-s2
    N-light-gene-ablation-n1-anno-only
    N-light-gene-ablation-n1-anno-only-s2
    N-light-gene-ablation-n2-reg-permuted
    N-light-gene-ablation-n2-reg-permuted-s2
    N-light-gene-ablation-n3-reg-zero
    N-light-gene-ablation-n3-reg-zero-s2
)

for ARM_ID in "${ABLATION_ARM_IDS[@]}"; do
    log "  ablation: $ARM_ID"
    uv run python "$RUN_SCRIPT" \
        --config "$MAIN_CFG" --device "$DEVICE" \
        --arm "$ARM_ID" --fold 0 \
    || log "  WARNING: arm $ARM_ID failed — continuing"
done

log "Ablation grid done"
refresh_report

# ── 5. Final report ───────────────────────────────────────────────────────────
log "=== 5/6  final report ==="
refresh_report

# ── 6. Commit & push ─────────────────────────────────────────────────────────
log "=== 6/6  commit and push ==="

# Stage only text/JSON/md artifacts — not large .npy/.zarr/.pt binaries
git add \
    "configs/experiment/stage0_7g_gene_only_probe.yaml" \
    "reports/inspection/stage0_7g_gene_only_probe/analysis.md" \
    "reports/inspection/stage0_7g_gene_only_probe/summary.json" \
    "reports/inspection/stage0_7g_gene_only_probe/task_comparison.json" \
    "reports/inspection/stage0_7g_gene_only_probe/gene_panel_manifest.json" \
    2>/dev/null || true

git add "reports/inspection/stage0_7g_gene_only_probe/per_arm/" 2>/dev/null || true

# Stage staging sub-dirs (arm-level staging JSONs, not model files)
find reports/inspection/stage0_7g_gene_only_probe/_staging_n_light_gene_* \
     -name "*.json" -o -name "*.md" 2>/dev/null | xargs git add 2>/dev/null || true

git diff --cached --quiet && {
    log "Nothing new to commit — done"
} || {
    git commit -m "feat(7g-prime): N-light-gene retrain (v2 contract), mbs_enet, ablation grid

- N-light-gene-max folds 1+2 retrained post fc8cd6f (orientation v2 +
  checkpoint-selection fix).
- N-light-gene-mean all 3 folds retrained clean.
- Post-hoc mbs_enet on both light arms.
- Annotation ablation grid: A0 m_only / A1 m_role / A2 m_context /
  A3 m_role_context / A4 full / N0 obs_only / N1 anno_only /
  N2 reg_permuted / N3 reg_zero — fold 0, two seeds each.
- analysis.md regenerated with all valid arm results.
- stage0_7g_gene_only_probe.yaml: added 18 ablation arm entries."
    git push && log "Pushed to origin/main"
}

log "=== All steps complete ==="
