#!/usr/bin/env bash
# Track B.4: genuine seed-43 2×2 cascade pooling confirmation on ATS.
# Uses experiment.seed=43 so the seed-offset fix is exercised.
# Sequential on GPU 0; skip-if-done safe to resume.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

LOG=scratch/logs/7h_ats_pooling_s2.log
mkdir -p scratch/logs configs/experiment/_generated
exec >>"$LOG" 2>&1

log() { echo "[$(date -Is)] $*"; }

log "=== generating seed-43 configs ==="
python3 - <<'PY'
from pathlib import Path
import yaml

arms = {
    "p2": ("configs/experiment/stage0_7g_gene_only_probe_p2.yaml", "stage0-7h-ats-P2-G-s2"),
    "mean_max": (
        "configs/experiment/stage0_7g_gene_only_probe_scalar_mean_max.yaml",
        "stage0-7h-ats-scalar-mean-max-s2",
    ),
    "max_mean": (
        "configs/experiment/stage0_7g_gene_only_probe_scalar_max_mean.yaml",
        "stage0-7h-ats-scalar-max-mean-s2",
    ),
    "vector_mean_max": (
        "configs/experiment/stage0_7g_gene_only_probe_vector_mean_max.yaml",
        "stage0-7h-ats-vector-mean-max-s2",
    ),
}
out = Path("configs/experiment/_generated")
out.mkdir(parents=True, exist_ok=True)
manifest = []
for name, (src, run_id) in arms.items():
    p = Path(src)
    if not p.is_file():
        print(f"skip missing {src}")
        continue
    cfg = yaml.safe_load(p.read_text())
    cfg.setdefault("experiment", {})["seed"] = 43
    cfg["experiment"]["name"] = str(cfg["experiment"].get("name", name)) + "_s2"
    # Keep ATS matrix/split — do not accidentally point at nine-pack.
    pilot = cfg.setdefault("pilot", {})
    pilot.setdefault("matrix_id", "matrix-hub-age-tissue-sex-full-v1")
    cfg.setdefault("split_id", "hub-ats-7e-3fold-v1")
    dest = out / f"stage0_7h_ats_{name}_s2.yaml"
    dest.write_text(yaml.safe_dump(cfg, sort_keys=False))
    manifest.append((str(dest), run_id))
    print("wrote", dest, "run_id", run_id)
Path("configs/experiment/_generated/ats_pooling_s2_manifest.txt").write_text(
    "\n".join(f"{c}\t{r}" for c, r in manifest) + "\n"
)
PY

REPORT=reports/inspection/stage0_7h_ats_pooling_s2
mkdir -p "$REPORT"

while IFS=$'\t' read -r cfg run_id; do
  [[ -z "${cfg:-}" ]] && continue
  log "=== train $run_id ==="
  uv run mbs train cascade \
    --config "$cfg" \
    --run-id "$run_id" \
    --device cuda \
    --skip-if-done \
    --report-dir "$REPORT/_staging_${run_id}"
done < configs/experiment/_generated/ats_pooling_s2_manifest.txt

log "=== write compact comparison ==="
uv run python - <<'PY'
import json, statistics
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path('.')
ART = ROOT / 'artifacts' / 'runs'
arms = {
    'P2-G-s2': 'stage0-7h-ats-P2-G-s2',
    'scalar-mean-max-s2': 'stage0-7h-ats-scalar-mean-max-s2',
    'scalar-max-mean-s2': 'stage0-7h-ats-scalar-max-mean-s2',
    'vector-mean-max-s2': 'stage0-7h-ats-vector-mean-max-s2',
}
rows = []
for label, rid in arms.items():
    folds = []
    for i in range(3):
        mp = ART / rid / f'fold_{i}' / 'metrics.json'
        if not mp.is_file():
            continue
        b = json.loads(mp.read_text())
        e = (b.get('evaluations') or {}).get('mbs_e2e') or {}
        m = e.get('metrics') or {}
        folds.append({
            'fold': i,
            'eval_split': e.get('eval_split'),
            'tissue_f1': (m.get('tissue') or {}).get('macro_f1'),
            'age_mae': (m.get('age') or {}).get('mae'),
            'sex_auroc': (m.get('sex') or {}).get('auroc'),
            'best_epoch': b.get('best_epoch') or (b.get('checkpoint_selection') or {}).get('best_epoch'),
        })
    def mean(key):
        vals = [float(f[key]) for f in folds if isinstance(f.get(key), (int, float))]
        return float(statistics.mean(vals)) if vals else None
    def sd(key):
        vals = [float(f[key]) for f in folds if isinstance(f.get(key), (int, float))]
        return float(statistics.stdev(vals)) if len(vals) >= 2 else None
    rows.append({
        'arm': label,
        'run_id': rid,
        'n_folds': len(folds),
        'folds': folds,
        'tissue_f1_mean': mean('tissue_f1'),
        'tissue_f1_sd': sd('tissue_f1'),
        'age_mae_mean': mean('age_mae'),
        'sex_auroc_mean': mean('sex_auroc'),
    })

out_dir = ROOT / 'reports' / 'inspection' / 'stage0_7h_ats_pooling_s2'
out_dir.mkdir(parents=True, exist_ok=True)
summary = {
    'generated_at': datetime.now(timezone.utc).isoformat(),
    'seed': 43,
    'split_id': 'hub-ats-7e-3fold-v1',
    'matrix_id': 'matrix-hub-age-tissue-sex-full-v1',
    'note': 'Genuine second-seed confirmation of cascade 2×2 pooling (Track B.4).',
    'arms': rows,
}
(out_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
lines = [
    '# 7H Track B.4 — ATS cascade pooling seed-43',
    '',
    f"Generated: `{summary['generated_at']}`",
    '',
    'Seed **43** on `hub-ats-7e-3fold-v1` (exercises the seed-offset fix).',
    '',
    '| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |',
    '|-----|------:|----------:|--------:|----------:|',
]
for r in rows:
    t = r['tissue_f1_mean']; a = r['age_mae_mean']; s = r['sex_auroc_mean']
    td = r['tissue_f1_sd']
    t_s = '—' if t is None else (f'{t:.3f}' if td is None else f'{t:.3f} (±{td:.3f})')
    a_s = '—' if a is None else f'{a:.3f}'
    s_s = '—' if s is None else f'{s:.3f}'
    lines.append(f"| {r['arm']} | {r['n_folds']} | {t_s} | {a_s} | {s_s} |")
lines += ['', 'Compare to seed-42 matched 16-ep screen in `stage0_7g_gene_only_probe/analysis.md`.', '']
(out_dir / 'analysis.md').write_text('\n'.join(lines) + '\n')
print('wrote', out_dir / 'analysis.md')
PY

log "=== B.4 done ==="
