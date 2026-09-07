#!/usr/bin/env bash
# Track B.4: genuine seed-43 2×2 cascade pooling confirmation on ATS.
# Do NOT run while nine-pack full refs own GPU 0
# (scratch/logs/7h_nine_pack_full.pid). Uses experiment.seed=43 so the
# flat/cascade seed-offset fix is exercised.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1

LOG=scratch/logs/7h_ats_pooling_s2.log
mkdir -p scratch/logs

# Prefer existing P2/scalar/vector YAMLs with seed override via temp copies.
python3 - <<'PY'
from pathlib import Path
import yaml
root = Path('.')
arms = {
    'p2': 'configs/experiment/stage0_7g_gene_only_probe_p2.yaml',
    'mean_max': 'configs/experiment/stage0_7g_gene_only_probe_scalar_mean_max.yaml',
    'max_mean': 'configs/experiment/stage0_7g_gene_only_probe_scalar_max_mean.yaml',
    'vector_mean_max': 'configs/experiment/stage0_7g_gene_only_probe_vector_mean_max.yaml',
}
out = Path('configs/experiment/_generated')
out.mkdir(parents=True, exist_ok=True)
for name, src in arms.items():
    p = Path(src)
    if not p.is_file():
        print(f'skip missing {src}')
        continue
    cfg = yaml.safe_load(p.read_text())
    cfg.setdefault('experiment', {})['seed'] = 43
    cfg['experiment']['name'] = str(cfg['experiment'].get('name', name)) + '_s2'
    dest = out / f'stage0_7h_ats_{name}_s2.yaml'
    dest.write_text(yaml.safe_dump(cfg, sort_keys=False))
    print('wrote', dest)
PY

echo "[ats-s2] configs ready under configs/experiment/_generated/"
echo "[ats-s2] launch example (after GPU free):"
echo "  uv run mbs train cascade --config configs/experiment/_generated/stage0_7h_ats_p2_s2.yaml --run-id stage0-7h-ats-P2-G-s2 --device cuda --skip-if-done"
