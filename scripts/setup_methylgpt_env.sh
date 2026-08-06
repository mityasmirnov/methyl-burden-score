#!/usr/bin/env bash
# Create / refresh the MethylGPT-compatible project environment.
#
# MethylGPT requires torchtext, which only ships ABI-matched wheels through
# torch 2.4. The main .venv tracks CpGPT / MBS (torch>=2.3, currently 2.13)
# and cannot import MethylGPT. Use this dedicated env for MethylGPT export work.
#
# Usage:
#   bash scripts/setup_methylgpt_env.sh
#   source scripts/activate_data_environment.sh
#   .venv-methylgpt/bin/python -c "from methylgpt import MethylGPTModel; print('ok')"
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

VENV="$MBS_ROOT/.venv-methylgpt"
PYTHON_BIN="$VENV/bin/python"

if [[ ! -x "$PYTHON_BIN" ]]; then
  uv venv --python 3.11 "$VENV"
fi

# Match vendor/methylgpt/requirements.txt torch/torchtext pins from the PyTorch index
# so the native torchtext extension links correctly.
uv pip install --python "$PYTHON_BIN" \
  --index-url https://download.pytorch.org/whl/cu121 \
  torch==2.1.0 torchtext==0.16.0

uv pip install --python "$PYTHON_BIN" \
  'numpy>=1.24,<2.0' 'pandas>=2.0,<2.3' 'scipy>=1.10,<1.14' \
  'scikit-learn>=1.2,<1.7' 'lightning>=2.0,<2.5' 'tqdm>=4.65' 'datasets>=2.14' \
  gdown scib anndata scanpy ipython matplotlib seaborn

uv pip install --python "$PYTHON_BIN" -e "$MBS_ROOT/vendor/methylgpt" --no-deps

"$PYTHON_BIN" - <<'PY'
from methylgpt import MethylGPTModel, MethylVocab, __version__
import torch
import torchtext

print(f"methylgpt {__version__}")
print(f"torch {torch.__version__}")
print(f"torchtext {torchtext.__version__}")
print(f"python {__import__('sys').executable}")
PY

printf '\nMethylGPT env ready: %s\n' "$VENV"
printf 'Activate data env, then run: %s ...\n' "$PYTHON_BIN"
