#!/bin/bash
set -euo pipefail

if ! command -v conda >/dev/null 2>&1; then
    echo "conda: MISSING"
    exit 1
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate clean_sashimi

python - <<'PY'
mods = ["numpy", "pandas", "matplotlib", "pysam"]
for m in mods:
    mod = __import__(m)
    print(f"{m}: {getattr(mod, '__version__', 'OK')}")
PY

samtools --version | head -1
