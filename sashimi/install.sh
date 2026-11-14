#!/bin/bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available."
    echo "Install Miniconda/Anaconda first, then rerun:"
    echo "  bash install.sh"
    exit 1
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -Fxq clean_sashimi; then
    echo "Environment clean_sashimi already exists."
    echo "Updating from environment.yml..."
    conda env update \
        -n clean_sashimi \
        -f "$HERE/environment.yml" \
        --prune
else
    echo "Creating clean_sashimi environment..."
    conda env create \
        -f "$HERE/environment.yml"
fi

conda activate clean_sashimi

python - <<'PY'
import numpy
import pandas
import matplotlib
import pysam

print("Environment OK")
print("numpy      ", numpy.__version__)
print("pandas     ", pandas.__version__)
print("matplotlib ", matplotlib.__version__)
print("pysam      ", pysam.__version__)
PY

echo
echo "Installation finished."
echo "Next:"
echo "  1. edit config.sh if paths differ"
echo "  2. bash run.sh"
