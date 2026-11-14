#!/bin/bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"

mkdir -p "$OUT"

# Conda activation
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        source "$CONDA_BASE/etc/profile.d/conda.sh"
    fi
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available. Run ./install.sh first or install Miniconda/Anaconda."
    exit 1
fi

conda activate clean_sashimi

PY="$HERE/make_clean_event_sashimi_v2.py"
GENES="${GENES:-$GENES_DEFAULT}"

# Validate inputs
for f in \
    "$MCF7_CTR1" "$MCF7_CTR2" "$MCF7_TIS1" "$MCF7_TIS2" \
    "$T47D_CTR1" "$T47D_CTR2" "$T47D_TIS1" "$T47D_TIS2"
do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing BAM:"
        echo "  $f"
        exit 1
    fi
done

for d in "$MCF7_RMATS" "$T47D_RMATS"
do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: missing rMATS folder:"
        echo "  $d"
        exit 1
    fi
done

echo -e "cellline\tgene\tevent_type\tregion\tdeltaPSI\tFDR\toutput" \
    > "$OUT/index.tsv"

for CELL in MCF7 T47D
do
    if [[ "$CELL" == "MCF7" ]]; then
        RMATS="$MCF7_RMATS"
        CTR1="$MCF7_CTR1"
        CTR2="$MCF7_CTR2"
        TIS1="$MCF7_TIS1"
        TIS2="$MCF7_TIS2"
    else
        RMATS="$T47D_RMATS"
        CTR1="$T47D_CTR1"
        CTR2="$T47D_CTR2"
        TIS1="$T47D_TIS1"
        TIS2="$T47D_TIS2"
    fi

    for GENE in $GENES
    do
        echo
        echo "========================================"
        echo "$CELL $GENE"
        echo "========================================"

        python "$PY" \
            --cellline "$CELL" \
            --gene "$GENE" \
            --rmats "$RMATS" \
            --ctr "$CTR1" "$CTR2" \
            --tis "$TIS1" "$TIS2" \
            --out "$OUT/${CELL}_${GENE}_clean_v2" \
            --flank "$FLANK" \
            --min-junction "$MIN_JUNCTION" \
            --max-junctions "$MAX_JUNCTIONS" \
            | tee -a "$OUT/index.tsv" || true
    done
done

echo
echo "DONE"
echo "Output:"
echo "  $OUT"
