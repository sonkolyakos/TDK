#!/usr/bin/env bash
# 03a_leafcutter.sh: junctions -> clustering -> differential splicing (per cell line)
# Requires: regtools, samtools, python3, R (4.2) with the leafcutter package
set -euo pipefail

LEAFCUTTER_DIR=${LEAFCUTTER_DIR:-tools/leafcutter}     # git clone davidaknowles/leafcutter
STAR_DIR=${STAR_DIR:-work/star}                        # <sample>_Aligned.sortedByCoord.out.bam
OUT=${OUT:-work/leafcutter}
THREADS=${THREADS:-10}
export PATH="$LEAFCUTTER_DIR/scripts:$PATH"

# sample:group  (T47D SRR32118608 excluded after QC; its junctions are still extracted, as in the original run)
MCF7="ELTE-TTK-Onko-1_S123_L001:CTR ELTE-TTK-Onko-2_S124_L001:CTR ELTE-TTK-Onko-3_S125_L001:CTR \
ELTE-TTK-Onko-17_S139_L001:TIS ELTE-TTK-Onko-18_S140_L001:TIS ELTE-TTK-Onko-19_S141_L001:TIS \
ELTE-TTK-Onko-34_S156_L001:REPOP ELTE-TTK-Onko-35_S157_L001:REPOP SRR32118588:REPOP"
T47D="SRR32118607:CTR SRR32118609:CTR SRR32118575:TIS SRR32118576:TIS SRR32118577:TIS \
SRR32118583:REPOP SRR32118584:REPOP SRR32118585:REPOP"
T47D_EXTRA_JUNC="SRR32118608"

run_cell () {
  local cell=$1 samples=$2 extra=${3:-}
  local dir=$OUT/$cell; mkdir -p "$dir"; cd "$dir"
  : > groups.txt; : > junc_files.txt

  for sg in $samples; do printf '%s\t%s\n' "${sg%%:*}" "${sg##*:}" >> groups.txt; done
  for s in $(cut -f1 groups.txt) $extra; do
    bam2junc.sh "$STAR_DIR/${s}_Aligned.sortedByCoord.out.bam" "${s}.junc"
    echo "${s}.junc" >> junc_files.txt
  done

  # clustering: Methods 4.5 (-r 0.001 is the default, set explicitly)
  python "$LEAFCUTTER_DIR/clustering/leafcutter_cluster.py" \
    -j junc_files.txt -m 50 -r 0.001 -l 500000 -o "$cell"

  # differential splicing: >=3 samples per group, >=1 sample per intron (Methods 4.5)
  Rscript "$LEAFCUTTER_DIR/scripts/leafcutter_ds.R" \
    --num_threads "$THREADS" \
    --min_samples_per_intron 1 --min_samples_per_group 3 \
    -o "${cell}_ds_results" \
    "${cell}_perind_numers.counts.gz" groups.txt
  cd - > /dev/null
}

run_cell mcf7 "$MCF7"
run_cell t47d "$T47D" "$T47D_EXTRA_JUNC"
