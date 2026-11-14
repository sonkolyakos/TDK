# config.sh
# Edit these paths for your environment.

# Root folder containing BAMs and splicing_analysis/rmats/
ROOT="/v/projects/scanalysis/aligned_data/rsem"

# Output folder
RMATS_BASE="/v/projects/scanalysis/rsem_quant"
OUT="$ROOT/clean_event_sashimi_v2"
# rMATS result folders 
T47D_RMATS="$RMATS_BASE/sejtvonal_A/rmats/CTR_vs_TIS/output"
MCF7_RMATS="$RMATS_BASE/sejtvonal_B/rmats/CTR_vs_TIS/output"

# Bulk BAM files
MCF7_CTR1="$ROOT/ELTE-TTK-Onko-1_S123_L001_Aligned.sortedByCoord.out.bam"
MCF7_CTR2="$ROOT/ELTE-TTK-Onko-2_S124_L001_Aligned.sortedByCoord.out.bam"
MCF7_TIS1="$ROOT/ELTE-TTK-Onko-17_S139_L001_Aligned.sortedByCoord.out.bam"
MCF7_TIS2="$ROOT/ELTE-TTK-Onko-18_S140_L001_Aligned.sortedByCoord.out.bam"

T47D_CTR1="$ROOT/SRR32118607_Aligned.sortedByCoord.out.bam"

# Corrected replacement BAM for T47D CTR2
T47D_CTR2="$ROOT/SRR32118609_Aligned.sortedByCoord.out.bam"

T47D_TIS1="$ROOT/SRR32118575_Aligned.sortedByCoord.out.bam"
T47D_TIS2="$ROOT/SRR32118576_Aligned.sortedByCoord.out.bam"

# Default genes. Override at runtime:
#   GENES="RPS24 ZNF207" bash run.sh
GENES_DEFAULT="RPS24 ZNF207 SRSF1 PSMB3 THUMPD2"

# Plotting defaults
FLANK=150
MIN_JUNCTION=5
MAX_JUNCTIONS=4
