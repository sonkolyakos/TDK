#!/usr/bin/env bash
# 01_preprocessing.sh: FastQC -> Trimmomatic -> STAR (2-pass) -> samtools index -> RSEM
# Expects paired FASTQ files: ${RAW_DIR}/<sample>_1.fastq.gz and <sample>_2.fastq.gz
set -euo pipefail

# ---- config (edit) ----------------------------------------------------------
RAW_DIR=${RAW_DIR:-data/fastq}
WORK_DIR=${WORK_DIR:-work}
STAR_INDEX=${STAR_INDEX:-ref/hg38/star_index}              # built from GRCh38 + GENCODE v45
GTF=${GTF:-ref/hg38/gencode.v45.annotation.gtf}
RSEM_REF=${RSEM_REF:-ref/hg38/rsem_index/hg38_rsem}        # rsem-prepare-reference output
ADAPTERS=${ADAPTERS:-$(find "$CONDA_PREFIX" -name TruSeq3-PE.fa | head -n 1)}
THREADS=${THREADS:-15}
# NOTE: the original run used RSEM's default (none). The libraries are strand-specific.
RSEM_STRAND=${RSEM_STRAND:-none}
# -----------------------------------------------------------------------------

QC_DIR=$WORK_DIR/qc;  TRIM_DIR=$WORK_DIR/trimmed
STAR_DIR=$WORK_DIR/star; RSEM_DIR=$WORK_DIR/rsem
mkdir -p "$QC_DIR/raw" "$QC_DIR/trimmed" "$TRIM_DIR" "$STAR_DIR" "$RSEM_DIR"

for r1 in "$RAW_DIR"/*_1.fastq.gz; do
  s=$(basename "$r1" _1.fastq.gz); r2=$RAW_DIR/${s}_2.fastq.gz

  # 1. QC of raw reads
  fastqc -t "$THREADS" "$r1" "$r2" -o "$QC_DIR/raw"

  # 2. Trimming (parameters as in the original run; check against Methods 4.2)
  trimmomatic PE -threads "$THREADS" "$r1" "$r2" \
    "$TRIM_DIR/${s}_1_paired.fastq.gz" "$TRIM_DIR/${s}_1_unpaired.fastq.gz" \
    "$TRIM_DIR/${s}_2_paired.fastq.gz" "$TRIM_DIR/${s}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:25:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25 \
    2> "$TRIM_DIR/${s}_trim.log"

  # 3. QC of trimmed reads
  fastqc -t "$THREADS" "$TRIM_DIR/${s}_1_paired.fastq.gz" "$TRIM_DIR/${s}_2_paired.fastq.gz" \
    -o "$QC_DIR/trimmed"

  # 4. STAR alignment (skip if finished)
  if [ ! -f "$STAR_DIR/${s}_Log.final.out" ]; then
    rm -rf "$STAR_DIR/${s}__STARtmp"
    STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" \
      --readFilesIn "$TRIM_DIR/${s}_1_paired.fastq.gz" "$TRIM_DIR/${s}_2_paired.fastq.gz" \
      --readFilesCommand zcat --outFileNamePrefix "$STAR_DIR/${s}_" \
      --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 25000000000 \
      --twopassMode Basic --outSAMunmapped Within \
      --quantMode TranscriptomeSAM --sjdbGTFfile "$GTF"
  fi
  samtools index -@ "$THREADS" "$STAR_DIR/${s}_Aligned.sortedByCoord.out.bam"

  # 5. RSEM gene/isoform quantification from the transcriptome BAM
  if [ ! -f "$RSEM_DIR/${s}.genes.results" ]; then
    rsem-calculate-expression --alignments --paired-end --no-bam-output \
      --strandedness "$RSEM_STRAND" --num-threads "$THREADS" \
      "$STAR_DIR/${s}_Aligned.toTranscriptome.out.bam" "$RSEM_REF" "$RSEM_DIR/$s"
  fi
done
