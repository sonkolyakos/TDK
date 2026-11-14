# Event-centric sashimi plots

Compact, publication-style sashimi plots for individual alternative splicing events, built from bulk RNA-seq BAM files and rMATS results. Used for the representative-event figures (Sashimi panels) of the TDK study *From Senescence to Repopulation: Alternative Splicing Dynamics in Breast Cancer Cells*.

## What the plot contains

- mean CTR coverage across replicates
- mean TIS coverage across replicates
- junction arcs with counts (skipped-exon events only, see below)
- minimal event-specific exon schematic
- intron compression
- rMATS event type, ΔPSI and FDR in the title

Supported event types (for event selection, coverage and schematic): `SE`, `A5SS`, `A3SS`, `MXE`, `RI`.

## Files

| File | Purpose |
|---|---|
| `make_clean_event_sashimi_v3.py` | plotting and rMATS event-selection code |
| `run.sh` | main runner |
| `config.sh` | input paths and plotting parameters |
| `environment.yml` | conda environment |
| `install.sh` | one-command installation |
| `check_installation.sh` | dependency check |

## Installation

```bash
bash install.sh
bash check_installation.sh
```

No separate `ggsashimi` installation is needed; the package uses its own Python plotting code.

Environment: Python 3.11, numpy, pandas, matplotlib, pysam, samtools.

## Configuration

Edit `config.sh` and set `ROOT`, the BAM paths and the rMATS result folders to match your directory layout:

```bash
nano config.sh
```

Use the BAM files of the replicates of each condition (CTR, TIS) and the rMATS output folder of the **CTR vs TIS** comparison, with CTR as `--b1`. The sign of ΔPSI shown in the title assumes this: the script multiplies rMATS `IncLevelDifference` by −1, so the plotted value is TIS − CTR only if CTR was the first group.

## Run

```bash
bash run.sh
```

Default genes: `RPS24`, `ZNF207`, `HNRNPA2B1`, `SRSF1`.

Run selected genes only:

```bash
GENES="RPS24 ZNF207" bash run.sh
```

Direct call of the plotting script:

```bash
python make_clean_event_sashimi_v3.py \
  --cellline MCF7 --gene ZNF207 \
  --rmats <path/to/rmats/CTR_vs_TIS/output> \
  --ctr <CTR_rep1.bam> <CTR_rep2.bam> <CTR_rep3.bam> \
  --tis <TIS_rep1.bam> <TIS_rep2.bam> <TIS_rep3.bam> \
  --preferred-type SE \
  --out figs/MCF7_ZNF207
```

## Output

Each plot is written as PDF and PNG to the directory given by `--out` (by default a subfolder of `$ROOT`). `index.tsv` records cell line, gene, selected rMATS event type, genomic region, ΔPSI, FDR and output path.

## Event selection

For each gene the script reads the rMATS `*.MATS.JC.txt` files and keeps events that satisfy all of:

```text
mean IJC+SJC >= MIN_READS in both conditions   (default 20)
FDR          <= 0.05
|ΔPSI|       >= 0.10
```

Among qualifying events it selects the one with

1. the lowest FDR, then
2. the largest |ΔPSI|.

`--preferred-type` restricts the search to one event type first (e.g. `SE`).

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--flank` / `FLANK` | 150 | bases of flanking sequence around the event |
| `--min-reads` | 20 | minimum mean IJC+SJC reads per condition (removes low-depth events) |
| `--preferred-type` | none | event type to search first |
| `--min-junction`, `--max-junctions` | – | accepted for backward compatibility with older runner scripts, **ignored** in v3 |

## Notes on what is drawn

- **Coverage** is the mean per-base read depth across the replicates of each condition, computed from the BAM files. It is **not** normalised to library size (no RPM scaling), so the y-axis is not directly comparable between CTR and TIS.
- **Junction arcs and counts** are taken directly from the rMATS `IJC`/`SJC` columns (mean across replicates), not recounted from the BAM files. The inclusion count is split equally between the upstream and downstream inclusion junctions. Arcs are drawn for **SE events only**; for other event types only coverage and the schematic are shown.
- The selected event is the single best-supported rMATS event of the gene, which may differ from other significant events of the same gene.

## Interpretation

The plots summarise bulk RNA-seq splice-junction evidence and rMATS event statistics. They are not transcript-isoform abundance estimates.

## Reproducibility

Raw data: NCBI GEO, accession [GSE287953](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287953). Alignments were generated with STAR 2.7.11b (GRCh38, GENCODE v45) and splicing was quantified with rMATS 4.1.2 (see the main analysis repository).
