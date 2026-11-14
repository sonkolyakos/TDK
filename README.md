# Clean event-centric sashimi plotting package

This package creates compact, publication-style sashimi plots from bulk RNA-seq BAM files and rMATS results.

## What the plot contains

- averaged CTR coverage from two replicates
- averaged TIS coverage from two replicates
- only the strongest local splice junctions
- thin, light junction arcs
- staggered junction-count labels to avoid overlap
- minimal event-specific exon schematic
- intron compression
- rMATS event type, ΔPSI and FDR in the title

The script supports:

- SE
- A5SS
- A3SS
- MXE
- RI

## Files

- `make_clean_event_sashimi_v2.py` — plotting and rMATS event-selection code
- `run.sh` — main runner
- `config.sh` — input paths and plotting parameters
- `environment.yml` — conda environment
- `install.sh` — one-command installation
- `check_installation.sh` — dependency check

## Installation

Unpack the archive:

```bash
tar -xzf clean_event_sashimi_package.tar.gz
cd clean_event_sashimi_package
```

Install the environment:

```bash
bash install.sh
```

No separate `ggsashimi` installation is required. This package uses its own event-centric Python plotting code.

## Configuration

Open:

```bash
nano config.sh
```

Change `ROOT`, BAM paths, and rMATS folders if your directory structure differs.

The supplied config reflects the original analysis layout:

```text
ROOT=/spisaknas/kisscs/spisaknas_data/senescence_bulk
```

Important: the supplied T47D CTR2 path points to the corrected S128 replacement BAM rather than the older S127 file.

## Run

```bash
bash run.sh
```

Default genes:

```text
RPS24
ZNF207
HNRNPA2B1
SRSF1
```

Run only selected genes:

```bash
GENES="RPS24 ZNF207" bash run.sh
```

## Output

By default:

```text
$ROOT/clean_event_sashimi_v2/
```

Each plot is written as both PDF and PNG.

An `index.tsv` file records:

- cell line
- gene
- selected rMATS event type
- genomic region
- ΔPSI
- FDR
- output path

## Event selection

For each gene, the script searches rMATS `*.MATS.JC.txt` files and requires:

```text
FDR <= 0.05
|ΔPSI| >= 0.10
```

Among qualifying events it prioritizes:

1. lower FDR
2. larger absolute ΔPSI

## Plot parameters

Defaults in `config.sh`:

```text
FLANK=150
MIN_JUNCTION=5
MAX_JUNCTIONS=4
```

To reduce clutter further:

```bash
MIN_JUNCTION=10
MAX_JUNCTIONS=3
```

can be set in `config.sh`.

## Interpretation

The plots summarize bulk RNA-seq splice-junction evidence and rMATS event statistics.

They should not be interpreted as transcript-isoform abundance estimates.

## Reproducibility

The environment includes:

- Python 3.11
- numpy
- pandas
- matplotlib
- pysam
- samtools
