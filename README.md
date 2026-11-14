# From Senescence to Repopulation: Alternative Splicing Dynamics in Breast Cancer Cells

Analysis code for the TDK study of the same title (Ákos Gábor Sonkoly, Molecular Bionics Engineering BSc, Pázmány Péter Catholic University, Faculty of Information Technology and Bionics; Institute of Molecular Life Sciences, HUN-REN RCNS; 2026).

The study characterises alternative splicing (AS) in the ER+ breast cancer cell lines **MCF7** and **T47D** across three states: untreated control (**CTR**), therapy-induced senescence (**TIS**) and repopulation (**REPOP**). Gene expression (RSEM + DESeq2) and splicing (rMATS and LeafCutter) are analysed from published bulk RNA-seq data, integrated, and followed across the CTR → TIS → REPOP trajectory.

The study is descriptive. It generates hypotheses about splicing in senescence escape and does not test the function of individual events.

## Data

Bulk RNA-seq from Bajtai et al., *Mol. Cancer* 24, 128 (2025): NCBI GEO [GSE287953](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287953).

- 2 cell lines × 3 conditions (CTR, TIS, REPOP) × 3 biological replicates
- strand-specific, poly(A)-selected, Illumina NovaSeq 6000, ~50 M 150 bp paired-end reads per library
- one T47D control replicate failed QC and was excluded

No raw or processed sequencing data are stored in this repository. Download the reads from GEO/SRA and run the pipeline below.
