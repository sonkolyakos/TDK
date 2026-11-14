# From Senescence to Repopulation: Alternative Splicing Dynamics in Breast Cancer Cells

Analysis code for the TDK study of the same title (Ákos Gábor Sonkoly, Molecular Bionics Engineering BSc, Pázmány Péter Catholic University, Faculty of Information Technology and Bionics; Epigenetics and Genome Editing Laboratory, Institute of Molecular Life Sciences, HUN-REN RCNS; 2026). The thesis is included in this repository: [`TDK_Sonkoly_Ákos_Gábor_2026.pdf`](TDK_Sonkoly_Ákos_Gábor_2026.pdf).

The study characterises alternative splicing (AS) in the ER+ breast cancer cell lines **MCF7** and **T47D** across three states: untreated control (**CTR**), therapy-induced senescence (**TIS**) and repopulation (**REPOP**). Gene expression (RSEM + DESeq2) and splicing (rMATS and LeafCutter) are analysed from published bulk RNA-seq data, integrated, and followed across the CTR → TIS → REPOP trajectory.

The study is descriptive. It generates hypotheses about splicing in senescence escape and does not test the function of individual events.

## Data

Bulk RNA-seq from Bajtai et al., *Mol. Cancer* 24, 128 (2025): NCBI GEO [GSE287953](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287953).

- 2 cell lines × 3 conditions (CTR, TIS, REPOP) × 3 biological replicates
- strand-specific, poly(A)-selected, Illumina NovaSeq 6000, ~50 M 150 bp paired-end reads per library
- one T47D control replicate (SRR32118608) failed QC and was excluded from the differential analyses

[`metadata_bulk.csv`](metadata_bulk.csv) is the sample sheet (sample identifiers, cell line, condition, replicate).

No raw or processed sequencing data are stored here. Download the reads from GEO/SRA and run the pipeline below.

## Repository layout

```
README.md
TDK_Sonkoly_Ákos_Gábor_2026.pdf     thesis
metadata_bulk.csv                    sample sheet
common_genes_multiDB_ORA.R           multi-database over-representation analysis (Figure 11)
bulk_analysis/                       preprocessing, DGE, differential splicing, integration
sashimi/                             event-centric sashimi plots (own README)
```

Contents of `bulk_analysis/` by analysis step (file names may differ):

| Step | Content |
|---|---|
| Preprocessing | FastQC, Trimmomatic, STAR 2-pass alignment, RSEM quantification |
| DGE / GSEA | DESeq2 (apeglm), PCA and correlation QC, KEGG GSEA |
| Differential splicing | rMATS runs; LeafCutter junction extraction, clustering and differential splicing |
| `rmats_overlap.Rmd` | global splicing landscape, event-type proportions, overlap between comparisons and cell lines |
| global ΔPSI controls | ΔPSI distribution and sign asymmetry of all tested events |
| `integration.Rmd` | rMATS ∩ LeafCutter, consensus genes, trajectories, concordance, cross-cell-line overlap |
| DAS vs DGE | overlap of DAS and DGE genes, ΔPSI vs log2FC, cell-line concordance |

## Pipeline overview

1. **Preprocessing and quantification.** Reads from the SRA Toolkit, FastQC before and after trimming, Trimmomatic (`ILLUMINACLIP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`), STAR 2.7.11b against GRCh38 (GENCODE v45), RSEM 1.3.3 for gene-level expression.
2. **Differential gene expression.** DESeq2 per cell line (`~condition`), all three contrasts (CTR vs TIS, TIS vs REPOP, CTR vs REPOP), apeglm-shrunken log2FC. DGE: BH-adjusted p < 0.05 and |log2FC| > 1. GSEA (KEGG) with clusterProfiler.
3. **Differential splicing.**
   - **rMATS** (event level: SE, RI, A5SS, A3SS, MXE; junction-count output `JC`). Events need a mean inclusion+skipping junction count ≥ 30 and are significant at FDR < 0.05 and |ΔPSI| > 0.1.
   - **LeafCutter** (intron level): `bam2junc.sh` (regtools), `leafcutter_cluster.py -m 50 -r 0.001 -l 500000`, `leafcutter_ds.R`. Introns are significant at FDR < 0.05 and |ΔPSI| > 0.1. The FDR is a cluster-level value from a single model with the three conditions, so it is shared by the three comparisons; only the ΔPSI differs between them.
4. **Integration.** The tools are matched by genomic coordinates. A gene is *differentially alternatively spliced (DAS)* if it is significant in rMATS and in LeafCutter for the same comparison (not necessarily through the same event). Directional concordance is evaluated on rMATS-significant SE events against unfiltered LeafCutter values.
5. **Trajectories.** Each CTR–TIS-significant event is classified by r = ΔPSI(REPOP−CTR) / ΔPSI(TIS−CTR): *reversed* (r < −0.33), *transient* (−0.33 ≤ r < 0.33), *partial* (0.33 ≤ r < 0.67), *persistent* (0.67 ≤ r ≤ 1.33), *progressive* (r > 1.33). Events not significant in CTR–TIS but significant in CTR–REPOP are *REPOP-specific*.
6. **Functional context.** Over-representation analysis (MSigDB C2:CP canonical pathways, GO:BP, Hallmark) of genes with shared splicing changes, DAS vs DGE comparison, and representative loci (RPS24, ZNF207, SRSF1, PSMB3, THUMPD2).

### Conventions

- **ΔPSI = second condition − first condition**: CTR_vs_TIS = TIS − CTR, TIS_vs_REPOP = REPOP − TIS, CTR_vs_REPOP = REPOP − CTR.
- rMATS reports `IncLevelDifference` = group 1 − group 2, so the loaders multiply it by −1 (`invert_dpsi = TRUE`). This is correct only if the first-named condition is `--b1` in every rMATS run (check the run commands or `b1.txt`).
- The LeafCutter ΔPSI of a skipping intron is multiplied by −1 when compared with rMATS inclusion levels.


## Figure → script map

Figure numbers follow the thesis.

| Figure | Content | Script |
|---|---|---|
| 4–5 | DGE volcano plots, GSEA dot plots | `bulk_analysis/` (DGE / GSEA) |
| 6 | ΔPSI distributions by event type | `rmats_overlap.Rmd` |
| 7A–C | event-type proportions, event overlap | `rmats_overlap.Rmd` |
| 8 | trajectory heatmaps | `integration.Rmd` |
| 9 | UpSet plots, rMATS vs LeafCutter | `integration.Rmd` |
| 10 | rMATS vs LeafCutter directional concordance | `integration.Rmd` |
| 11 | multi-database ORA | `common_genes_multiDB_ORA.R` |
| 12–14 | DAS vs DGE, ΔPSI vs log2FC, cross-cell-line concordance | `bulk_analysis/` (DAS vs DGE) |
| 15–20 | sashimi plots and splice-composition bar plots | `sashimi/` |
| 17 | mean \|ΔPSI\| trajectory profiles | `integration.Rmd` |

## Running the analysis

The analysis scripts read rMATS and LeafCutter results from directories set at the top of each file (`RMATS_BASE_DIR`, `LEAF_BASE_DIR`). Expected layout:

```
rmats/{mcf7,t47d}/{CTR_vs_TIS,TIS_vs_REPOP,CTR_vs_REPOP}/output/{SE,RI,MXE,A5SS,A3SS}.MATS.JC.txt
leafcutter/{mcf7,t47d}_ds_results_{cluster_significance,effect_sizes}.txt, mcf7_annotated_effect_sizes.txt
```

Change these paths to your own before running.

Suggested order: preprocessing → DGE/GSEA → rMATS and LeafCutter → `rmats_overlap.Rmd` → `integration.Rmd` → ORA:

```bash
Rscript common_genes_multiDB_ORA.R common_genes_MCF7_T47D_core_pan.csv ORA_out tested_splicing_universe.tsv
```

The ORA inputs are produced by `integration.Rmd`: the consensus genes shared by both cell lines (`core_pan`) and the genes testable for splicing in both cell lines (`univ_pan`, used as the background universe).

## Software versions

| Tool | Version |
|---|---|
| SRA Toolkit | 3.4.1 |
| FastQC | 0.12.1 |
| Trimmomatic | 0.41 |
| STAR | 2.7.11b |
| RSEM | 1.3.3 |
| R | 4.2.3 |
| DESeq2 | 1.38.3 (apeglm shrinkage) |
| clusterProfiler | 4.6.2 |
| rMATS | 4.1.2 (`-t paired --readLength 150 --variable-read-length --cstat 0.0001`, annotated splice sites only) |
| regtools | 1.0.0 |
| LeafCutter | `leafcutter_cluster.py -m 50 -r 0.001 -l 500000`; `leafcutter_ds.R` (≥ 3 samples per group, ≥ 1 per intron) |
| Python | 3.11 (numpy, pandas, matplotlib, pysam) for the sashimi plots |

Other R packages: ComplexHeatmap, circlize, eulerr, ggplot2, patchwork, msigdbr, org.Hs.eg.db, data.table.

## Notes and limitations

- Bulk RNA-seq with three replicates per condition from two ER+ luminal cell lines; results may not generalise to other subtypes or to tumours in vivo.
- Short reads report relative exon/intron usage, not full-length isoforms or protein levels.
- Trajectory classes are based on ΔPSI point estimates. A non-significant CTR–REPOP difference is not proof of identity, and the *reversed* class is sensitive to small initial changes.
- Bulk data cannot distinguish within-cell changes from changes in population composition.
- Observational study: no causal role of any splicing event is claimed. Proposed follow-up: RT-PCR validation, long-read sequencing, single-cell data and perturbation of candidate regulators such as SRSF1.

## Citation

> Sonkoly, Ákos G. *From Senescence to Repopulation: Alternative Splicing Dynamics in Breast Cancer Cells.* TDK study, Pázmány Péter Catholic University, 2026.
>
> Bajtai, E. et al. Therapy-induced senescence is a transient drug resistance mechanism in breast cancer. *Mol. Cancer* 24, 128 (2025). doi:10.1186/s12943-025-02310-0

