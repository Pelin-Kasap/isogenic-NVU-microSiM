
Code to reproduce the RNA-seq analyses (Fig. S1, S2) and the statistical analyses (Fig. 3D, 4B, 4D, 4F, S4) of the manuscript.
Title and DOI will be added upon publication.

## Repository structure

```
rnaseq/
├── 20261014_BMEC_FigS1_final.R
├── 20261014_astrocyte_FigS2_final.R
└── data/
    ├── counts_astrocytes.tsv
    └── counts_primary_BMEC_line.tsv
statistics/
├── statistics_main.R
├── statistics_supplement.R
└── data/
    └── Fig3D_long.csv
```

| Script | Reproduces |
|---|---|
| `rnaseq/20261014_BMEC_FigS1_final.R` | Fig. S1: EECM-BMEC RNA-seq (marker heatmaps, differential expression, volcano plot, top-50 DEG heatmap, GO enrichment) |
| `rnaseq/20261014_astrocyte_FigS2_final.R` | Fig. S2: astrocyte RNA-seq (same analyses) |
| `statistics/statistics_main.R` | Statistics for Fig. 3D, 4B, 4D, 4F |
| `statistics/statistics_supplement.R` | Statistics for Fig. S4A–G |

## Input data

### RNA-seq (`rnaseq/data/`)
Gene-level counts were generated with snakePipes v3.4.0 (mRNAseq workflow): adapter and quality trimming with cutadapt v4.6, alignment with STAR v2.7.10b to GRCh38 (Ensembl primary assembly) with GENCODE release 31 annotation, duplicate marking with sambamba v0.8.0, and gene-level counting with featureCounts (Subread v2.1.1; `-C -Q 10 --primary`, plus `-p -B` for paired-end data). Quality control: MultiQC v1.30.

| File | Samples | Source |
|---|---|---|
| `counts_astrocytes.tsv` | hiPSC-derived astrocytes (5 donor lines) and primary fetal astrocytes (3) | reprocessed from GEO GSE120411 (runs SRR7901104, SRR7901110, SRR7901116, SRR7901128, SRR7901140, SRR7901146, SRR7901148, SRR7901150) |
| `counts_primary_BMEC_line.tsv` | primary human BMEC reference | reprocessed from GEO GSE97575 (GSM2572306, SRR5441866) |
| `counts_EECM-BMECs.tsv` | EECM-BMECs (5 donor lines) | GEO GSE346643; **will be added to this repository upon publication** |

### Statistics (`statistics/data/`)
All values are provided in the Source Data file of the manuscript, where technical replicates (devices) are grouped by independent experiment or donor line.
`Fig3D_long.csv` is included as an example of the input format (one row per device, with its experimental group and independent experiment or donor line). The corresponding files for the remaining figures (`Fig4B_long.csv`, `Fig4D_long.csv`, `Fig4F_long.csv`, `FigS4AB_long.csv`, `FigS4C_long.csv`, `FigS4D_long.csv`, `FigS4FG_long.csv`) will be added upon publication.

## How to run
In R, set the working directory to the `rnaseq` or `statistics` folder and source the script, e.g.

```r
setwd("path/to/repository/statistics")
source("statistics_main.R")
```

RNA-seq scripts write figures and DEG tables to `rnaseq/results/<date>_<celltype>/`, together with a `sessionInfo.txt` listing all package versions.

## Software
- RNA-seq: R v4.5.2, DESeq2 v1.50.2, clusterProfiler v4.18.4, GO.db v3.22.0, GOSemSim v2.36.0, org.Hs.eg.db v3.22.0, together with dplyr, ggplot2, ggrepel, patchwork, pheatmap, dendextend, stringr, openxlsx and magick.
- Statistics: R v4.5.2, lme4 v2.0-1, lmerTest v3.2-1, pbkrtest v0.5.5, emmeans v2.0.4, car v3.1-5, rstatix v1.1.0, dplyr v1.2.1. 

## Contact
Pelin Kasap 
pelin.kasap@unibas.ch
