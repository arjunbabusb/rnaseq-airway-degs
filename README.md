# RNA-Seq Differential Gene Expression Analysis in R

A beginner-to-intermediate computational biology project analyzing real RNA-seq
data in R — from raw counts to differential expression to biological
interpretation. Built as a learning project while picking up R programming.

## Overview

This project uses the `airway` dataset (Himes et al., 2014, *PLoS ONE*), which
contains RNA-seq counts from human airway smooth muscle cell lines, some
treated with dexamethasone (a glucocorticoid) and some untreated. The goal is
to identify genes whose expression changes significantly with treatment, and
to interpret those genes biologically.

## Learning Objectives

- Core R programming: vectors, data frames, functions, control flow
- Data wrangling and visualization with `dplyr` / `ggplot2`
- RNA-seq data structures, normalization, and quality control
- Differential expression analysis using `DESeq2`
- Statistical concepts: hypothesis testing, p-values, multiple testing correction
- Functional enrichment analysis (GO / KEGG) with `clusterProfiler`
- Reproducible reporting with R Markdown
- Version control and project organization with Git/GitHub

## Dataset

**Package:** [`airway`](https://bioconductor.org/packages/release/data/experiment/html/airway.html) (Bioconductor)
**Samples:** 8 (4 dexamethasone-treated, 4 untreated), 4 cell lines
**Reference:** Himes BE, et al. (2014) *RNA-Seq Transcriptome Profiling
Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates
Cytokine Function in Airway Smooth Muscle Cells.* PLoS ONE 9(6):e99625.

## Project Timeline

### Week 1 — R Foundations & Exploratory Data Analysis
- R basics: vectors, data frames, functions, `dplyr`, `ggplot2`
- Load the `airway` dataset, inspect sample metadata
- Plot library sizes and raw count distributions
- **Output:** `scripts/01_eda.R`, exploratory plots

### Week 2 — Normalization & Quality Control
- Why raw RNA-seq counts aren't directly comparable across samples
- CPM normalization, log transformation
- PCA and sample-correlation heatmap for QC
- **Output:** `scripts/02_normalization_qc.R`, PCA plot, heatmap

### Week 3 — Differential Expression Analysis
- Statistical foundations: hypothesis testing, FDR correction
- Full `DESeq2` workflow
- Volcano plot, MA plot, ranked table of significant genes
- **Output:** `scripts/03_differential_expression.R`, `results/tables/DE_results.csv`

### Week 4 — Functional Enrichment & Reporting
- GO / KEGG enrichment analysis on significant genes
- Heatmap of top differentially expressed genes
- Compile the full analysis into one reproducible R Markdown report
- **Output:** `scripts/04_enrichment_analysis.R`, `report/final_report.Rmd` (+ knitted HTML)

### Week 5 (Optional Extension)
Pick one:
- Build a small Shiny app to interactively explore the volcano plot / gene list
- Repeat the pipeline on a disease-relevant dataset from GEO/TCGA
- Gene co-expression network analysis (WGCNA)

## Repository Structure

```
.
├── README.md
├── setup.R
├── scripts/
│   ├── 01_eda.R
│   ├── 02_normalization_qc.R
│   ├── 03_differential_expression.R
│   └── 04_enrichment_analysis.R
├── report/
│   └── final_report.Rmd
└── results/
    ├── figures/
    └── tables/
```
*(Data is not stored in the repo — it's loaded directly from the `airway`
Bioconductor package via `setup.R`.)*

## How to Run

1. Install R (≥ 4.2) and RStudio
2. Clone this repo and open it as an RStudio project
3. Run `source("setup.R")` once to install all required packages
4. Run the scripts in `scripts/` in order (01 → 04)
5. Knit `report/final_report.Rmd` to generate the full HTML report

## Key R Packages

`tidyverse`, `DESeq2`, `airway`, `pheatmap`, `EnhancedVolcano`,
`clusterProfiler`, `org.Hs.eg.db`, `rmarkdown`

## References

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, et al. (2014)
RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid
Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle
Cells. PLoS ONE 9(6): e99625.
