# ============================================================
# setup.R
# Run this ONCE before starting Week 1.
# Installs every package used across the whole project.
# ============================================================

installed <- rownames(installed.packages())

# ---- 1. CRAN packages -----------------------------------------------------
cran_packages <- c("tidyverse", "pheatmap", "ggrepel", "knitr", "rmarkdown")
to_install <- setdiff(cran_packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install)
}

# ---- 2. Bioconductor packages ----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("DESeq2", "airway", "clusterProfiler", "org.Hs.eg.db", "EnhancedVolcano")
bioc_to_install <- setdiff(bioc_packages, installed)
if (length(bioc_to_install) > 0) {
  BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
}

# ---- 3. Verify everything loads --------------------------------------------
library(tidyverse)
library(DESeq2)
library(airway)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

cat("\n✅ All packages installed and loaded successfully. You're ready for Week 1!\n")

# ---- Notes ------------------------------------------------------------------
# - This can take 15-30 minutes the first time (Bioconductor packages are large).
# - If a package fails to install, re-run just that line individually -
#   error messages usually tell you exactly what's missing.
# - Ignore the earlier "dplyr/tidyr not installed" warning you saw in RStudio -
#   installing "tidyverse" above covers dplyr, tidyr, ggplot2, and more at once.
install.packages("gprofiler2")
library(gprofiler2)
sig_gene <- res_annotated[res_annotated$sig %in% c("Up", "Down"), ]
