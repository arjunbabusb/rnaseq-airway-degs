# =========================================================
# RNA-seq Differential Expression Analysis
# Dataset: airway (Himes et al., 2014, PLoS ONE)
# Run this top to bottom (Source) any time your R session
# resets — it regenerates every object from scratch.
# =========================================================

library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)
library(gprofiler2)

# ---- Load data ----
data(airway)

# ---- Week 2: Build DESeqDataSet, filter low-count genes ----
dds <- DESeqDataSetFromMatrix(
  countData = assay(airway),
  colData   = colData(airway),
  design    = ~ cell + dex
)

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]̥

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# ---- Week 2: Variance-stabilizing transform + QC plots ----
vsd <- vst(dds, blind = TRUE)

plotPCA(vsd, intgroup = "dex")

sampleDists       <- dist(t(assay(vsd)))
sampleDistMatrix  <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix)

# ---- Week 3: Differential expression test ----
dds$dex <- relevel(dds$dex, ref = "untrt")
dds     <- DESeq(dds)
res     <- results(dds)
summary(res)

# ---- Week 3: Rank by significance + volcano plot ----
res <- res[order(res$padj), ]

res_df      <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df$sig <- "Not significant"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey")) +
  theme_minimal()

# ---- Week 3: Annotate with real gene symbols ----
gene_symbols <- as.data.frame(rowData(airway))[, c("gene_id", "gene_name")]

res_annotated <- merge(res_df, gene_symbols, by.x = "gene", by.y = "gene_id")
res_annotated <- res_annotated[order(res_annotated$padj), ]

# Sanity check: CRISPLD2 is the headline gene from the original publication
res_annotated[res_annotated$gene_name == "CRISPLD2", ]

write.csv(res_annotated, "results_DE_annotated.csv", row.names = FALSE)

# ---- Week 4: Heatmap of top 20 DE genes ----
top_hits <- head(res_annotated, 20)
mat      <- assay(vsd)[top_hits$gene, ]
rownames(mat) <- top_hits$gene_name

pheatmap(mat,
         annotation_col = as.data.frame(colData(vsd))["dex"],
         scale = "row")

# ---- Week 4: Pathway enrichment (gprofiler2 — no org.Hs.eg.db needed) ----
sig_genes <- res_annotated[res_annotated$sig %in% c("Up", "Down"), ]
gene_list <- sig_genes$gene_name

gostres <- gost(query = gene_list, organism = "hsapiens", sources = c("GO:BP", "KEGG"))
head(gostres$result)
install.packages("gprofiler2")
library(gprofiler2)

sig_genes <- res_annotated[res_annotated$sig %in% c("Up", "Down"), ]
gene_list <- sig_genes$gene_name

gostres <- gost(query = gene_list, organism = "hsapiens", sources = c("GO:BP", "KEGG"))
head(gostres$result)
# More specific GO:BP terms (excluding huge, generic categories)
go_bp <- gostres$result[gostres$result$source == "GO:BP", ]
go_bp_specific <- go_bp[go_bp$term_size < 500, ]
go_bp_specific <- go_bp_specific[order(go_bp_specific$p_value), ]
head(go_bp_specific[, c("term_name", "p_value", "term_size", "intersection_size")])

# KEGG pathways specifically
kegg_hits <- gostres$result[gostres$result$source == "KEGG", ]
kegg_hits <- kegg_hits[order(kegg_hits$p_value), ]
head(kegg_hits[, c("term_name", "p_value", "term_size", "intersection_size")])
