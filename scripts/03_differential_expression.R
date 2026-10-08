# TCGA-BRCA differential expression analysis
# Compare primary tumor and solid tissue normal samples using DESeq2

library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)

# Load processed counts and metadata
load("data/brca_filtered_counts.RData")
load("data/tcga_brca_100samples.RData")

sample_metadata <- as.data.frame(colData(brca_se))

sample_metadata$condition <- ifelse(
  sample_metadata$shortLetterCode == "TP",
  "Tumor",
  "Normal"
)

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("Normal", "Tumor")
)

# Verify sample order
stopifnot(
  all(rownames(sample_metadata) == colnames(counts_filtered))
)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~ condition
)

# Differential expression analysis
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

res <- res[order(res$padj), ]

# Annotate genes
gene_annotation <- as.data.frame(rowData(brca_se))

gene_annotation <- gene_annotation[
  ,
  c("gene_id", "gene_name", "gene_type")
]

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_annotated <- merge(
  res_df,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

# Save differential expression results
dir.create("data", showWarnings = FALSE)

write.csv(
  res_df,
  "data/BRCA_differential_expression.csv",
  row.names = FALSE
)

write.csv(
  res_annotated,
  "data/BRCA_differential_expression_annotated.csv",
  row.names = FALSE
)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA
pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = condition)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of TCGA-BRCA RNA-seq Samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Sample Type"
  ) +
  theme_classic()

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "results/figures/BRCA_PCA.png",
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Volcano plot
res_df$significance <- "Not significant"

res_df$significance[
  res_df$padj < 0.05 &
  res_df$log2FoldChange >= 1
] <- "Upregulated"

res_df$significance[
  res_df$padj < 0.05 &
  res_df$log2FoldChange <= -1
] <- "Downregulated"

res_df$neg_log10_padj <- -log10(res_df$padj)

res_df$neg_log10_padj[
  is.infinite(res_df$neg_log10_padj)
] <- NA

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Gene Expression: TCGA-BRCA",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Category"
  ) +
  theme_classic()

ggsave(
  "results/figures/BRCA_volcano.png",
  volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Differential expression analysis complete.\n")