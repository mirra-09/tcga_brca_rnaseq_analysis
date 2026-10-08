# TCGA-BRCA pathway enrichment analysis
# GO Biological Process and KEGG enrichment of differentially expressed genes

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load annotated differential expression results
res_annotated <- read.csv(
  "data/BRCA_differential_expression_annotated.csv"
)

# Significant genes
sig_genes <- res_annotated[
  res_annotated$padj < 0.05 &
  abs(res_annotated$log2FoldChange) >= 1 &
  !is.na(res_annotated$gene_name),
]

gene_list <- unique(sig_genes$gene_name)

# GO Biological Process enrichment
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(ego),
  "data/BRCA_GO_enrichment.csv",
  row.names = FALSE
)

go_plot <- dotplot(
  ego,
  showCategory = 15
) +
  ggtitle("GO Biological Process Enrichment - TCGA-BRCA") +
  theme_classic()

ggsave(
  "results/figures/BRCA_GO_enrichment.png",
  go_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# KEGG enrichment
kegg_genes <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ekegg <- enrichKEGG(
  gene = kegg_genes$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

write.csv(
  as.data.frame(ekegg),
  "data/BRCA_KEGG_enrichment.csv",
  row.names = FALSE
)

kegg_plot <- dotplot(
  ekegg,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment - TCGA-BRCA") +
  theme_classic()

ggsave(
  "results/figures/BRCA_KEGG_enrichment.png",
  kegg_plot,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Enrichment analysis complete.\n")