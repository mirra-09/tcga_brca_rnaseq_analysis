# TCGA-BRCA RNA-seq data query
# Query all STAR-count RNA-seq samples for tumor and normal tissue

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
