# TCGA-BRCA RNA-seq preprocessing
# Prepare a balanced cohort of 50 primary tumors and 50 solid tissue normals

library(TCGAbiolinks)
library(SummarizedExperiment)

set.seed(42)

# Query TCGA-BRCA STAR-count RNA-seq data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

results <- getResults(query)

# Select a reproducible balanced cohort
tumor_results <- results[results$sample_type == "Primary Tumor", ]
normal_results <- results[results$sample_type == "Solid Tissue Normal", ]

tumor_selected <- tumor_results[
  sample(seq_len(nrow(tumor_results)), size = 50),
]

normal_selected <- normal_results[
  sample(seq_len(nrow(normal_results)), size = 50),
]

selected_results <- rbind(tumor_selected, normal_selected)

# Download selected samples
query_selected <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  barcode = selected_results$cases
)

GDCdownload(query_selected)

# Prepare SummarizedExperiment
brca_se <- GDCprepare(query_selected)

# Extract raw unstranded counts
counts <- assay(brca_se, "unstranded")

# Filter low-expression genes
keep <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep, ]

# Create sample metadata
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

# Save processed objects locally
dir.create("data", showWarnings = FALSE)

save(
  counts_filtered,
  file = "data/brca_filtered_counts.RData"
)

save(
  brca_se,
  file = "data/tcga_brca_100samples.RData"
)

write.csv(
  sample_metadata,
  "data/sample_metadata_100samples.csv",
  row.names = TRUE
)

cat("Preprocessing complete.\n")
cat("Samples:", ncol(counts_filtered), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
