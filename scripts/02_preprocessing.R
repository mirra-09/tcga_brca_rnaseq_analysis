library(SummarizedExperiment)
library(dplyr)

# Load queried object (change filename if needed)
load("data/tcga_brca_raw.RData")  # replace with your saved object name

# Assume object name is 'se'
# Extract raw counts
counts <- assay(se)

# Extract metadata
metadata <- colData(se)

# Filter low-count genes
keep <- rowSums(counts >= 10) >= 10
filtered_counts <- counts[keep, ]

# Save filtered data
if (!dir.exists("data")) dir.create("data")

save(filtered_counts, metadata,
     file = "data/tcga_brca_filtered.RData")

cat("Preprocessing complete.\n")
