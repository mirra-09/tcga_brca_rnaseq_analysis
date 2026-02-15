library(dplyr)

# Load data
counts <- read.csv("data/raw_counts_matrix.csv", row.names = 1)
metadata <- read.csv("data/sample_metadata.csv")

# Filter low-expression genes
keep <- rowSums(counts >= 10) >= 10
filtered_counts <- counts[keep, ]

# Save filtered data
save(filtered_counts, metadata,
     file = "data/tcga_brca_filtered.RData")

cat("Preprocessing complete.\n")
