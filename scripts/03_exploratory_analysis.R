library(ggplot2)

# Load filtered data
load("data/tcga_brca_filtered.RData")

# Log transform
log_counts <- log2(filtered_counts + 1)

# PCA
pca <- prcomp(t(log_counts), scale.=TRUE)

# Extract sample type (adjust column name if needed)
sample_type <- metadata$sample_type

# Create results folder
if (!dir.exists("results")) dir.create("results")

# PCA plot
p <- ggplot(data.frame(PC1=pca$x[,1],
                       PC2=pca$x[,2],
                       Type=sample_type),
            aes(x=PC1, y=PC2, color=Type)) +
  geom_point(size=2) +
  theme_minimal() +
  ggtitle("PCA of TCGA BRCA Samples")

ggsave("results/pca_plot.png", plot=p, width=8, height=6)

cat("Exploratory analysis complete.\n")
