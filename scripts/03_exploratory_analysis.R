library(ggplot2)

load("data/tcga_brca_filtered.RData")

log_counts <- log2(filtered_counts + 1)

pca <- prcomp(t(log_counts), scale.=TRUE)

sample_type <- metadata$shortLetterCode  # adjust if needed

if (!dir.exists("results")) {
  dir.create("results")
}

p <- ggplot(data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Type = sample_type
), aes(x = PC1, y = PC2, color = Type)) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("PCA of TCGA BRCA Samples")

ggsave("results/pca_plot.png", plot = p, width = 8, height = 6)
