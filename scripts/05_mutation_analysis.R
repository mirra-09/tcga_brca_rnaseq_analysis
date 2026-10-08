# TCGA-BRCA mutation-expression analysis
# Compare gene expression between PIK3CA- and TP53-mutant and wild-type tumors

library(DESeq2)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load processed RNA-seq data
load("data/brca_filtered_counts.RData")
load("data/tcga_brca_100samples.RData")

# Load mutation data
# Mutation MAF files are downloaded from the GDC and are intentionally
# excluded from this repository.
maf_files <- list.files(
  "GDCdata",
  recursive = TRUE,
  pattern = "maf|masked_somatic",
  ignore.case = TRUE,
  full.names = TRUE
)

if (length(maf_files) == 0) {
  stop(
    "No mutation files found. Download the TCGA-BRCA ",
    "Masked Somatic Mutation files from the GDC first."
  )
}

# Read MAF files
maf_list <- lapply(
  maf_files,
  function(x) {
    tryCatch(
      read.delim(
        x,
        comment.char = "#",
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )
  }
)

maf_list <- maf_list[
  !sapply(maf_list, is.null)
]

maf <- do.call(rbind, maf_list)

# Keep non-silent mutations
non_silent <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site"
)

maf <- maf[
  maf$Variant_Classification %in% non_silent,
]

# Patient-level identifiers
maf$patient_id <- substr(
  maf$Tumor_Sample_Barcode,
  1,
  12
)

# RNA-seq tumor samples
tumor_metadata <- as.data.frame(
  colData(brca_se)
)

tumor_metadata <- tumor_metadata[
  tumor_metadata$shortLetterCode == "TP",
]

tumor_metadata$patient_id <- substr(
  rownames(tumor_metadata),
  1,
  12
)

tumor_patient_ids <- tumor_metadata$patient_id

# Restrict mutations to the RNA-seq cohort
maf_cohort <- maf[
  maf$patient_id %in% tumor_patient_ids,
]

# Mutation status function
get_mutation_status <- function(
  gene,
  mutation_data,
  patient_ids
) {

  mutated <- unique(
    mutation_data$patient_id[
      mutation_data$Hugo_Symbol == gene
    ]
  )

  data.frame(
    patient_id = patient_ids,
    status = ifelse(
      patient_ids %in% mutated,
      "Mutant",
      "Wild-type"
    )
  )
}

# PIK3CA status
pik3ca_status <- get_mutation_status(
  "PIK3CA",
  maf_cohort,
  tumor_patient_ids
)

# TP53 status
tp53_status <- get_mutation_status(
  "TP53",
  maf_cohort,
  tumor_patient_ids
)

mutation_status <- data.frame(
  patient_id = tumor_patient_ids,
  PIK3CA_status = pik3ca_status$status,
  TP53_status = tp53_status$status
)

# Save mutation status table
write.csv(
  mutation_status,
  "data/mutation_status_matched_cohort.csv",
  row.names = FALSE
)

# Match mutation status to expression samples
tumor_columns <- rownames(tumor_metadata)

tumor_patient_from_expression <- substr(
  tumor_columns,
  1,
  12
)

matched <- tumor_patient_from_expression %in%
  mutation_status$patient_id

expression_samples <- tumor_columns[matched]

matched_ids <- tumor_patient_from_expression[matched]

mutation_status_matched <- mutation_status[
  match(
    matched_ids,
    mutation_status$patient_id
  ),
]

# Create matched expression matrix
tumor_counts <- counts_filtered[
  ,
  expression_samples
]

tumor_metadata_matched <- tumor_metadata[
  expression_samples,
  ,
  drop = FALSE
]

# PIK3CA analysis
coldata_pik3ca <- tumor_metadata_matched

coldata_pik3ca$mutation_status <-
  factor(
    mutation_status_matched$PIK3CA_status,
    levels = c("Wild-type", "Mutant")
  )

dds_pik3ca <- DESeqDataSetFromMatrix(
  countData = tumor_counts,
  colData = coldata_pik3ca,
  design = ~ mutation_status
)

dds_pik3ca <- DESeq(
  dds_pik3ca
)

res_pik3ca <- results(
  dds_pik3ca,
  contrast = c(
    "mutation_status",
    "Mutant",
    "Wild-type"
  )
)

res_pik3ca <- res_pik3ca[
  order(res_pik3ca$padj),
]

# PIK3CA annotation
gene_annotation <- as.data.frame(
  rowData(brca_se)
)

gene_annotation <- gene_annotation[
  ,
  c("gene_id", "gene_name", "gene_type")
]

res_pik3ca_df <- as.data.frame(
  res_pik3ca
)

res_pik3ca_df$gene_id <-
  rownames(res_pik3ca_df)

res_pik3ca_annotated <- merge(
  res_pik3ca_df,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

write.csv(
  res_pik3ca_annotated,
  "data/PIK3CA_mutant_vs_wildtype_DEG.csv",
  row.names = FALSE
)

# PIK3CA volcano plot
res_pik3ca_annotated$significance <-
  "Not significant"

res_pik3ca_annotated$significance[
  res_pik3ca_annotated$padj < 0.05 &
  res_pik3ca_annotated$log2FoldChange >= 1
] <- "Upregulated"

res_pik3ca_annotated$significance[
  res_pik3ca_annotated$padj < 0.05 &
  res_pik3ca_annotated$log2FoldChange <= -1
] <- "Downregulated"

res_pik3ca_annotated$neg_log10_padj <-
  -log10(res_pik3ca_annotated$padj)

pik3ca_volcano <- ggplot(
  res_pik3ca_annotated,
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
    title =
      "Differential Expression by PIK3CA Mutation Status",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Category"
  ) +
  theme_classic()

ggsave(
  "results/figures/PIK3CA_mutant_vs_wildtype_volcano.png",
  pik3ca_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# TP53 analysis
coldata_tp53 <- tumor_metadata_matched

coldata_tp53$mutation_status <-
  factor(
    mutation_status_matched$TP53_status,
    levels = c("Wild-type", "Mutant")
  )

dds_tp53 <- DESeqDataSetFromMatrix(
  countData = tumor_counts,
  colData = coldata_tp53,
  design = ~ mutation_status
)

dds_tp53 <- DESeq(
  dds_tp53
)

res_tp53 <- results(
  dds_tp53,
  contrast = c(
    "mutation_status",
    "Mutant",
    "Wild-type"
  )
)

res_tp53 <- res_tp53[
  order(res_tp53$padj),
]

# TP53 annotation
res_tp53_df <- as.data.frame(
  res_tp53
)

res_tp53_df$gene_id <-
  rownames(res_tp53_df)

res_tp53_annotated <- merge(
  res_tp53_df,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

write.csv(
  res_tp53_annotated,
  "data/TP53_mutant_vs_wildtype_DEG.csv",
  row.names = FALSE
)

# TP53 volcano plot
res_tp53_annotated$significance <-
  "Not significant"

res_tp53_annotated$significance[
  res_tp53_annotated$padj < 0.05 &
  res_tp53_annotated$log2FoldChange >= 1
] <- "Upregulated"

res_tp53_annotated$significance[
  res_tp53_annotated$padj < 0.05 &
  res_tp53_annotated$log2FoldChange <= -1
] <- "Downregulated"

res_tp53_annotated$neg_log10_padj <-
  -log10(res_tp53_annotated$padj)

tp53_volcano <- ggplot(
  res_tp53_annotated,
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
    title =
      "Differential Expression by TP53 Mutation Status",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Category"
  ) +
  theme_classic()

ggsave(
  "results/figures/TP53_mutant_vs_wildtype_volcano.png",
  tp53_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Mutation-expression analysis complete.\n")
cat(
  "Matched RNA/mutation samples:",
  length(expression_samples),
  "\n"
)