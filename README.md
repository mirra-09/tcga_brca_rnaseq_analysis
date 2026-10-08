# TCGA-BRCA RNA-seq and Cancer Genomics Analysis

An end-to-end cancer genomics analysis of TCGA Breast Invasive Carcinoma (TCGA-BRCA) RNA-seq data, integrating differential gene expression, pathway enrichment, and recurrent somatic mutation analysis.

## Overview

This project investigates transcriptional differences between breast tumor and normal tissue and explores how recurrent somatic mutations are associated with distinct gene-expression programs.

The workflow connects:

RNA-seq data → differential expression → pathway analysis → mutation-expression analysis

A reproducible balanced cohort of 100 TCGA-BRCA samples was analyzed:

- 50 Primary Tumor samples
- 50 Solid Tissue Normal samples

RNA-seq data were obtained from the NCI Genomic Data Commons using TCGAbiolinks and the STAR - Counts workflow.

## Objectives

- Process and analyze TCGA-BRCA RNA-seq data.
- Identify genes differentially expressed between tumor and normal tissue.
- Characterize enriched biological pathways and processes.
- Integrate somatic mutation data with tumor RNA-seq profiles.
- Compare transcriptional programs between PIK3CA-mutant and wild-type tumors.
- Compare transcriptional programs between TP53-mutant and wild-type tumors.

## Analysis Workflow

TCGA-BRCA RNA-seq
        |
        v
100-sample balanced cohort
50 Tumor + 50 Normal
        |
        v
Gene filtering
        |
        v
DESeq2 differential expression
        |
        +------------------+
        |                  |
        v                  v
       PCA          Tumor vs Normal
                           |
                           v
                  GO + KEGG enrichment

Tumor RNA-seq
        |
        +------------------+
        |                  |
        v                  v
    PIK3CA              TP53
 Mutation status      Mutation status
        |                  |
        v                  v
 Mutant vs WT         Mutant vs WT
        |                  |
        v                  v
       GSEA               GSEA

## Data and Preprocessing

TCGA-BRCA RNA-seq data were queried using:

- Project: TCGA-BRCA
- Data category: Transcriptome Profiling
- Data type: Gene Expression Quantification
- Workflow: STAR - Counts
- Sample types: Primary Tumor and Solid Tissue Normal

A fixed random seed (set.seed(42)) was used to select 50 tumor and 50 normal samples.

Genes were filtered using the criterion:

At least 10 counts in at least 10 samples

After filtering, 26,344 genes were retained for downstream analysis.

## Differential Expression

Differential expression between tumor and normal samples was performed using DESeq2.

The experimental design was:

~ condition

with Normal as the reference group.

Using:

- adjusted p-value < 0.05
- absolute log2 fold change ≥ 1

the analysis identified 6,718 differentially expressed genes.

### PCA

Variance-stabilizing transformation followed by principal component analysis showed clear separation between tumor and normal samples.

- PC1 explained 31% of the variance.
- PC2 explained 22% of the variance.
- Together, PC1 and PC2 explained approximately 53% of the variance.

![PCA of TCGA-BRCA samples](results/figures/BRCA_PCA.png)

### Differential Expression Volcano Plot

![Differential expression volcano plot](results/figures/BRCA_volcano.png)

## Pathway Enrichment

Differentially expressed genes were analyzed using:

- Gene Ontology Biological Process enrichment
- KEGG pathway enrichment

Major biological themes included extracellular matrix organization, extracellular structure organization, immune-related processes, cell adhesion, cytokine signaling, calcium signaling, and cell-cycle-associated pathways.

### GO Enrichment

![GO enrichment](results/figures/BRCA_GO_enrichment.png)

### KEGG Enrichment

![KEGG enrichment](results/figures/BRCA_KEGG_enrichment.png)

## Somatic Mutation Analysis

Somatic mutation data were obtained from the NCI Genomic Data Commons and matched to tumor RNA-seq samples using patient-level identifiers.

Of the 50 tumor samples, 46 had matching mutation records and were used for mutation-expression analysis.

Two recurrent breast cancer genes were investigated:

- PIK3CA
- TP53

Within the original 50-tumor cohort:

| Gene | Mutant | Wild-type |
|------|--------|-----------|
| PIK3CA | 17 | 33 |
| TP53 | 15 | 35 |

## PIK3CA Mutation-Associated Expression

PIK3CA-mutant and wild-type tumors were compared using DESeq2.

GO-based GSEA identified enrichment of immune-related programs including:

- Antigen processing and presentation
- B-cell-mediated immunity
- Lymphocyte-mediated immunity
- T-cell activation
- Leukocyte-mediated immune processes

These findings represent exploratory associations within the selected cohort and do not establish causality.

![PIK3CA differential expression](results/figures/PIK3CA_mutant_vs_wildtype_volcano.png)

![PIK3CA GO GSEA](results/figures/PIK3CA_GO_GSEA.png)

## TP53 Mutation-Associated Expression

TP53-mutant and wild-type tumors showed broader transcriptional differences.

GO-based GSEA highlighted programs associated with:

- Cell-cycle progression
- Chromosome segregation
- Nuclear division
- Organelle fission
- Epithelial differentiation
- Keratinization

Some enriched GO categories are annotated as "meiotic" because they contain genes involved in general chromosome-segregation and cell-division machinery. These terms do not imply that breast tumors are undergoing meiosis.

![TP53 differential expression](results/figures/TP53_mutant_vs_wildtype_volcano.png)

![TP53 GO GSEA](results/figures/TP53_GO_GSEA.png)

## Key Findings

1. Tumor and normal breast tissue showed strong transcriptional separation in PCA.
2. 6,718 genes met the differential-expression thresholds.
3. Differentially expressed genes were enriched in extracellular matrix, immune, cell-adhesion, signaling, and cell-cycle-related pathways.
4. PIK3CA mutation status was associated with differences in immune-related transcriptional programs.
5. TP53 mutation status was associated with cell-cycle, chromosome-segregation, and epithelial-differentiation programs.
6. The mutation-expression analysis demonstrates how somatic genomic alterations can be connected to transcriptional phenotypes.

## Repository Structure

tcga_brca_rnaseq_analysis/
│
├── README.md
├── .gitignore
│
├── scripts/
│   ├── 01_tcga_query.R
│   ├── 02_preprocessing.R
│   ├── 03_differential_expression.R
│   ├── 04_enrichment.R
│   └── 05_mutation_analysis.R
│
├── data/
│   ├── BRCA_differential_expression.csv
│   ├── BRCA_differential_expression_annotated.csv
│   ├── BRCA_GO_enrichment.csv
│   ├── BRCA_KEGG_enrichment.csv
│   ├── PIK3CA_mutant_vs_wildtype_DEG.csv
│   ├── TP53_mutant_vs_wildtype_DEG.csv
│   └── TP53_GO_GSEA.csv
│
├── results/
│   └── figures/
│       ├── BRCA_PCA.png
│       ├── BRCA_volcano.png
│       ├── BRCA_GO_enrichment.png
│       ├── BRCA_KEGG_enrichment.png
│       ├── PIK3CA_mutant_vs_wildtype_volcano.png
│       ├── PIK3CA_GO_GSEA.png
│       ├── TP53_mutant_vs_wildtype_volcano.png
│       └── TP53_GO_GSEA.png
│
└── environment/
    └── sessionInfo.txt

## Tools and Packages

The analysis was performed using R and Bioconductor packages including:

- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- ggplot2

Session information is provided in:

environment/sessionInfo.txt

## Limitations

- This analysis uses a balanced 100-sample subset rather than the complete TCGA-BRCA cohort.
- The cohort was selected for reproducibility and computational practicality.
- Mutation-expression analyses were limited to 46 tumor samples with matching mutation records.
- Clinical follow-up information retrieved for this analysis was incomplete, so survival analysis was not performed.
- Mutation-associated expression analyses are exploratory and do not establish causality.
- Differentially expressed genes should not automatically be interpreted as cancer driver genes.
- Raw GDC downloads and large serialized R objects are excluded from the repository.

## Reproducibility

The analysis scripts document the workflow used to query, process, analyze, and visualize the TCGA-BRCA data.

A fixed random seed was used for cohort selection to make the 100-sample analysis reproducible.

Raw TCGA/GDC files are not included in this repository because of their size. They can be retrieved directly from the NCI Genomic Data Commons using the workflow documented in the scripts.

Software and package information used for the analysis is provided in environment/sessionInfo.txt.