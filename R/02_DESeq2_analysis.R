# ============================================
# Step 03: Differential expression analysis
# Arabidopsis Chloroplast Stress Networks
# ============================================


# Load libraries

library(DESeq2)
library(ggplot2)


# Load prepared data

load("data/processed_counts.RData")


# Create DESeq2 object

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)


# Run differential expression analysis

dds <- DESeq(dds)


# Extract results

results_deseq <- results(
  dds,
  contrast = c(
    "condition",
    "stress",
    "control"
  )
)


# Order results by adjusted p-value

results_deseq <- results_deseq[
  order(results_deseq$padj),
]


# Show top genes

head(results_deseq)

write.csv(
  as.data.frame(results_deseq),
  "results/tables/DESeq2_all_results.csv"
)
