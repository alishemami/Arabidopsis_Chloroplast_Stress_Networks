# Prepare RNA-seq data for downstream analysis

library(DESeq2)
library(ggplot2)
library(pheatmap)


# Import raw RNA-seq count matrix

rna_data <- read.table(
  "data/GSE265858_RNAseq_counts.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# Inspect imported data

head(rna_data)

dim(rna_data)

colnames(rna_data)



# Extract count data from experimental samples

counts <- rna_data[, c(
  "WT-1_count",
  "WT-2_count",
  "#5-1_count",
  "#5-2_count"
)]


# Rename samples according to experimental conditions

colnames(counts) <- c(
  "WT_1",
  "WT_2",
  "Stress_1",
  "Stress_2"
)


# Store Arabidopsis gene identifiers as row names

rownames(counts) <- rna_data$AGI


head(counts)

dim(counts)



# Define sample metadata

sample_info <- data.frame(
  condition = c(
    "control",
    "control",
    "stress",
    "stress"
  )
)


rownames(sample_info) <- colnames(counts)

sample_info$condition <- factor(sample_info$condition)


sample_info



# Save processed data for downstream analysis

save(
  counts,
  sample_info,
  file = "data/processed_counts.RData"
)



# Create DESeq2 object for quality control

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)


# Normalize counts for sample-level visualization

dds <- estimateSizeFactors(dds)


vsd <- vst(
  dds,
  blind = TRUE
)



# PCA analysis

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)


percent_var <- round(
  100 * attr(pca_data, "percentVar")
)


pca_plot <- ggplot(
  pca_data,
  aes(
    PC1,
    PC2,
    color = condition,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_classic()


ggsave(
  "results/figures/PCA_plot.png",
  pca_plot,
  width = 7,
  height = 5
)



# Sample correlation heatmap

sample_cor <- cor(
  assay(vsd),
  method = "pearson"
)


pheatmap(
  sample_cor,
  main = "Sample correlation",
  filename = "results/figures/sample_correlation_heatmap.png"
)