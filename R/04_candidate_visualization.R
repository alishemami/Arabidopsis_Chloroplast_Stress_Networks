# Visualization of differential expression and chloroplast candidates


library(ggplot2)
library(pheatmap)
library(dplyr)



# Load DESeq2 results

deg_results <- read.csv(
  "results/tables/DESeq2_all_results.csv",
  row.names = 1
)



# Prepare DEG classification for volcano plot

deg_results <- deg_results[
  !is.na(deg_results$padj),
]


deg_results$DEG <- "Not significant"


deg_results$DEG[
  deg_results$padj < 0.05 &
    deg_results$log2FoldChange > 1
] <- "Up"


deg_results$DEG[
  deg_results$padj < 0.05 &
    deg_results$log2FoldChange < -1
] <- "Down"



# Volcano plot

volcano_plot <- ggplot(
  deg_results,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = DEG
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  theme_classic() +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    title = "Differential expression under chloroplast stress"
  )


ggsave(
  "results/figures/volcano_plot.png",
  volcano_plot,
  width = 7,
  height = 5
)



# Load chloroplast candidate table

final_candidates <- read.csv(
  "results/tables/final_candidate_table.csv"
)



# Candidate expression plot

final_candidates$label <- ifelse(
  is.na(final_candidates$SYMBOL),
  final_candidates$TAIR,
  final_candidates$SYMBOL
)


candidate_plot <- ggplot(
  final_candidates,
  aes(
    x = log2FoldChange,
    y = reorder(label, log2FoldChange),
    color = direction
  )
) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "log2 Fold Change",
    y = "Candidate genes",
    title = "Chloroplast stress candidate genes"
  )


ggsave(
  "results/figures/chloroplast_candidate_plot.png",
  candidate_plot,
  width = 8,
  height = 7
)



# Candidate heatmap

load(
  "data/processed_counts.RData"
)


candidate_ids <- final_candidates$TAIR


candidate_counts <- counts[
  rownames(counts) %in% candidate_ids,
]


candidate_expression <- log2(
  candidate_counts + 1
)


candidate_expression <- t(
  scale(t(candidate_expression))
)


pheatmap(
  candidate_expression,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = "results/figures/chloroplast_candidate_heatmap.png",
  main = "Chloroplast stress candidate genes"
)