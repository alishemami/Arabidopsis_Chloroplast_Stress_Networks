# Functional analysis and chloroplast candidate identification


library(gprofiler2)
library(AnnotationDbi)
library(org.At.tair.db)
library(dplyr)



# Load differential expression results

deg_results <- read.csv(
  "results/tables/DESeq2_all_results.csv",
  row.names = 1
)



# Remove genes without adjusted p-values

deg_results <- deg_results[
  !is.na(deg_results$padj),
]



# Define significant differentially expressed genes

DEG_list <- rownames(
  deg_results[
    deg_results$padj < 0.05 &
      abs(deg_results$log2FoldChange) > 1,
  ]
)



# Separate up- and downregulated genes

up_list <- rownames(
  deg_results[
    deg_results$padj < 0.05 &
      deg_results$log2FoldChange > 1,
  ]
)


down_list <- rownames(
  deg_results[
    deg_results$padj < 0.05 &
      deg_results$log2FoldChange < -1,
  ]
)



# Perform GO enrichment analysis

GO_all <- gost(
  query = DEG_list,
  organism = "athaliana",
  significant = TRUE
)


GO_up <- gost(
  query = up_list,
  organism = "athaliana",
  significant = TRUE
)


GO_down <- gost(
  query = down_list,
  organism = "athaliana",
  significant = TRUE
)


# Save GO enrichment results

GO_up_table <- GO_up$result

GO_down_table <- GO_down$result


# Keep only simple columns that can be exported as CSV

GO_up_table <- GO_up_table[
  sapply(GO_up_table, function(x) !is.list(x))
]


GO_down_table <- GO_down_table[
  sapply(GO_down_table, function(x) !is.list(x))
]


write.csv(
  GO_up_table,
  "results/tables/GO_up_results.csv",
  row.names = FALSE
)


write.csv(
  GO_down_table,
  "results/tables/GO_down_results.csv",
  row.names = FALSE
)



# Annotate significant genes using Arabidopsis annotation database

annotate_genes <- function(gene_list) {
  
  annotation <- AnnotationDbi::select(
    org.At.tair.db,
    keys = gene_list,
    columns = c(
      "SYMBOL",
      "GENENAME",
      "GO"
    ),
    keytype = "TAIR"
  )
  
  
  annotation %>%
    group_by(
      TAIR,
      SYMBOL,
      GENENAME
    ) %>%
    summarise(
      GO_terms = paste(
        unique(GO),
        collapse = "; "
      ),
      .groups = "drop"
    )
}



DEG_annotation <- annotate_genes(DEG_list)



# Identify genes associated with chloroplast functions

chloroplast_candidates <- DEG_annotation[
  grepl(
    "chloroplast|plastid|thylakoid|plastidial|stroma",
    paste(
      DEG_annotation$GENENAME,
      DEG_annotation$GO_terms
    ),
    ignore.case = TRUE
  ),
]



# Add differential expression information

final_candidates <- merge(
  chloroplast_candidates,
  as.data.frame(deg_results),
  by.x = "TAIR",
  by.y = "row.names"
)



# Keep one annotation entry per gene

final_candidates <- final_candidates %>%
  distinct(
    TAIR,
    .keep_all = TRUE
  )



# Add regulation direction

final_candidates$direction <- ifelse(
  final_candidates$log2FoldChange > 0,
  "Up",
  "Down"
)



# Save final chloroplast candidate list

write.csv(
  final_candidates,
  "results/tables/final_candidate_table.csv",
  row.names = FALSE
)



cat(
  "Chloroplast candidates:",
  nrow(final_candidates),
  "\n"
)