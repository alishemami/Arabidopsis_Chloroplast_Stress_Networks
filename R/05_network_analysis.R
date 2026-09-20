# Build chloroplast stress network tables


library(dplyr)


# Load candidate genes

network_nodes <- read.csv(
  "results/tables/final_candidate_table.csv"
)



# Add CLPD as predefined central candidate

clpd_node <- data.frame(
  TAIR = "AT5G51070",
  SYMBOL = "CLPD",
  GENENAME = "ATP-dependent Clp protease regulatory subunit",
  log2FoldChange = 3.564841,
  padj = 8.478012e-144,
  direction = "Up",
  category = "Clp protease system",
  evidence = "Strong induction of CLPD"
)


network_nodes <- bind_rows(
  network_nodes,
  clpd_node
) %>%
  distinct(
    TAIR,
    .keep_all = TRUE
  )



# Define functional categories

network_nodes$category[
  grepl(
    "RNA editing|editing factor",
    network_nodes$GENENAME,
    ignore.case = TRUE
  )
] <- "Chloroplast RNA editing"


network_nodes$category[
  grepl(
    "thylakoid|plastid|chloroplast",
    network_nodes$GENENAME,
    ignore.case = TRUE
  )
] <- "Chloroplast maintenance"



# Save network nodes

write.csv(
  network_nodes,
  "results/tables/chloroplast_network_nodes.csv",
  row.names = FALSE
)



# Define conceptual biological relationships

edges <- data.frame(
  source = c(
    "CLPD",
    "CLPD",
    "CLPD",
    "RNA editing",
    "Photosynthesis",
    "Chloroplast transcription"
  ),
  
  target = c(
    "RNA editing",
    "Photosynthesis",
    "Chloroplast maintenance",
    "Plastid transcript regulation",
    "Light harvesting machinery",
    "Photosynthetic regulation"
  ),
  
  relationship = c(
    "associated with increased RNA maintenance response",
    "associated with reduced photosynthetic program",
    "linked to chloroplast stress adaptation",
    "increased chloroplast RNA processing activity",
    "reduced photosynthetic capacity",
    "reduced transcriptional regulation of chloroplast programs"
  )
)



write.csv(
  edges,
  "results/tables/chloroplast_network_edges.csv",
  row.names = FALSE
)



cat(
  "Network nodes:",
  nrow(network_nodes),
  "\n"
)

cat(
  "Network edges:",
  nrow(edges),
  "\n"
)