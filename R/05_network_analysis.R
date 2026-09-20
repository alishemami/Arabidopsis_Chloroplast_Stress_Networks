# Build chloroplast stress network tables

library(dplyr)


# Load candidate genes

network_nodes <- read.csv(
  "results/tables/final_candidate_table.csv"
)


# Define functional categories
# Order matters: specific categories first

network_nodes$category <- "Chloroplast maintenance"


network_nodes$category[
  grepl(
    "RNA editing|editing factor|RNA-binding domain|transcript",
    network_nodes$GENENAME,
    ignore.case = TRUE
  )
] <- "Chloroplast RNA processing"


network_nodes$category[
  grepl(
    "photosystem|photosynthetic|light-harvesting|thylakoid",
    network_nodes$GENENAME,
    ignore.case = TRUE
  )
] <- "Photosynthetic machinery"



# Add CLPD as predefined central candidate

clpd_node <- data.frame(
  
  TAIR = "AT5G51070",
  SYMBOL = "CLPD",
  GENENAME = "ATP-dependent Clp protease regulatory subunit",
  log2FoldChange = 3.564841,
  padj = 8.478012e-144,
  direction = "Up",
  category = "Clp proteostasis",
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



# Save nodes

write.csv(
  network_nodes,
  "results/tables/chloroplast_network_nodes.csv",
  row.names = FALSE
)



# Conceptual biological relationships

edges <- data.frame(
  
  source = c(
    "CLPD",
    "CLPD",
    "CLPD",
    "Chloroplast RNA processing",
    "Photosynthetic machinery"
  ),
  
  target = c(
    "Clp proteostasis",
    "Chloroplast RNA processing",
    "Photosynthetic machinery",
    "Plastid transcript regulation",
    "Light harvesting machinery"
  ),
  
  relationship = c(
    "associated with increased chloroplast proteostasis response",
    "associated with chloroplast RNA maintenance",
    "associated with altered photosynthetic program",
    "increased chloroplast RNA processing activity",
    "reduced photosynthetic capacity"
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