

library(dplyr)



# Load candidate genes

network_nodes <- read.csv(
  "results/tables/final_candidate_table.csv"
)



# Define functional categories

network_nodes$category <- "Chloroplast maintenance"


# Chloroplast RNA processing module

network_nodes$category[
  grepl(
    "RNA editing|editing factor|RNA-binding|transcript",
    network_nodes$GENENAME,
    ignore.case = TRUE
  )
] <- "Chloroplast RNA processing"



# Photosynthetic machinery module

network_nodes$category[
  network_nodes$category == "Chloroplast maintenance" &
    grepl(
      "photosystem|photosynthetic|light-harvesting|thylakoid",
      network_nodes$GENENAME,
      ignore.case = TRUE
    )
] <- "Photosynthetic machinery"



# Add CLPD as central stress candidate

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
    "CLPD",
    "Chloroplast RNA processing",
    "Photosynthetic machinery"
  ),
  
  target = c(
    "Clp proteostasis",
    "Chloroplast RNA processing",
    "Photosynthetic machinery",
    "Chloroplast maintenance",
    "Plastid transcript regulation",
    "Light harvesting machinery"
  ),
  
  relationship = c(
    "associated with increased chloroplast proteostasis response",
    "associated with chloroplast RNA maintenance",
    "associated with altered photosynthetic program",
    "linked to chloroplast stress adaptation",
    "increased chloroplast RNA processing activity",
    "reduced photosynthetic capacity"
  )
  
)



# Save network edges

write.csv(
  edges,
  "results/tables/chloroplast_network_edges.csv",
  row.names = FALSE
)



# Summary

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