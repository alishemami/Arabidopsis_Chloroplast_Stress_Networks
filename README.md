# Arabidopsis Chloroplast Stress Response: RNA-seq Analysis Pipeline

A reproducible computational workflow to investigate chloroplast-associated stress responses in *Arabidopsis thaliana* using RNA-seq data analysis.

## Overview

Chloroplasts are dynamic organelles that continuously adapt their activity in response to environmental changes. These adaptations involve coordinated regulation of photosynthesis, chloroplast protein homeostasis, RNA processing, and organelle maintenance pathways.

This project analyzes RNA-seq data from control and stress-associated *Arabidopsis thaliana* samples to identify differentially expressed genes, characterize chloroplast-related responses, and prioritize candidate genes for further investigation.

The workflow includes:

- RNA-seq preprocessing
- Differential expression analysis using DESeq2
- Functional annotation and GO enrichment analysis
- Chloroplast candidate prioritization
- Candidate visualization
- Construction of a conceptual chloroplast stress response model

---

## Dataset

**Organism:** *Arabidopsis thaliana*

**Dataset:** GSE265858

**Input file:**

```
data/GSE265858_RNAseq_counts.txt
```

The input consists of gene-level RNA-seq count data used for downstream differential expression analysis.

The original dataset sample labels were simplified into control and stress-associated conditions for downstream analysis.

---

## Analysis workflow

```
R/
├── 01_preprocess_data.R
├── 02_DESeq2_analysis.R
├── 03_functional_analysis.R
├── 04_candidate_visualization.R
└── 05_network_analysis.R
```

---

## Differential expression analysis

Differential expression analysis was performed using **DESeq2**.

Comparison:

```
stress vs control
```

Significance criteria:

```
adjusted p-value < 0.05
|log2FoldChange| > 1
```

The analysis identified:

- 919 significant differentially expressed genes
- 453 upregulated genes
- 466 downregulated genes

Example candidate:

```
AT5G51070 (CLPD)

log2FoldChange: 3.56
adjusted p-value: 8.48e-144
```

---

## Functional analysis and candidate prioritization

Functional annotation and GO enrichment analysis were performed to characterize chloroplast-associated responses.

Candidate prioritization integrates:

- Differential expression strength
- Chloroplast-associated annotation
- Biological relevance to stress adaptation

GO enrichment was used as an independent functional characterization step.

---

## Conceptual chloroplast stress network

The final model summarizes relationships between major chloroplast stress response modules:

- Clp proteostasis
- Chloroplast RNA processing
- Photosynthetic machinery
- Chloroplast maintenance

The network represents a conceptual biological model and is not a measured protein-protein interaction or co-expression network.

---

## Outputs

Main generated outputs include:

Differential expression results:

```
results/tables/DESeq2_all_results.csv
```

Candidate gene table:

```
results/tables/final_candidate_table.csv
```

Network tables:

```
results/tables/chloroplast_network_nodes.csv
results/tables/chloroplast_network_edges.csv
```

Figures:

```
results/figures/
```

---

## Experimental limitation

The dataset contains two biological replicates per condition.

Therefore, identified differentially expressed genes and candidate genes should be interpreted as hypothesis-generating findings requiring further experimental validation.

---

## Reproducibility

The complete workflow can be executed sequentially:

```r
source("R/01_preprocess_data.R")
source("R/02_DESeq2_analysis.R")
source("R/03_functional_analysis.R")
source("R/04_candidate_visualization.R")
source("R/05_network_analysis.R")
```

---

## Software environment

Main R packages:

- DESeq2
- ggplot2
- pheatmap
- gprofiler2
- dplyr

Session information is provided in:

```
sessionInfo.txt
```

---

## License

This project is released under the MIT License.
