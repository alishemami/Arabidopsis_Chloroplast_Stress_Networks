# Arabidopsis Chloroplast Stress RNA-seq Analysis

A reproducible computational workflow to investigate chloroplast-associated stress responses in Arabidopsis thaliana using RNA-seq data analysis.

## Overview

Chloroplasts are dynamic organelles that respond to environmental stress through coordinated regulation of photosynthesis, protein homeostasis, RNA processing, and 
organelle maintenance pathways.

This project analyzes RNA-seq data from stressed and control Arabidopsis thaliana samples to identify stress-responsive genes, characterize chloroplast-associated 
transcriptional responses, and generate candidate genes for further investigation.

The workflow combines differential expression analysis, functional annotation, candidate gene identification, and a conceptual chloroplast stress network model.

The analysis was designed as a reproducible pipeline where each step can be independently executed and traced.

---

## Dataset

Organism: Arabidopsis thaliana

Dataset: GSE265858

Input file:

data/GSE265858_RNAseq_counts.txt

The dataset contains gene-level RNA-seq count data used for downstream differential expression analysis.

---

## Analysis workflow

The complete analysis is organized into five main R scripts:

R/

01_preprocess_data.R  
02_DESeq2_analysis.R  
03_functional_analysis.R  
04_candidate_visualization.R  
05_network_analysis.R  

---

## 01. Data preprocessing

01_preprocess_data.R

This step prepares the raw RNA-seq count matrix for downstream analysis.

Main steps:

- Import raw count data
- Define experimental groups
- Prepare count matrix and sample metadata
- Generate inputs required for DESeq2 analysis

---

## 02. Differential expression analysis

02_DESeq2_analysis.R

Differential expression analysis was performed using DESeq2.

Comparison:

Stress vs Control

Significance criteria:

adjusted p-value < 0.05

|log2FoldChange| > 1

Differential expression summary:

- 919 significant differentially expressed genes
- 453 upregulated genes
- 466 downregulated genes

Output:

results/tables/DESeq2_all_results.csv

---

## Functional annotation and candidate identification

03_functional_analysis.R

Functional annotation and Gene Ontology enrichment analysis were performed using Arabidopsis annotation resources.

Chloroplast-associated candidates were identified based on:

- chloroplast-related gene descriptions
- plastid-associated functions
- thylakoid-related processes
- chloroplast-associated GO annotations

The analysis identified:

31 chloroplast-associated candidate genes

Output:

results/tables/final_candidate_table.csv

---

## CLPD integration

The chloroplast proteostasis regulator CLPD (AT5G51070) was incorporated into the final candidate model as a predefined candidate based on its biological relevance.

CLPD showed strong transcriptional induction under stress conditions:

log2FoldChange ≈ +3.56

adjusted p-value ≈ 8.5e-144

CLPD was not selected solely through automated annotation filtering. Instead, it was included as a hypothesis-driven candidate to connect transcriptomic evidence with 
chloroplast proteostasis biology.

---

## Candidate visualization

04_candidate_visualization.R

This step generates visual summaries of the identified candidates.

Generated figures:

results/figures/volcano_plot.png

results/figures/chloroplast_candidate_plot.png

results/figures/chloroplast_candidate_heatmap.png

---

## Chloroplast stress network model

05_network_analysis.R

A conceptual chloroplast stress network was constructed by integrating:

- annotation-derived chloroplast candidates
- CLPD as a central candidate
- functional relationships related to chloroplast maintenance, RNA processing, and photosynthetic responses

Final network:

Network nodes: 32

Network edges: 6

Network tables:

results/tables/chloroplast_network_nodes.csv

results/tables/chloroplast_network_edges.csv

The network represents a computational hypothesis model and requires experimental validation.

---

## Reproducibility

The complete workflow can be reproduced by running:

source("R/01_preprocess_data.R")

source("R/02_DESeq2_analysis.R")

source("R/03_functional_analysis.R")

source("R/04_candidate_visualization.R")

source("R/05_network_analysis.R")

Required R packages:

- DESeq2
- ggplot2
- pheatmap
- gprofiler2
- AnnotationDbi
- org.At.tair.db
- dplyr

---

## Project structure

Arabidopsis_Chloroplast_Stress_Networks/

R/

01_preprocess_data.R

02_DESeq2_analysis.R

03_functional_analysis.R

04_candidate_visualization.R

05_network_analysis.R


data/

GSE265858_RNAseq_counts.txt


results/

figures/

tables/


README.md

---

## Future directions

This project provides a computational framework for identifying chloroplast stress-associated candidates and generating hypotheses for further experimental investigation.

Future directions include:

- experimental validation of candidate genes
- investigation of CLPD-dependent chloroplast proteostasis responses
- integration with additional transcriptomic and proteomic datasets

---

## Project status

This project is currently under development.

All computational predictions and network relationships represent hypotheses generated from transcriptomic analysis and require experimental validation.
