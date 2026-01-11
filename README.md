# ProtAN (**Prot**eome**AN**alyzer)

ProtAN is a bioinformatics tool dedicated for analysis of large proteomics datasets. It identifies the differentially expressed proteins (DEPs) mapping them onto molecular pathways thus allowing comprehensive interpretation and analysis of biologically relevant data.

## Introduction

Mass spectrometry-based Proteomics is the large-scale study of proteins to understand biological systems.
The output data from a mass spectrometer is analyzed in two parts: pre- and post-processing.

In pre-processing, the proteins in the samples are identified by aligning them against a known database followed by data normalization, missing value imputation and log2 transformation. Software such as MaxQuant (DDA-based) or Perseus (DIA-based) can be used for this.

The post-processing steps include statistically identifying the dysregulated proteins and their functional role in biological networks. Downstream analysis also includes protein-protein interactions, motif enrichment, peptide selection and many more.

## Motivation 

Working in the Proteomics field I realized that single, unified tools are available to perform data pre-processing (eg. MaxQuant/ Perseus) but very few to perform robust post-processing data analysis. 

Researchers have to browse through multiple tools and software, each having their compatibility issues, different file formats and handling multiple files. This makes it frustrating for the researchers.

Thus, ProtAN is a combined tool dedicated for analyzing post-processed proteomics data with a single file input.

# ProtAN Pipeline Overview

Load data (.csv file)
           ↓

Fold-change analysis (Statistical analysis)
         
          ↓

Volcano plot

          ↓

Select DEPs

          ↓

PCA (DEPs)

          ↓
Heatmap (DEPs)
 ↓
Pathway enrichment (Reactome + KEGG)
 ↓
PPI network (STRING)
 ↓
Top-hit prioritization
 ↓
Motif enrichment analysis
 ↓
PDF report generation



## 1. Input file

 `.csv` file in the following format-

Values are log2-transformed, normalized, imputed from the pre-processing steps

## 2. Statistical analysis

Includes Student t-test, Fold change= 1.5, p-value < 0.05, Mann–Whitney

Volcano plot, selection of top differentially regulated proteins, PCA analysis and heatmap

## 3. Pathway Enrichment

Top DEPs are mapped onto their pathways and sub-pathways using Reactome and KEGG. 
Usage of `.json` files, `API calls` and `requests`

## 4.Protein-Protein Interaction (PPI) Network Building

Build interaction networks for deregulated proteins
Identify hubs and bottlenecks

Usage of `networkx` and STRING API

## 5. Motif Enrichment

Use of Uniprot API, Biopython

## 6. Output files

results/

├── volcano.png

├── pca.png

├── heatmap.png

├── reactome_enrichment.csv

├── kegg_enrichment.csv

├── ppi_top_hits.csv

├── motif_enrichment.csv

├── proteomics_analysis_report.pdf

└── motif_and_network_report.pdf



## Installations

Python v3.9+
pip install pandas numpy scipy statsmodels matplotlib seaborn scikit-learn
pip install gseapy networkx requests biopython reportlab


## Dataset used

Investigating Cisplatin Resistance in Squamous Cervical Cancer: Proteomic Insights into DNA Repair Pathways and Omics-Based Drug Repurposing
Amrita Mukherjee, Sayan Manna, Avinash Singh, Adrija Ray, and Sanjeeva Srivastava
Journal of Proteome Research 2025 24 (6), 2628-2642
DOI: 10.1021/acs.jproteome.4c00885

The dataset used here is log2 transformed and normalized Table S4.

The format to be used is-


## NOTE

This code works for paired datasets eg-
Cisplatin-sensitive vs resistant TNBC

Drug A vs control

Hypoxia vs normoxia

Knockdown vs WT

This code does not generalize-
Paired vs unpaired design
Column naming inconsistency
Gene/Protein ID inconsistency
