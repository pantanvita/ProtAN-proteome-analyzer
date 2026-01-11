# ProtAN (Proteome ANalyzer)

ProtAN is a Python-based bioinformatics pipeline dedicated for the comparative analysis of large-scale proteomics datasets. It analyzes the differentially expressed proteins (DEPs) mapping them onto molecular pathways thus allowing comprehensive interpretation and analysis of biologically relevant data.

## 📊 Introduction

Drug sensitivity and resistance in cancer research remains a major clinical challenge, necessitating systematic approaches to identify deregulated proteins, affected pathways, and altered protein–protein interaction (PPI) networks.

High-throughput mass spectrometry-based Proteomics is widely used to investigate proteins and molecular mechanisms underlying drug sensitivity and resistance in cancer.

The output data from a mass spectrometer is analyzed in two parts: pre- and post-processing.

In pre-processing, the proteins in the samples are identified by aligning against a known database followed by data normalization, missing value imputation and log2- transformation. Software such as MaxQuant (DDA-based) or Perseus (DIA-based) can be used for this.

The post-processing steps include statistically identifying the dysregulated proteins and their functional role in biological networks. Downstream analysis includes protein-protein interactions, motif enrichment, peptide selection and many more.

## 💡 Motivation for the project

While working on the analysis of drug resistivity proteomics data I realized that single, unified tools are available for data pre-processing (eg. MaxQuant/ Perseus) but very few to perform robust post-processing data analysis. 

Researchers have to browse through multiple tools and software, each having their compatibility issues, different file formats and handling multiple files. This makes it cumbersome for the researchers.

Thus, ProtAN provides a reproducible, command-line–driven Python pipeline for comparative analysis of drug-sensitive vs resistant proteomics data.
The workflow assumes that proteomics data have already undergone missing value imputation, normalization, and log-transformation, allowing the pipeline to focus on downstream statistical analysis, functional interpretation, and network biology.

## 🔎 ProtAN Pipeline Overview

The pipeline performs the following major steps:

1. Input data parsing and validation
2. Differential expression analysis (fold change and statistics)
3. Data visualization (volcano plot, heatmap, PCA)
4. Pathway enrichment analysis (Reactome and KEGG)
5. Protein–protein interaction (PPI) network construction using STRING
6. Network-level statistics and hub identification
7. Motif enrichment analysis among deregulated proteins
8. Automated PDF report generation

All steps are executed through a single CLI tool, ensuring reproducibility and scalability to other datasets with similar structure.

## 📥 Input data format

 `.csv` or `.xlsx` file in the following format-
 
 ![image](prot-data-str.png)

S1, S2, S3: drug-sensitive biological replicates

R1, R2, R3: drug-resistant biological replicates

Assumptions-
* Values are log-transformed (log2 intensity)
* Data are normalized across samples
* Missing values have been imputed prior to analysis

## 📝 Methods

**1. Data Parsing and Quality Control**

Used of pandas to load csv or Excel files. Column names are stripped of whitespace and validated to ensure correct group assignment. Numerical columns are coerced into floating-point format to prevent downstream statistical errors.

**2. Differential Expression Analysis**

**Fold Change Calculation**

For each protein:
* Mean expression is calculated across biological replicates for each condition
* Log2 fold change is computed as:

  log2FC = mean(R) − mean(S)

Positive values indicate upregulation in resistant cells, while negative values indicate downregulation.

**Statistical Testing**

* An unpaired two-sample t-test is applied between sensitive and resistant replicates for each protein. P-values are adjusted for multiple testing using the Benjamini–Hochberg false discovery rate (FDR) method.
* Proteins are classified as significantly deregulated based on:

  Adjusted p-value threshold (default: 0.05)

  Absolute log2 fold change threshold (1.5)

**3. Data Visualization**

* **Volcano Plot**

   * A volcano plot is generated to visualize statistical significance versus magnitude of change:

   x-axis: log2 fold change

   y-axis: −log10(adjusted p-value)

   * Significant proteins are highlighted and optionally labeled.

* **Heatmap**

   * A heatmap of the DEPs is generated using hierarchical clustering, enabling visualization of expression patterns across all samples.

* **Principal Component Analysis (PCA)**

   * PCA is performed on the expression matrix of significant proteins to assess global separation between sensitive and resistant samples.

**4. Pathway Enrichment Analysis**

**Reactome Pathway Analysis**

* Significant proteins are mapped to Reactome pathways using the Reactome Content Service REST API.
* Enrichment is calculated by comparing observed protein counts per pathway against background expectations.

**KEGG Pathway Analysis**

* KEGG pathway enrichment is performed using KEGG gene mappings and over-representation analysis.
* Pathways are ranked by enrichment significance.

Both pathway analyses help identify:

* DNA damage response pathways
* Drug metabolism and transport pathways
* Associated signaling pathways

**5. Protein–Protein Interaction Network Analysis**

**STRING Database Integration**

* Protein–protein interactions are retrieved using the STRING API, filtered by a confidence score threshold.
* The resulting interaction network is constructed using networkx.

**Network Statistics**

The following network-level metrics are computed:

* Number of nodes and edges
* Average node degree
* Network density
* Connected components
* Degree centrality
* Highly connected hub proteins are identified as potential key regulators or therapeutic targets.

**6. Motif Enrichment Analysis**

* Sequence motifs enriched among deregulated proteins are identified using motif-scanning outputs (e.g., MEME, FIMO, or curated motif databases).

Motifs are:

* Tested for statistical enrichment
* Annotated with known functional roles (e.g., phosphorylation sites, interaction motifs)
* Linked to potential regulatory mechanisms or protein–protein interactions

## 6. 📤 Output files

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

## 🏃‍♀️ Installations

Python v3.9+

pip install pandas numpy scipy statsmodels matplotlib seaborn scikit-learn

pip install gseapy networkx requests biopython reportlab

## 📖 Dataset Used

Investigating Cisplatin Resistance in Squamous Cervical Cancer: Proteomic Insights into DNA Repair Pathways and Omics-Based Drug Repurposing

Amrita Mukherjee, Sayan Manna, Avinash Singh, Adrija Ray, and Sanjeeva Srivastava

Journal of Proteome Research 2025 24 (6), 2628-2642

DOI:10.1021/acs.jproteome.4c00885

The dataset used here is Table S4.

## 💻🔨 Applicability and Extensibility

Although developed for drug-sensitive vs resistant cell lines, this pipeline can be applied to:

* Other drug resistance models (cell lines/ tissue datasets)
* Paired cancer proteomics datasets
* Label-free or isobaric-tag–based quantitative proteomics
* Knockdown vs Wildtype
* By modifying column names and thresholds, the pipeline can be adapted to a wide range of comparative proteomics studies.

This code does not generalize-
* Paired vs unpaired design
* Column naming inconsistency
* Gene/Protein ID inconsistency

## 🔓 Summary

This pipeline provides an end-to-end, reproducible framework for integrating statistical proteomics, pathway biology, network analysis, and motif-level interpretation. It bridges wet-lab proteomics with computational systems biology, enabling mechanistic insights into drug resistance.
