# GSE58095-Data-Normalization-and-Subtype-Analysis
GSE58095 Data Normalization and Subtype Analysis
# GSE58095: Systemic Sclerosis Data Normalization and Subtype Analysis

## Overview
This repository contains an R workflow for processing and analyzing the **GSE58095** gene expression dataset. The study investigates gene expression signatures in patients with Systemic Sclerosis (SSc), specifically distinguishing between **Diffuse Cutaneous Scleroderma (dcSSc)**, **Limited Cutaneous Scleroderma (lcSSc)**, and **Healthy Controls**.

The script performs the entire preprocessing pipeline, starting from raw data retrieval from GEO (Gene Expression Omnibus) to dimensionality reduction using PCA.

**Platform:** Illumina HumanHT-12 V4.0 expression beadchip.

## Workflow Steps
The `GSE58095 Data Normalization and Subtype Analysis.R` script executes the following steps:

1.  **Data Acquisition:** - Automatically downloads the dataset using `GEOquery`.
    - Extracts expression matrix and clinical metadata.
2.  **Quality Control (Pre-Normalization):**
    - Checks for missing values.
    - Generates boxplots of raw expression data.
3.  **Normalization:**
    - Applies **Quantile Normalization** using the `limma` package to correct for technical variations between arrays.
    - Generates boxplots of normalized data for verification.
4.  **Filtering:**
    - Removes low-expressed genes (genes falling below the 25th percentile of mean expression) to reduce noise.
5.  **Phenotype Extraction:**
    - Parses metadata to classify samples into three distinct subtypes:
        - Healthy Control
        - Diffuse Scleroderma
        - Limited Scleroderma
6.  **Dimensionality Reduction:**
    - Performs **Principal Component Analysis (PCA)**.
    - Visualizes sample clustering based on disease subtypes.
7.  **Data Persistence:**
    - Saves the cleaned expression data, metadata, and subtype annotations for downstream analysis.

## Prerequisites
To run this script, you need **R** installed along with the following Bioconductor and CRAN packages:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))
install.packages(c("dplyr", "ggplot2"))
