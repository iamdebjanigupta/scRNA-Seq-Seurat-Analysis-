#Computational Exploration of Single-Cell RNA-seq Data Using Seurat

##Author:
Debjani Gupta
MSc. Forensic Science
National Forensic Sciences University

##Project Context:
This mini-project was completed as part of a 6-Week Hands-on Training Program on Next-Generation Sequencing (NGS) conducted by BIOTECHTREK.
This work is not an academic submission to NFSU but an independent computational training project.

##Project Overview
Single-cell RNA sequencing (scRNA-seq) enables transcriptomic profiling at the resolution of individual cells, allowing the study of cellular heterogeneity in complex tissues. In this project, an end-to-end scRNA-seq workflow was implemented using the Seurat package in R to analyze human immune cell data. The analysis included:
-Quality control (QC) and filtering of low-quality cells
-Normalization and identification of highly variable genes
-Dimensionality reduction using Principal Component Analysis (PCA)
-Graph-based clustering and UMAP visualization
-Cell type annotation using canonical marker genes
-Differential gene expression analysis between biologically relevant groups
This project demonstrates a reproducible computational pipeline for analyzing single-cell transcriptomic data.

##Dataset:
Nature (2025) Study: https://www.nature.com/articles/s41586-025-09686-5
CellxGene Dataset Link: https://datasets.cellxgene.cziscience.com/36a22a49-8532-43d4-afa2-d90c894a6c7a.h5ad

##Software and Tools Used:
-R version 4.5.2 (UCRT)
-Seurat
-zellkonverter
-dplyr
-ggplot2
-EnhancedVolcano

##Repository Structure:
figures/        → QC plots, PCA diagnostics, UMAP, marker heatmaps, DE plots  
tables/         → Cluster marker genes and differential expression results (CSV)  
scripts/        → R analysis script (analysis.R)  
report/         → Final project report (PDF)  
session_info/   → R session information for reproducibility  

##Key Results:
-QC filtering removed low-quality cells and potential doublets.
-15 principal components were selected based on the Elbow plot.
-Clustering at resolution 0.5 identified biologically meaningful immune cell populations.
-Canonical markers enabled confident cell type annotation.
-Differential expression analysis revealed functionally relevant genes.
-A rare NK cell subset and mild expression gradients were observed.

##Reproducibility:
-All analysis steps are documented in scripts/analysis.R.
-Session information is provided in session_info/ to ensure reproducibility.

##References:
-Seurat Tutorial: https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/Seurat_QC_to_Clustering/
-Nature (2025) Study
-CellxGene Data Portal

