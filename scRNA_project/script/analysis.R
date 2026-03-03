
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install zellkonverter from Bioconductor
BiocManager::install("zellkonverter")
?readH5AD
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
library(zellkonverter)  # load package
library(Seurat)         # load Seurat
# Read the dataset (replace path if different)
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# Convert to Seurat object
seurat_obj <- as.Seurat(sce, counts = "X")
# Load all libraries manually
library(Seurat)
library(SeuratDisk)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
# Optional: EnhancedVolcano if you want the volcano plot
# library(EnhancedVolcano)
PackageCheck("ggplot2")
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
seurat_obj <- as.Seurat(sce, counts = "X")
library(Seurat)
library(SeuratDisk)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
setwd("C:/Users/debja/OneDrive/Documents/scRNA_project")
if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("tables")) dir.create("tables")
if(!dir.exists("session_info")) dir.create("session_info")
sce <- readH5AD("dataset.h5ad")
seurat_obj <- as.Seurat(sce, counts = "X")
# 1️⃣ Load packages
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(SeuratObject)  # optional, but safe
# 2️⃣ Read h5ad
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# 3️⃣ Convert manually to Seurat without triggering PackageCheck()
seurat_obj <- CreateSeuratObject(
  counts = counts(sce),          # raw counts
  meta.data = as.data.frame(colData(sce))
)
# Load zellkonverter and read dataset
library(zellkonverter)
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# List assays
assayNames(sce)
# Check first few rows of first assay
head(assay(sce, 1))
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
# 1️⃣ Load libraries
library(zellkonverter)   # to read h5ad
library(Seurat)          # for Seurat object
library(SingleCellExperiment)
# 2️⃣ Read your h5ad file
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# 3️⃣ Check what assays exist in the dataset
assayNames(sce)
# Example output might be: "X" "logcounts"
# Use the one that contains raw counts ("X" usually)
# 4️⃣ Create Seurat object safely
seurat_obj <- CreateSeuratObject(
  counts = assay(sce, "X"),           # use correct assay name from previous step
  meta.data = as.data.frame(colData(sce))
)
# 5️⃣ Add mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
# 1️⃣ Load h5ad
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# 2️⃣ Check what assays exist
assayNames(sce)
# Example output: "X" "logcounts"
# "X" is usually raw counts
# 3️⃣ Create Seurat object manually
seurat_obj <- CreateSeuratObject(
  counts = assay(sce, "X"),          # use assay with raw counts
  meta.data = as.data.frame(colData(sce))
)
# 4️⃣ Add mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
# Read the dataset
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# Check what assays exist
assayNames(sce)
# Example output: "X" "raw" "logcounts"
# Use the one with actual raw counts, usually "X"
# Extract counts matrix safely
counts_matrix <- as.matrix(assay(sce, "X"))
# Optional: check dimensions
dim(counts_matrix)
# Create Seurat object from the raw counts
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = as.data.frame(colData(sce))
)
# Add percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
library(zellkonverter)
library(Seurat)
# 1️⃣ Read h5ad using zellkonverter
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# 2️⃣ Check assays
assayNames(sce)
# Example output: "X" "raw" "logcounts"
# 3️⃣ Extract raw counts as a matrix
counts_matrix <- as.matrix(assay(sce, "X"))  # Replace "X" with correct assay if different
library(zellkonverter)
library(Seurat)
library(Matrix)  # for sparse matrices
# 1️⃣ Load h5ad
sce <- readH5AD("C:/Users/debja/OneDrive/Documents/scRNA_project/dataset.h5ad")
# 2️⃣ Check assays
assayNames(sce)
# Suppose raw counts are in "X"
# 3️⃣ Extract as sparse matrix
counts_sparse <- as(assay(sce, "X"), "dgCMatrix")  # memory-efficient
# 4️⃣ Extract metadata
meta_df <- as.data.frame(colData(sce))
# 5️⃣ Create Seurat object using sparse counts
seurat_obj <- CreateSeuratObject(counts = counts_sparse, meta.data = meta_df)
# 6️⃣ Add mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
# QC
VlnPlot(seurat_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
