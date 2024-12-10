library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")
name_sc_data <- "days23"
path_to_data <- "Data/expression_23days"
# Specify the genes of interest OR use gene_oI()
#genes_of_interest <- c("SRCIN1") 

name_new_dir <- paste(getwd(),"/Results/",name_sc_data,"_Replicate_paper",sep="")
if (dir.exists(name_new_dir)==FALSE) {
  dir.create(name_new_dir)
} 

cores_ram(4,6)

# _______________________________________________________________________________


# Load data
raw_sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 

# Create Seurat object
sc_data_SeuObj <- CreateSeuratObject(counts = raw_sc_data, min.cells = 3, min.features = 500, project = name_sc_data, names.delim = "-", names.field = 2)

# Normalize the data
sc_data <- NormalizeData(sc_data_SeuObj, normalization.method = "LogNormalize", scale.factor = 1e6)

# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)

# # Scale the data
# sc_data_SD <- ScaleData(sc_data)
# sc_data <- sc_data_SD

# Score cells for cell cycle phases
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
sc_data <- CellCycleScoring(sc_data, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

# Regress out unwanted sources of variation (total UMI count and cell cycle effects)
sc_data <- ScaleData(sc_data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

# Perform dimensionality reduction (PCA)
sc_data <- RunPCA(sc_data, npcs = 30, verbose = FALSE)

# Cluster the cells
sc_data <- FindNeighbors(sc_data, dims = 1:30)
sc_data <- FindClusters(sc_data, resolution = 2)

# Perform non-linear dimensionality reduction (UMAP/t-SNE)
sc_data_UMAP <- RunUMAP(sc_data, dims = 1:30)

plot_UMAP <- DimPlot(sc_data_UMAP, reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(paste("UMAP of Clusters -",name_sc_data))
plot_UMAP






