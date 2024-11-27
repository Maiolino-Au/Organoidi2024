library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

gene_oI <- function() {
  genes_of_interest <- readline(prompt="Gene of interest: ")
  de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) # Filter markers for genes of interest
  print(de_genes) # View the differential expression results
  goi_plot_UMAP <- FeaturePlot(sc_data_UMAP, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(UMAP)"))
  goi_plot_tSNE <- FeaturePlot(sc_data_tSNE, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(tSNE)"))
  CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2) 
}

# Set working directory
cat("Current directory:", getwd())
i <- readline(prompt="Change directory? [y/n] ")
if (i=="y") {
  work_dir <- readline(prompt="Enter your working directory: ")
  setwd(work_dir)
  cat("Current directory:", getwd())
}

name_sc_data <- readline(prompt="Project/session name? ")
path_to_data <- readline(prompt="path = ")

# OR

cores_ram(,)
setwd("")
name_sc_data <- ""
path_to_data <- ""

# _______________________________________________________________________________

raw_sc_data <- Read10X(data.dir = path_to_data) # gene.column = 1 -> Ensembl Gene IDs; = 2 -> Gene Symbols (default)

# Create Seurat object
sc_data_SeuObj <- CreateSeuratObject(counts = raw_sc_data, min.cells = 3, min.features = 200, project = name_sc_data, names.delim = "-", names.field = 2)


# add the percentage of these mitochondrial or ribosomal genes to the meta.data
sc_data_SeuObj[["percent.mito"]] <- PercentageFeatureSet(object = sc_data_SeuObj, pattern = "^mt-")
sc_data_SeuObj[["percent.ribo"]] <- PercentageFeatureSet(object = sc_data_SeuObj, pattern = "Rps|Rpl|Mrpl|Mrps")

sc_data_SeuObj <- subset(x = sc_data_SeuObj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40 )


# Normalize the data
sc_data <- NormalizeData(sc_data_SeuObj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)

# Scale the data
sc_data_SD <- ScaleData(sc_data)
sc_data <- sc_data_SD

# Perform dimensionality reduction (PCA)
sc_data <- RunPCA(sc_data, npcs = 50, verbose = FALSE)

ElbowPlot(object = sc_data, ndims = 50)

# Cluster the cells
sc_data <- FindNeighbors(sc_data, dims = 1:40)
sc_data <- FindClusters(sc_data, resolution = 0.6)

# Perform non-linear dimensionality reduction (UMAP/t-SNE)
sc_data_UMAP <- RunUMAP(sc_data, dims = 1:40)
sc_data_tSNE <- RunTSNE(sc_data, dims = 1:40)

# Visualization of clusters
plot_UMAP <- DimPlot(sc_data_UMAP, reduction = "umap", label = TRUE, pt.size = 1) + ggtitle("UMAP of Clusters")
plot_tSNE <- DimPlot(sc_data_tSNE, reduction = "tsne", label = TRUE, pt.size = 1) + ggtitle("tSNE of Clusters")
CombinePlots(plot=list(plot_UMAP,plot_tSNE), ncol=2)

# Perform differential expression analysis
cluster_markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Specify the genes of interest OR use gene_oI()
genes_of_interest <- c("") 

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) 

# View the differential expression results
print(de_genes)

# Visualize expression of specific genes
goi_plot_UMAP <- FeaturePlot(sc_data_UMAP, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(UMAP)"))
goi_plot_tSNE <- FeaturePlot(sc_data_tSNE, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(tSNE)"))
CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2)








# END
plan("multisession", workers = 1)
options(future.globals.maxSize = 0.5 * 1024^3)
