library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

# gene_oI <- function() {
#   genes_of_interest <- readline(prompt="Gene of interest: ")
#   de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) # Filter markers for genes of interest
#   print(de_genes) # View the differential expression results
#   goi_plot_UMAP <- FeaturePlot(sc_data_UMAP, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(UMAP)",name_sc_data))
#   goi_plot_tSNE <- FeaturePlot(sc_data_tSNE, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(tSNE)",name_sc_data))
#   CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2)
#   jpeg(paste("Results/",name_sc_data,"_",genes_of_interest,"_ClusterPlot",sep=""),width=1920, height=1080)
#   CombinePlots(plot=list(plot_UMAP,plot_tSNE), ncol=2)
#   dev.off()
# }

# # Set working directory
# cat("Current directory:", getwd())
# i <- readline(prompt="Change directory? [y/n] ")
# if (i=="y") {
#   work_dir <- readline(prompt="Enter your working directory: ")
#   setwd(work_dir)
#   cat("Current directory:", getwd())
# }
# 
# name_sc_data <- readline(prompt="Project/session name? ")
# path_to_data <- readline(prompt="path = ")

# OR

cores_ram(,)
setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")
name_sc_data <- ""
path_to_data <- ""
genes_of_interest <- c("SRCIN1") 

# _______________________________________________________________________________

if (dir.exists(paste(getwd(),"/Results/",name_sc_data,"_results",sep=""))==FALSE) {
  dir.create(paste(getwd(),"/Results/",name_sc_data,"_results",sep=""))
}

raw_sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 

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
table(Idents(sc_data))

# Save
save(sc_data, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_clustered.Robj", sep=""))

# Perform non-linear dimensionality reduction (UMAP/t-SNE)
sc_data_UMAP <- RunUMAP(sc_data, dims = 1:40)
sc_data_tSNE <- RunTSNE(sc_data, dims = 1:40)

# Save
save(sc_data_UMAP, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_UMAP.Robj", sep=""))
save(sc_data_tSNE, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_tSNE.Robj", sep=""))

# Visualization of clusters
plot_UMAP <- DimPlot(sc_data_UMAP, reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(paste("UMAP of Clusters -",name_sc_data))
plot_tSNE <- DimPlot(sc_data_tSNE, reduction = "tsne", label = TRUE, pt.size = 1) + ggtitle(paste("tSNE of Clusters -",name_sc_data))
CombinePlots(plot=list(plot_UMAP,plot_tSNE), ncol=2)

# Save
jpeg(paste("Results/",name_sc_data,"_results/",name_sc_data,"_ClusterPlot.jpg",sep=""),width=1920, height=1080)
CombinePlots(plot=list(plot_UMAP,plot_tSNE), ncol=2)
dev.off()

cores_ram(2,13)

# Perform differential expression analysis
cluster_markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

# Save
save(cluster_markers, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_cluster_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_top_30_markers.csv",sep=""))




# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) 

# View the differential expression results
print(de_genes)
# Save
write.csv(de_genes, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_",genes_of_interest,"ClusterExpression.csv",sep=""))

# Visualize expression of specific genes
goi_plot_UMAP <- FeaturePlot(sc_data_UMAP, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(UMAP)",name_sc_data))
goi_plot_tSNE <- FeaturePlot(sc_data_tSNE, features = genes_of_interest) + ggtitle(paste(genes_of_interest,"(tSNE)",name_sc_data))
CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2)

# Save
jpeg(paste("Results/",name_sc_data,"_results/",name_sc_data,"_",genes_of_interest,"_ClusterPlot.jpg",sep=""),width=1920, height=1080)
CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2)
dev.off()

# Save all
# save(sc_data, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_clustered.Robj", sep=""))
# save(sc_data_UMAP, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_UMAP.Robj", sep=""))
# save(sc_data_tSNE, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_tSNE.Robj", sep=""))
# jpeg(paste("Results/",name_sc_data,"_results/",name_sc_data,"_ClusterPlot.jpg",sep=""),width=1920, height=1080)
# CombinePlots(plot=list(plot_UMAP,plot_tSNE), ncol=2)
# dev.off()
# write.csv(top_markers_all, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_top_30_markers.csv",sep=""))
# save(cluster_markers, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_cluster_markers.Robj", sep=""))
# write.csv(de_genes, file = paste("Results/",name_sc_data,"_results/",name_sc_data,"_",genes_of_interest,"ClusterExpression.csv",sep=""))
# jpeg(paste("Results/",name_sc_data,"_results/",name_sc_data,"_",genes_of_interest,"_ClusterPlot.jpg",sep=""),width=1920, height=1080)
# CombinePlots(plot=list(goi_plot_UMAP,goi_plot_tSNE), ncol=2)
# dev.off()


# END
plan("multisession", workers = 1)
options(future.globals.maxSize = 0.5 * 1024^3)

