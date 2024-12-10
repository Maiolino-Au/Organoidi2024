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


setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")
name_sc_data <- "days23"
path_to_data <- "Data/expression_23days"
# Specify the genes of interest OR use gene_oI()
genes_of_interest <- c("SRCIN1") 

name_new_dir <- paste(getwd(),"/Results/",name_sc_data,"_results",sep="")
if (dir.exists(name_new_dir)==FALSE) {
  dir.create(name_new_dir)
} 

cores_ram(4,6)
# _______________________________________________________________________________


umap_data_table <- read.delim("Data/umap_23days.txt")
umap_data_table <- umap_data_table[-1,]
for (i in 2:7) {if (i!=4) {umap_data_table[,i]<-as.numeric(umap_data_table[,i])}}

raw_sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data_SeuObj <- CreateSeuratObject(counts = raw_sc_data, min.cells = 3, min.features = 200, project = name_sc_data, names.delim = "-", names.field = 2)
# add the percentage of these mitochondrial or ribosomal genes to the meta.data
sc_data_SeuObj[["percent.mito"]] <- PercentageFeatureSet(object = sc_data_SeuObj, pattern = "^mt-")
sc_data_SeuObj[["percent.ribo"]] <- PercentageFeatureSet(object = sc_data_SeuObj, pattern = "Rps|Rpl|Mrpl|Mrps")
sc_data_SeuObj <- subset(x = sc_data_SeuObj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40 )

