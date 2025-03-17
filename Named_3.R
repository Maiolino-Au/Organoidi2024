setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

# setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

# install.packages('devtools')
# library(devtools)
# devtools::install_github('immunogenomics/presto')

library(presto)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

name_new_dir_0 <- paste(getwd(),"/Results/Cluster_from_paper",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,64)




i <- 1

print(paste(timepoints[i],"Starting"))

path_to_data <- paste("Data/expression_",timepoints[i],sep="")

# ___                     _________               ______               Prova |
#                                                                            V 

# name_new_dir_1 <- paste(name_new_dir_0,"/Named_CellType_",timepoints[i],"_Prova",sep="")
name_new_dir_1 <- paste(name_new_dir_0,"/Named_CellType_",timepoints[i],sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 

sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)


UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[i],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (j in 2:7) {if (j!=4) {UMAP_coordinates[,j]<-as.numeric(UMAP_coordinates[,j])}}

# Create a named vector of clusters
cluster_map <- setNames(UMAP_coordinates$CellType, UMAP_coordinates$NAME)

# Assign external clusters to the Seurat object
sc_data$external_cluster <- cluster_map[colnames(sc_data)]

# Set 'external_cluster' as the active identity
Idents(sc_data) <- sc_data$external_cluster

# Ensure 'external_cluster' is a factor with ordered levels
sc_data$external_cluster <- factor(sc_data$external_cluster, levels = sort(unique(sc_data$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data) <- sc_data$external_cluster

# Save
save(sc_data, file = paste(name_new_dir_1,"/Named_CellType",timepoints[i],".Robj", sep=""))

print(paste(timepoints[i],"FindAllMarkers"))

cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,   # Considera solo i marker espressi positivamente
                                  min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25)  # Soglia minima di LogFC

# Save
save(cluster_markers, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_CellType_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_top_100_markers_CellType.csv",sep=""))

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
# Save
write.csv(de_genes, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_",genes_of_interest,"_CellTypeExpression.csv",sep=""))

print(paste(timepoints[i],"End"))


i <- 6

print(paste(timepoints[i],"Starting"))

path_to_data <- paste("Data/expression_",timepoints[i],sep="")

name_new_dir_1 <- paste(name_new_dir_0,"/Named_CellType_",timepoints[i],sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 

sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)


UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[i],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (j in 2:7) {if (j!=4) {UMAP_coordinates[,j]<-as.numeric(UMAP_coordinates[,j])}}

# Create a named vector of clusters
cluster_map <- setNames(UMAP_coordinates$CellType, UMAP_coordinates$NAME)

# Assign external clusters to the Seurat object
sc_data$external_cluster <- cluster_map[colnames(sc_data)]

# Set 'external_cluster' as the active identity
Idents(sc_data) <- sc_data$external_cluster

# Ensure 'external_cluster' is a factor with ordered levels
sc_data$external_cluster <- factor(sc_data$external_cluster, levels = sort(unique(sc_data$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data) <- sc_data$external_cluster

# Save
save(sc_data, file = paste(name_new_dir_1,"/Named_CellType",timepoints[i],".Robj", sep=""))

print(paste(timepoints[i],"FindAllMarkers"))

cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,   # Considera solo i marker espressi positivamente
                                  min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25)  # Soglia minima di LogFC

# Save
save(cluster_markers, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_CellType_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_top_30_markers_CellType.csv",sep=""))

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
# Save
write.csv(de_genes, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_",genes_of_interest,"_CellTypeExpression.csv",sep=""))

print(paste(timepoints[i],"End"))


i <- 7

print(paste(timepoints[i],"Starting"))

path_to_data <- paste("Data/expression_",timepoints[i],sep="")

name_new_dir_1 <- paste(name_new_dir_0,"/Named_CellType_",timepoints[i],sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 

sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)


UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[i],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (j in 2:7) {if (j!=4) {UMAP_coordinates[,j]<-as.numeric(UMAP_coordinates[,j])}}

# Create a named vector of clusters
cluster_map <- setNames(UMAP_coordinates$CellType, UMAP_coordinates$NAME)

# Assign external clusters to the Seurat object
sc_data$external_cluster <- cluster_map[colnames(sc_data)]

# Set 'external_cluster' as the active identity
Idents(sc_data) <- sc_data$external_cluster

# Ensure 'external_cluster' is a factor with ordered levels
sc_data$external_cluster <- factor(sc_data$external_cluster, levels = sort(unique(sc_data$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data) <- sc_data$external_cluster

# Save
save(sc_data, file = paste(name_new_dir_1,"/Named_CellType",timepoints[i],".Robj", sep=""))

print(paste(timepoints[i],"FindAllMarkers"))

cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,   # Considera solo i marker espressi positivamente
                                  min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25)  # Soglia minima di LogFC

# Save
save(cluster_markers, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_CellType_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_top_30_markers_CellType.csv",sep=""))

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
# Save
write.csv(de_genes, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_",genes_of_interest,"_CellTypeExpression.csv",sep=""))

print(paste(timepoints[i],"End"))


i <- 8

print(paste(timepoints[i],"Starting"))

path_to_data <- paste("Data/expression_",timepoints[i],sep="")

name_new_dir_1 <- paste(name_new_dir_0,"/Named_CellType_",timepoints[i],sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 

sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)


UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[i],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (j in 2:7) {if (j!=4) {UMAP_coordinates[,j]<-as.numeric(UMAP_coordinates[,j])}}

# Create a named vector of clusters
cluster_map <- setNames(UMAP_coordinates$CellType, UMAP_coordinates$NAME)

# Assign external clusters to the Seurat object
sc_data$external_cluster <- cluster_map[colnames(sc_data)]

# Set 'external_cluster' as the active identity
Idents(sc_data) <- sc_data$external_cluster

# Ensure 'external_cluster' is a factor with ordered levels
sc_data$external_cluster <- factor(sc_data$external_cluster, levels = sort(unique(sc_data$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data) <- sc_data$external_cluster

# Save
save(sc_data, file = paste(name_new_dir_1,"/Named_CellType",timepoints[i],".Robj", sep=""))

print(paste(timepoints[i],"FindAllMarkers"))

cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,   # Considera solo i marker espressi positivamente
                                  min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25)  # Soglia minima di LogFC

# Save
save(cluster_markers, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_CellType_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_top_30_markers_CellType.csv",sep=""))

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
# Save
write.csv(de_genes, file = paste(name_new_dir_1,"/Named_",timepoints[i],"_",genes_of_interest,"_CellTypeExpression.csv",sep=""))

print(paste(timepoints[i],"End"))




