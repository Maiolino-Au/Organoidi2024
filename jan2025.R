setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

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

name_new_dir_0 <- paste(getwd(),"/Results/Cluster&CellTypes",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________



for (a in 1:length(timepoints)) {
  UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[a],".txt",sep=""))
  UMAP_coordinates <- UMAP_coordinates[-1,]
  for (i in 2:7) {if (i!=4) {UMAP_coordinates[,i]<-as.numeric(UMAP_coordinates[,i])}}
  
  CellType <- sort(unique(UMAP_coordinates$CellType))
  N_clusters <- max(as.numeric(UMAP_coordinates$Cluster))
  
  results <- matrix(0, ncol = N_clusters + 1, nrow = length(CellType))
  results <- data.frame(results)
  rownames(results) <- CellType
  colnames(results) <- 0:N_clusters
  results$Tot <- 0
  
  for (CT in 1:length(CellType)) {
    results[CT,21] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT],7])
    for (Cl in 0:19) {results[CT,Cl+1] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT] & UMAP_coordinates$Cluster==Cl,7])/results[CT,21]}
  }
  write.csv(results, file = paste(name_new_dir_0,"/",timepoints[a],"_Cluster&CellTypes.csv",sep=""))
}


# name <- timepoints[1]
# UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",name,".txt",sep=""))
# UMAP_coordinates <- UMAP_coordinates[-1,]
# for (i in 2:7) {if (i!=4) {UMAP_coordinates[,i]<-as.numeric(UMAP_coordinates[,i])}}
# 
# CellType <- sort(unique(UMAP_coordinates$CellType))
# N_clusters <- max(as.numeric(UMAP_coordinates$Cluster))
# 
# results <- matrix(0, ncol = N_clusters + 1, nrow = length(CellType))
# results <- data.frame(results)
# rownames(results) <- CellType
# colnames(results) <- 0:N_clusters
# results$Tot <- 0
# 
# for (CT in 1:length(CellType)) {
#   results[CT,21] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT],7])
#   for (Cl in 0:19) {results[CT,Cl+1] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT] & UMAP_coordinates$Cluster==Cl,7])/results[CT,21]}
# }













