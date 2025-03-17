# setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

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

# install.packages('enrichR')
library(enrichR)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

name_new_dir_0 <- paste(getwd(),"/Results/MARZO",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")


cores_ram(1,64)

# __________________________________________________________________________________________

y <- 1
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_CellType",timepoints[y],".Robj", sep = ""))
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_",timepoints[y],"_CellType_markers.Robj", sep = ""))

top100 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()

# __________________________________________________________________________________
# ___________________________________ Annotation ___________________________________
# __________________________________________________________________________________

# ___________________________________ ChatGPT ___________________________________

CellType <- unique(top100$cluster)
top_genes_per_cell <- data.frame()

for (i in 1:length(CellType)) {
  top_genes_per_cell[i,1] <- CellType[i]
  top_genes_per_cell[i,2:101] <- c(top100[top100$cluster==CellType[i],7])
}


for (j in 1:length(CellType)) {
  testo <- "Hi, this is a list of the 100 most diffrentially expressed genes of a cluster derived from the single cell RNA sequencing of cells taken from a brain organid. Can you give me your 3 best guesses wich cell types are they? Write giust the 3 guesses in decreasing order of probability. "
  b <- ""
  if (b!=testo) {b<-testo}
  for (i in 2:101) {
    b <- paste(b,top_genes_per_cell[j,i])
  }
  print(CellType[j])
  print(b)
  print("----------------------------------------------")
}


# ___________________________________ EnirchR ___________________________________

# List available databases
dbs <- listEnrichrDbs()
# print(dbs)


gene_list <- top100[top100$cluster==CellType[1],]
gene_list <- gene_list[order(gene_list$avg_log2FC, decreasing = TRUE),] 
gene_list <- c(gene_list[,7])



# Perform enrichment analysis
enriched <- enrichr(gene_list, databases = c("PanglaoDB_Augmented_2021", 
                                             "ARCHS4_Tissues", 
                                             "Allen_Brain_Atlas_10x_scRNA_2021",
                                             "KEGG_2021_Human",
                                             "CellMarker_Augmented_2021",
                                             "Azimuth_Cell_Types_2021"))

# View results
print(enriched$PanglaoDB_Augmented_2021)

# __________________________________________________________________________________
# ___________________________________ Comparison ___________________________________
# __________________________________________________________________________________









