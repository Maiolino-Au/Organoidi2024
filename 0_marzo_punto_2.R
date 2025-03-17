setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

# setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)
library(tibble)

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

name_new_dir_0 <- paste(getwd(),"/Results/Punto2",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________

#_____ top 30
y<-1
top30 <- read.csv(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_",timepoints[y],"_top_30_markers_CellType.csv", sep=""))
# b<-c(top30[top30$cluster=="aRG",8])
# save(b,file = paste(getwd(),"/top_genes_23d_aRG.csv",sep=""))

top100 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()

CellType <- unique(top100$cluster)
top_genes_per_cell <- data.frame()

for (i in 1:length(CellType)) {
  top_genes_per_cell[i,1] <- CellType[i]
  top_genes_per_cell[i,2:101] <- c(top100[top100$cluster==CellType[i],7])
}

question <- function(x) {
  testo <- "Hi, this is a list of the 100 most diffrentially expressed genes of a cluster derived from the single cell RNA sequencing of cells taken from a brain organid. Can you give me your 3 best guesses wich cell types are they? Write giust the 3 guesses in decreasing order of probability. "
  b <- ""
  if (b!=testo) {b<-testo}
  for (i in 2:101) {
    b <- paste(b,top_genes_per_cell[x,i])
  }
  print(CellType[x])
  print(b)
  print("----------------------------------------------")
}

for (j in 1:length(CellType)) {
  question (j)
}


#_____ load find all markers

# markers <- load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[1],"/Named_CellType",timepoints[1],".Robj", sep = "")) %>% FindAllMarkers(,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)





# _______________________________________________________________________________


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


# _______________________________________________________________________________


# library(devtools)
# install_github("ctlab/fgsea")
# install.packages("msigdbr")
library(fgsea)
library(msigdbr)

gene_ranks <- data.frame(c(a[a$cluster==CellType[i],8]),c(a[a$cluster==CellType[i],3]))
gene_ranks <- sort(gene_ranks, decreasing = TRUE)  # Ensure sorted ranking


# Extract gene names and log2FC values
gene_ranks <- setNames(a[a$cluster == CellType[i], 3],  # Log fold changes
                       a[a$cluster == CellType[i], 8])  # Gene names

# Sort in descending order (required for GSEA)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

fgseaRes <- fgsea(pathways = gene_sets, 
                  stats = gene_ranks, 
                  # nperm = 1000,
                  scoreType = "pos"
)


# Load a reference database (e.g., MSigDB Cell Type Signatures)
msigdb <- msigdbr(species = "Homo sapiens", category = "C8")  # C8 = cell type signatures

# Convert to list format
gene_sets <- split(msigdb$gene_symbol, msigdb$gs_name)

# Run fgsea
fgseaRes <- fgsea(pathways = gene_sets, 
                  stats = gene_list,  # If you have fold changes
                  minSize = 10, maxSize = 500, 
                  nperm = 1000)

# Show top enriched cell types
head(fgseaRes[order(fgseaRes$padj), ])
