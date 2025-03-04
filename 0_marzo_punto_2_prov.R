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

name_new_dir_0 <- paste(getwd(),"/Results/Punto2",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________

y<-1
a <- read.csv(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_",timepoints[y],"_top_30_markers_CellType.csv", sep=""))
b<-c(a[a$cluster=="aRG",8])
save(b,file = paste(getwd(),"/top_genes_23d_aRG.csv",sep=""))

CellType <- unique(a$cluster)
top_genes_per_cell <- data.frame()

for (i in 1:length(CellType)) {
  top_genes_per_cell[i,1] <- CellType[i]
  top_genes_per_cell[i,2:31] <- c(a[a$cluster==CellType[i],8])
}





# List available databases
dbs <- listEnrichrDbs()
print(dbs)


gene_list <- c(a[a$cluster==CellType[i],8])

# Perform enrichment analysis
enriched <- enrichr(gene_list, databases = c("PanglaoDB_Augmented_2021", "ARCHS4_Tissues"))

# View results
print(enriched$PanglaoDB_Augmented_2021)



library(Seurat)
install.packages("Azimuth")
library(Azimuth)

# Load reference for brain/neural cells
reference <- Azimuth::LoadReference("pbmc_reference")

# Run annotation
query <- MapQuery(
  anchorset = reference,
  query = my_seurat_object,
  reference = reference
)

# View results
query$predicted.celltype


install.packages("SingleR")
install.packages("celldex")
library(SingleR)
library(celldex)

# Load a reference dataset (e.g., Human Primary Cell Atlas)
ref <- celldex::HumanPrimaryCellAtlasData()

# Assume `my_gene_expression` is a matrix with gene expression data
pred <- SingleR(test = my_gene_expression, ref = ref, labels = ref$label.main)

# View predicted cell type
pred$labels

library(devtools)
install_github("ctlab/fgsea")
install.packages("msigdbr")
library(fgsea)
library(msigdbr)

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

