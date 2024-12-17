library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

cores_ram(2,6)
setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")
name_sc_data <- "days23"
path_to_data <- "Data/expression_23days"
# Specify the genes of interest OR use gene_oI()
genes_of_interest <- c("SRCIN1") 

# _______________________________________________________________________________

name_new_dir <- paste(getwd(),"/Results/",name_sc_data,"_Replicate_paper",sep="")
if (dir.exists(name_new_dir)==FALSE) {
  dir.create(name_new_dir)
} 

raw_sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data_SeuObj <- CreateSeuratObject(counts = raw_sc_data, min.cells = 3, min.features = 500, project = name_sc_data, names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data_SeuObj, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)

umap_23d <- read.delim("Data/umap_23days.txt")
umap_23d <- umap_23d[-1,]
for (i in 2:7) {if (i!=4) {umap_23d[,i]<-as.numeric(umap_23d[,i])}}


# seurat_barcodes <- colnames(sc_data)
# umap_barcodes <- umap_23d$NAME
# 
# all(umap_barcodes %in% seurat_barcodes)  # Should return TRUE


# my_clusters <- umap_23d$Cluster
# 
# sc_data_ex_cl <- sc_data
# sc_data_ex_cl$external_cluster <- my_clusters
# table(Idents(sc_data_ex_cl))


# Create a named vector of clusters
cluster_map <- setNames(umap_23d$Cluster, umap_23d$NAME)

# Assign external clusters to the Seurat object
sc_data_ex_cl <- sc_data
sc_data_ex_cl$external_cluster <- cluster_map[colnames(sc_data)]

table(sc_data_ex_cl$external_cluster)

# Set 'external_cluster' as the active identity
Idents(sc_data_ex_cl) <- sc_data_ex_cl$external_cluster

# # Verify the active identities are correct
# table(Idents(sc_data_ex_cl))  # Should match external_cluster


# Ensure 'external_cluster' is a factor with ordered levels
sc_data_ex_cl$external_cluster <- factor(sc_data_ex_cl$external_cluster, levels = sort(unique(sc_data_ex_cl$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data_ex_cl) <- sc_data_ex_cl$external_cluster

# Verify the order
table(Idents(sc_data_ex_cl))  # Should now match the increasing numeric order



cores_ram(2,6)
cluster_markers <- FindAllMarkers(sc_data_ex_cl, 
                          only.pos = TRUE,   # Considera solo i marker espressi positivamente
                          min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                          logfc.threshold = 0.25)  # Soglia minima di LogFC


# Save
save(cluster_markers, file = paste(name_new_dir,"/",name_sc_data,"_cluster_markers.Robj", sep=""))

top_markers_all <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# Save
write.csv(top_markers_all, file = paste(name_new_dir,"/",name_sc_data,"_top_30_markers.csv",sep=""))

# Filter markers for genes of interest
de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) 

# View the differential expression results
print(de_genes)
# Save
write.csv(de_genes, file = paste(name_new_dir,"/",name_sc_data,"_",genes_of_interest,"ClusterExpression.csv",sep=""))
