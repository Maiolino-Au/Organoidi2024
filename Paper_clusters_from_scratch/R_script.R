# LOAD LIBRARIES
library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)
library(presto)
#library(enrichR)
library(cowplot)


# SET UP NAMES
timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
housekeeping_genes <- c("ACTB","DLG4")
genes_of_interest <- c("SRCIN1","KIAA1217","CIT")
path_to_data <- "/sharedFolder/data_organoidi_velasco/Data/"

# sc_data <- ""

name_new_dir_0 <- paste(getwd(),"/Results",sep="")
if (!dir.exists(name_new_dir_0)) {
  dir.create(name_new_dir_0)
} 

name_new_dir_partial <- paste(getwd(),"/Partial",sep="")
if (!dir.exists(name_new_dir_partial)) {
  dir.create(name_new_dir_partial)
} 


load_data <- function(timepoint){
  # Load the data
  sc_data <- Read10X(data.dir = paste0(path_to_data, "expression_", timepoints[i],), gene.column = 1)
  
  # Create Seurat object
  sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
    
  # Normalize the data
  sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
  
  # Find variable features
  sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)
  
  return(sc_data)
}

# PCA&Clusterization 
PCA_cluster <- function(sc_data, res){
  # PCA
  sc_data <- RunPCA(sc_data, npcs = 50, verbose = FALSE)
  #print(ElbowPlot(object = sc_data, ndims = 50))
    
  # Cluster the cells
  sc_data <- FindNeighbors(sc_data, dims = 1:40)
  sc_data <- FindClusters(sc_data, resolution = res)
  
  #print(table(Idents(sc_data)))

  save(sc_data, file = paste(name_new_dir_partial, "/PCA_res_",res,"_",timepoints[i],".Robj", sep=""))
  return(sc_data)
}

# FIND ALL MARKERS
cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,   # Considera solo i marker espressi positivamente
                                  min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25)  # Soglia minima di LogFC

# Save the markers
save(cluster_markers, file = paste(name_new_dir_partial, "/cluster_markers_",timepoints[i],".Robj", sep=""))

de_genes_f <- function(genes_oi){
    # Find differentially expressed genes
    provv <- cluster_markers %>% filter(gene %in% genes_of_interest) 
    
    # Save the DE genes
    save(provv, file = paste(name_new_dir_partial, "/de_genes_",timepoints[i],".txt", sep=""))
    
    return(provv)
}

