# LOAD LIBRARIES
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(presto))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tictoc))

suppressPackageStartupMessages(library(enrichR))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

suppressPackageStartupMessages(library(SingleR))


total_time <- function(seconds) {
    d <- seconds %/% (86400)
    h <- (seconds %% 86400) %/% 3600
    m <- (seconds %% 3600) %/% 60
    s <- seconds %% 60
    
    cat(sprintf("Total Time: %f Days, %f Hours, %f Minutes and %f Seconds\n", d, h, m, s))
}

load.data <- function(
    file_name,
    data_path = path_to_data
) {
    # Load the data
    data <- Read10X(data.dir = paste(data_path, "expression_", file_name, sep = ""), gene.column = 1)

    return(data)
}

preprocessing <- function(
    data,
    file_name,
    output = F
) {
    # Create Seurat object
    sc_data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 500, project = file_name, names.delim = "-", names.field = 2)

    # Normalize the data
    sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = output)

    # Find variable features
    sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000, verbose = output)

    # Scale the data
    sc_data <- ScaleData(sc_data, verbose = output)

    return(sc_data)
}

PCA.cluster <- function(
    data, 
    file_name,
    res = 1, 
    n_dim = 20, 
    save_dir = F,
    save_add_on = "",
    output = F   
) {
    # PCA
    data <- RunPCA(data, npcs = max(n_dim), verbose = output)
    data <- RunUMAP(data, dims = n_dim, verbose = output)
    #print(ElbowPlot(object = data, ndims = max(50, n_dim)))

    # Cluster the cells
    data <- FindNeighbors(data, dims = n_dim, verbose = output)
    data <- FindClusters(data, resolution = res, verbose = output)
    
    #print(table(Idents(data)))

    # Save the parial
    if (save_dir != F) {
        save(data, file = paste0(save_dir, "/Clustered", file_name, save_add_on, ".Robj"))
    }
    
    return(data)
}

cluster.plot <- function(
    data,
    file_name = timepoints[time_point],
    pt.size = 1,
    print_plot = F
) { 
    plots <- list(
        UMAP = DimPlot(data, reduction = "umap", label = TRUE, pt.size = pt.size) + 
            ggtitle(paste("UMAP -",file_name)),
        PCA = DimPlot(data, reduction = "pca", label = TRUE, pt.size = pt.size) + 
            ggtitle(paste("PCA -",file_name))
    )

    if (print_plot) {print(plots)}
    
    return(plots)
}

cluster.markers <- function(
    data,
    output = F
) {
    # Find all markers for every cluster compared to all remaining cells
    markers <- FindAllMarkers(data,
                              only.pos = TRUE,   # Considera solo i marker espressi positivamente
                              min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                              logfc.threshold = 0.25,  # Soglia minima di LogFC
                              verbose = output)
        
    return(markers)
}

top.genes <- function(
    clmark = cluster_markers,
    num = 50
) {
    top_genes_by_cluster <- cluster_markers %>% group_by(cluster) %>% top_n(n = top_n_genes, wt = avg_log2FC) %>% as.data.frame()

    return(top_genes_by_cluster)
}
