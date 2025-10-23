# LOAD LIBRARIES
suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(future)
    library(ggplot2)
    library(dplyr)
    library(presto)
    library(cowplot)
    library(tictoc)
    library(SingleCellExperiment)
    library(SingleR)
    library(scuttle)
    library(zellkonverter)
    library(patchwork)
    library(SeuratDisk)
})


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
    print(paste(data_path, "expression_", file_name, sep = ""))
    data <- Read10X(data.dir = paste(data_path, "expression_", file_name, sep = ""), gene.column = 1)

    return(data)
}

preprocessing <- function(
    data,
    file_name,
    output = F
) {
    # Create Seurat object
    data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 500, project = file_name, names.delim = "-", names.field = 2)

    mito.genes <- grep(pattern = "^MT-", x = rownames(x = data), value = TRUE)
    percent.mito <- Matrix::colSums(GetAssayData(data, slot='counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(data, slot="counts"))
    data[['percent.mito']] <- percent.mito
    ribo.genes <- grep("^RP[S,L]",rownames(data), value = TRUE)
    percent.ribo <- Matrix::colSums(GetAssayData(data, slot='counts')[ribo.genes, ]) / Matrix::colSums(GetAssayData(data, slot="counts"))
    data[['percent.ribo']] <- percent.ribo

    # Normalize the data
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = output)

    # Find variable features
    data <- FindVariableFeatures(data, selection.method = "mvp", nfeatures = 2000, verbose = output)

    data = CellCycleScoring(data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    data$CC.Difference <- data$S.Score - data$G2M.Score

    # Scale the data
    data<-ScaleData(data, features=VariableFeatures(data),vars.to.regress=c("nCount_RNA","CC.Difference"), verbose = output)

    return(data)
}

PCA.cluster <- function(
    data, 
    file_name,
    res = 1, 
    n_dim = 1:20, 
    dir_save = F,
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
    if (dir_save != F) {
        saveRDS(data, file = paste0(dir_save, "/Clustered_", file_name, "_", save_add_on, ".rds"))
    }
    
    return(data)
}

plottamelo.tutto <- function(
    data,
    file_name,
    genes_of_interest = genes_of_interest,
    dir_save = dir_save,
    save_add_on = "",
    pt.size = 1,
    width = 1920*2,
    height = 1920*2,
    scale.dimensions = 1,
    print_plot = F
) { 
    umap <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = pt.size) + 
    ggtitle(paste("UMAP -",file_name))

    gene_plot <- FeaturePlot(data, features = genes_of_interest) &
    theme(aspect.ratio = 1)

    plot_blend <- umap | gene_plot 

    width = width*scale.dimensions
    height = height*scale.dimensions
    
    if(save_add_on != "") {file_name = paste(file_name, save_add_on, sep = "_")}
    
    ggsave(filename = paste0(dir_save, "Plot_", file_name,"_umap.png"), 
           plot = umap, 
           width = width, height = height, units = "px")
    ggsave(filename = paste0(dir_save, "Plot_", file_name, "_", genes_of_interest,"_umap.png"), 
           plot = gene_plot, 
           width = width, height = height, units = "px")
    ggsave(filename = paste0(dir_save, "Plot_", file_name, "_combined_umap&", genes_of_interest, ".png"), 
           plot = plot_blend, 
           width = width*2, height = height, units = "px")
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

annotate.singler <- fucntion(
    data, 
    file_name,
    ref,
    dir_save = F,
    output = F
) {
    results <- SingleR(
        test = as.SingleCellExperiment(data),
        ref = ref,
        labels = ref$cell_type
    )

    data$singler_labels <- results$labels

    saveRDS(ref, paste0(dir_save, "/Annotated_SingleR_", file_name, ".rds"))
    
    return(data)
}
