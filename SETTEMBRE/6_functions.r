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

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(patchwork))


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
        save(data, file = paste0(dir_save, "/Clustered_", file_name, "_", save_add_on, ".Robj"))
    }
    
    return(data)
}

plottamelo.tutto <- function(
    data,
    file_name,
    genes_of_interest = genes_of_interest,
    dir_save = dir_save,
    save_add_on,
    pt.size = 1,
    print_plot = F
) { 
    # umap
    umap <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = pt.size) + 
        ggtitle(paste("UMAP -",file_name))

    ggsave(filename = paste0(dir_save, "Plot_", file_name, "_", save_add_on,"_umap.png"), 
           plot = umap, 
           width = 1920*2, height = 1920*2, units = "px")

    plot_blend <- FeaturePlot(data, 
                              features = genes_of_interest, 
                              blend = T, 
                              ncol = 2, 
                              cols = c("gray90", "blue", "red")) &
      theme(aspect.ratio = 1)

    for (i in 1:4) {
        ggsave(filename = paste0(dir_save, "Plot_", file_name, "_", save_add_on,"_blend_", i, ".png"), 
               plot = plot_blend[[i]], 
               width = 1920*2, height = 1920*2, units = "px")
    }

    plot_blend <- (plot_blend[[1]] | plot_blend[[2]]) /
                  (plot_blend[[3]] | plot_blend[[4]]) +
                  plot_annotation(
                    title = paste(file_name, "-", genes_of_interest[1], "&", genes_of_interest[2], "Expression Blend"),
                    subtitle = "UMAP colored by expression levels",
                    caption = "Subset of clusters with avg_log2FC > 0.6 - Seurat FeaturePlot (blend = TRUE)"
                  )

    ggsave(
    filename = paste0(dir_save, "Plot_", file_name, "_", save_add_on, "_blend_0.png"),
    plot = plot_blend,
    width = 1920*4, height = 1920*2, units = "px"
)

    ggsave(
    filename = paste0(dir_save, "Plot_", file_name, "_", save_add_on,"_umap&blend.png"),
    plot = umap | plot_blend+
              plot_annotation(
                title = paste(file_name, "-", genes_of_interest[1], "&", genes_of_interest[2], "Expression Blend"),
                subtitle = "UMAP colored by expression levels",
                caption = "Subset of clusters with avg_log2FC > 0.6 - Seurat FeaturePlot (blend = TRUE)"
              ),
    width = 1920*4, height = 1920*2, units = "px"
)
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
