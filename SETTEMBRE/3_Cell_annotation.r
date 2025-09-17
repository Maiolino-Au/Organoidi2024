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

suppressPackageStartupMessages(library(GPTCelltype))
suppressPackageStartupMessages(library(openai))

timepoints = c("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month")
path_to_data = "/sharedFolder/Data/"
dir_save = "Results_september/"

housekeeping_genes = c("ACTB", "DLG4")
genes_of_interest = c("SRCIN1", "KIAA1217", "CIT")


Sys.setenv(OPENAI_API_KEY = "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")



for (file_name in timepoints[2]) {
    dir_results <- paste0(dir_save, file_name)     
    if (!dir.exists(dir_results)) {dir.create(dir_results)}  
    
    sc_data_clustered <- readRDS(paste0(dir_save, file_name, "/", file_name, "_SeuratObj.rds"))
    cluster_markers <- read.csv(paste0(dir_save, file_name, "/", file_name, "_cluster_markers.csv"))

    res <- gptcelltype(cluster_markers, model = "gpt-4-0613")
    res <- gsub("^[0-9]{1,2}\\.\\s*", "", res)

    sc_data_clustered@meta.data$celltype_GPT <- as.factor(res[as.character(Idents(sc_data_clustered))])
    cluster_markers <- cluster_markers %>%
      mutate(celltype_GPT = res[as.character(cluster)])

    umap_plot <- DimPlot(sc_data_clustered, group.by = "celltype_GPT")
    ggsave(filename = paste0(dir_results, "/", file_name,"_UMAP_GPTcelltype.png"), plot = umap_plot, width = 1920*2, height = 1980*2, units = "px")

    de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
    write.csv(de_genes, file = paste0(dir_results, "/", file_name,"_de_genes_GPTcelltype.csv"))
    
    ggsave(filename = paste0(dir_results, "/", file_name,"_FeaturePlot.png"),
           plot = FeaturePlot(sc_data_clustered, c(genes_of_interest, housekeeping_genes[1])),  
           width = 1920*2, height = 1980*2, units = "px")

    # Get UMAP coordinates
    umap <- Embeddings(sc_data_clustered, reduction = "umap")
    # Get expression values for your genes
    expr <- FetchData(sc_data_clustered, vars = c(genes_of_interest, housekeeping_genes[1]))
    # Combine them into one data frame
    expr$UMAP_1 <- umap[,1]
    expr$UMAP_2 <- umap[,2]
    a <- list()
    for (i in 1:4) {
        df <- data.frame(value = expr[, i])
        a[[i]] <- ggplot(df, aes(x = value)) +
          geom_density(fill = "skyblue", alpha = 0.5) +
          labs(title = colnames(expr[i]), x = "Value", y = "Density")
    }
    density_plot <- do.call(plot_grid, c(a, labels = "AUTO", ncol = 2))

    ggsave(filename = paste0(dir_results, "/", file_name,"_UMAP_GPTcelltype.png"), 
           plot = umap_plot, width = 1920*2, height = 1980*2, units = "px")

    # Save
    saveRDS(sc_data_clustered, file = paste0(dir_results, "/", file_name, "_SeuratObj.rds"))
    write.csv(sc_data_clustered@meta.data, file = paste0(dir_results, "/", file_name, "_SeuratObj_meatadata.csv"))
    write.csv(cluster_markers, file = paste0(dir_results, "/", file_name,"_cluster_markers.csv"))
}

