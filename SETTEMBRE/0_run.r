source(file = "0_functions.r")

timepoints = c("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month")
housekeeping_genes = c("ACTB", "DLG4")
genes_of_interest = c("SRCIN1", "KIAA1217", "CIT")
path_to_data = "/sharedFolder/Data/"
dir_save = "Results_september/"
res_list = c(1, 0.5, 0.2, 0.1)#seq(0.1, 1, by = 0.1)
dim_list = c(1:20)
top_n_genes = 50

Sys.setenv(OPENAI_API_KEY = "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")

tic("Clustering")

# Seurat
for (file_name in timepoints[5:8]) {
    tic(file_name)

    dir_results <- paste0(dir_save, "final_", file_name)     
    if (!dir.exists(dir_results)) {dir.create(dir_results)} 

    cat(file_name,  "loaded")
    raw_data <- load.data(file_name)

    sc_data <- preprocessing(raw_data)

    sc_data_clustered <- PCA.cluster(sc_data, res = res_list[2], n_dim= dim_list)

    cluster_plots <- cluster.plot(
        sc_data_clustered,
        file_name = file_name,
        pt.size = 0.5,
        print_plot = F
    )
    # save
    for (nm in names(cluster_plots)) {
      ggsave(
        filename = paste0(dir_results, "/", file_name,"_",nm, ".png"),
        plot = cluster_plots[[nm]], 
        width = 1920*2, height = 1980*2, units = "px"
      )
    }

    cluster_markers <- cluster.markers(sc_data_clustered, output = F)

    top_genes <- top.genes(top_genes_by_cluster, num = 50)

    # GPT cell type annotation
    res <- gptcelltype(cluster_markers, model = "gpt-4-0613")
    res <- gsub("^[0-9]{1,2}\\.\\s*", "", res)

    sc_data_clustered@meta.data$celltype_GPT <- as.factor(res[as.character(Idents(sc_data_clustered))])
    cluster_markers <- cluster_markers %>%
      mutate(celltype_GPT = res[as.character(cluster)])
    
    # UMAP with GPT cell types
    umap_plot <- DimPlot(sc_data_clustered, group.by = "celltype_GPT")
    ggsave(filename = paste0(dir_results, "/", file_name,"_UMAP_GPTcelltype.png"), plot = umap_plot, width = 1920*2, height = 1980*2, units = "px")

    # Save de genes for genes of interest
    de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
    write.csv(de_genes, file = paste0(dir_results, "/", file_name,"_de_genes_GPTcelltype.csv"))
    
    # Feature plot
    ggsave(filename = paste0(dir_results, "/", file_name,"_FeaturePlot.png"),
           plot = FeaturePlot(sc_data_clustered, c(genes_of_interest, housekeeping_genes[1])),  
           width = 1920*2, height = 1980*2, units = "px")

    # Density plot
    umap <- Embeddings(sc_data_clustered, reduction = "umap")
    expr <- FetchData(sc_data_clustered, vars = c(genes_of_interest, housekeeping_genes[1]))
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
    ggsave(filename = paste0(dir_results, "/", file_name,"_DensityPlot.png"),
           plot = density_plot,  
           width = 1920*2, height = 1980*2, units = "px")

    

    # Save
    # Seurat Obj
    saveRDS(sc_data_clustered, file = paste0(dir_results, "/", file_name, "_SeuratObj.rds"))

    # Cluster
    write.csv(sc_data_clustered@meta.data, file = paste0(dir_results, "/", file_name, "_SeuratObj_meatadata.csv"))

    # PLOTS
    for (nm in names(cluster_plots)) {
      ggsave(
        filename = paste0(dir_results, "/", file_name,"_",nm, ".png"),
        plot = cluster_plots[[nm]], 
        width = 1920*2, height = 1980*2, units = "px"
      )
    }

    # Markers
    write.csv(cluster_markers, file = paste0(dir_results, "/", file_name,"_cluster_markers",".csv"))
    
    # Top Genes
    write.csv(top_genes, file = paste0(dir_results, "/", file_name,"_top_genes",".csv"))

    elapsed_sub <- toc(log = TRUE, quiet = T)
    total_time(elapsed_sub$toc - elapsed_sub$tic)
}

elapsed <- toc(log = TRUE, quiet = T)
total_time(elapsed$toc - elapsed$tic)