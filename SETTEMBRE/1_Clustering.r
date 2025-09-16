source(file = "1_Clustering_Functions.r")

timepoints = c("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month")
housekeeping_genes = c("ACTB", "DLG4")
genes_of_interest = c("SRCIN1", "KIAA1217", "CIT")
path_to_data = "/sharedFolder/Data/"
dir_save = "Results_september/"
res_list = c(1, 0.5, 0.2, 0.1)#seq(0.1, 1, by = 0.1)
dim_list = c(1:20)
top_n_genes = 50

tic("Clustering")

# Seurat
for (file_name in timepoints) {
    tic(file_name)

    dir_results <- paste0(dir_save, file_name)     
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

    cluster_markers <- cluster.markers(sc_data_clustered, output = F)

    top_genes <- top.genes(top_genes_by_cluster, num = 50)

    # Save
    # Seurat Obj
    saveRDS(sc_data_clustered, file = paste0(dir_results, file_name, "_SeuratObj.rds"))

    # Cluster
    write.csv(sc_data_clustered@meta.data, file = paste0(dir_results, file_name, "_SeuratObj_meatadata.csv"))

    # PLOTS
    for (nm in names(cluster_plots)) {
      ggsave(
        filename = paste0(dir_results, "/", file_name,"_",nm, ".png"),
        plot = cluster_plots[[nm]]
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