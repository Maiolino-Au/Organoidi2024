source(file = "6_functions.r")

timepoints <- c("1month", "2month", "3month")
dir_save <- "/sharedFolder/Results/" 

run.month <- function(
    file_name,
    ref = "",
    genes_of_interest = c("SRCIN1", "SATB2"),
    path_to_data = "/sharedFolder/Data/",
    dir_save = dir_save,
    resolution = 1,
    dimensions = c(1:20),
    output = F
) {
    tic()
    # Create directory for specific file
    dir.create(dir_save, recursive = T, showWarnings = F)
    dir_results <- paste0(dir_save, "Res_", file_name, "/")
    dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)

    # Load data
    raw_data <- load.data(
        file_name = file_name,
        data_path = path_to_data
    )

    print("Preprocessing")

    # Preprocessing
    sc_data <- preprocessing(
        data = raw_data,
        file_name = file_name,
        output = output
    )

    print("Clustering")

    # Clustering
    sc_data <- PCA.cluster(data = sc_data,
        file_name = file_name,
        res = resolution,
        n_dim = dimensions,
        dir_save = dir_save,
        save_add_on = "",
        output = output
    )

    print("Markers")

    # Markers
    cluster_markers <- cluster.markers(data = sc_data,
        output = output)
    write.csv(cluster_markers, file = paste0(dir_results, "/markers_", file_name, ".csv"))

    print("de_genes")
    
    # Differentially expressed genes
    de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest) %>% arrange(avg_log2FC) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.6) 
    write.csv(de_genes, file = paste0(dir_results, "/de_genes_", file_name, ".csv"))

    print("Plots")
    
    # Plots, no annotation
    plottamelo.tutto(
        data = sc_data,
        file_name = file_name,
        genes_of_interest = genes_of_interest,
        dir_save = dir_results,
        save_add_on = "no_annotation",
        pt.size = 1,
        print_plot = F
    )

    # Annotation
    sc_data <- annotate.singler(
        data = sc_data, 
        file_name = file_name,
        ref = ref,
        dir_save = dir_save,
        output = F
    )

    plottamelo.tutto(
        data = sc_data,
        file_name = file_name,
        genes_of_interest = genes_of_interest,
        dir_save = dir_results,
        save_add_on = "SingleR_annotation",
        pt.size = 1,
        print_plot = F
    )

    elapsed <- toc(log = TRUE, quiet = T)
    return(total_time(elapsed$toc - elapsed$tic))
}

#library(reticulate)
#use_python("/usr/bin/python3", required = TRUE)
#library(zellkonverter)

#ref <- LoadH5Seurat("/SingleR_reference/HNOA_cleaned_seurat.rds")

library(zellkonverter)
library(SingleCellExperiment)
library(scuttle)
library(SingleR)
ref <- readH5AD("/SingleR_reference/HNOA_cleaned.h5ad")

# Convert to SingleCellExperiment
ref <- as.SingleCellExperiment(ref)

# Perform log-normalization (required by SingleR)
ref <- logNormCounts(ref)

for (i in timepoints) {
    run.month(
        file_name = i, 
        dir_save = "/sharedFolder/plots_oct/",
        ref = ref,
        genes_of_interest = "SRCIN1",
        path_to_data = "/sharedFolder/Data/",
    )

    cat("\n=============== Done -", i, "===============\n\n")
}

print("Done")