source(file = "6_functions.r")

timepoints = c("1month","2month", "3month")
dir_save = "/sharedFolder/Results/" %>% dir.create(recursive = T, showWarnings = F)

run.month <- function (
    file_name,
    ref = ref,
    genes_of_interest = c("SRCIN1", "SATB2"),
    path_to_data = "/Data/",
    dir_save = dir_save,
    resolution = 1,
    dimensions = c(1:20)
) {
    # Create directory for specific file
    dir_results = paste0(dir_save, "Res_", file_name, "/") %>% dir.create(recursive = T, showWarnings = F)

    # Load data
    raw_data <- load.data(file_name = file_name,
                         data_path = path_to_data)

    # Preprocessing
    sc_data <- preprocessing(data = raw_data,
                            file_name = file_name)

    # Clustering
    sc_data <- PCA.cluster(data, 
                           file_name, 
                           res = 1, 
                           n_dim = 20, 
                           save_dir = dir_save,
                           save_add_on = "",)

    # Markers

    # Differentially expressed genes

    # Plots, no annotation

    # Annotation

    # Subclustering

    # Markers of subclusters

    # Differentially expressed genes

    # Plots
}


ref <- readH5AD("/sharedFolder/public-anndata-project_570_1__published__cKUH0kM97x2oczRNPkqvk_1265.h5ad")

for (i in timepoints) {
    run.month(file_name = i)
}
