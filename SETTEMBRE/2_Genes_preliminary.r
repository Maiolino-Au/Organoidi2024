# LOAD LIBRARIES
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tictoc))

total_time <- function(seconds) {
    d <- seconds %/% (86400)
    h <- (seconds %% 86400) %/% 3600
    m <- (seconds %% 3600) %/% 60
    s <- seconds %% 60
    
    cat(sprintf("Total Time: %f Days, %f Hours, %f Minutes and %f Seconds\n", d, h, m, s))
}

timepoints = c("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month")
housekeeping_genes = c("ACTB", "DLG4")
genes_of_interest = c("SRCIN1", "KIAA1217", "CIT")
path_to_data = "/sharedFolder/Data/"
dir_save = "Results_september/"

de.genes <- function(
    file_name,
    goi,
    hk
) {
    cluster_markers <- read.csv(paste0(dir_save, file_name, "/", file_name, "_cluster_markers.csv"))

    de_list <- cluster_markers %>% filter(gene %in% goi)

    for (nm in names(de_list)) {
        if (nrow(de_list[[nm]]) > 0) de_list[[nm]]$type <- as.character(nm)
    }

    de_genes <- bind_rows(de_list)

    return(de_genes)
}

tic()
for (file_name in timepoints) {
    de_genes <- de.genes(
        file_name = file_name,
        goi = genes_of_interest,
        hk = housekeeping_genes
    )
    print(file_name)
    print(de_genes)

    dir_results <- paste0(dir_save, file_name)     
    if (!dir.exists(dir_results)) {dir.create(dir_results)} 
    
    write.csv(de_genes, file = paste0(dir_results, "/", file_name,"_de_genes",".csv"))
}
elapsed <- toc(log = TRUE, quiet = T)
total_time(elapsed$toc - elapsed$tic)