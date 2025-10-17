setwd('/sharedFolder/Paper_clusters_from_scratch/0_September')
source("6_run.r")
tic()
for (i in timepoints) {
    run.month(file_name = i, 
              path_to_data = "/sharedFolder/Data/",
              dir_save = "/sharedFolder/Paper_clusters_from_scratch/0_September/plots_oct/",
             output = T)
}
elapsed <- toc(log = TRUE, quiet = T)
total_time(elapsed$toc - elapsed$tic)