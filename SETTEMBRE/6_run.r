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
    
}


ref <- 

for (i in timepoints) {
    run.month(file_name = i)
}