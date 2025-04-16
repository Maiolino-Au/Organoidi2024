# setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

# install.packages('devtools')
# library(devtools)
# devtools::install_github('immunogenomics/presto')

library(presto)

# install.packages('enrichR')
library(enrichR)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

name_new_dir_0 <- paste(getwd(),"/Results/MARZO",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")


cores_ram(1,64)

# __________________________________________________________________________________________

y <- 2
a <- timepoints[y]
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_CellType",timepoints[y],".Robj", sep = ""))
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_",timepoints[y],"_CellType_markers.Robj", sep = ""))

top100 <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  as.data.frame() %>% 
  arrange(cluster, desc(avg_log2FC))

# __________________________________________________________________________________
# ___________________________________ Annotation ___________________________________
# __________________________________________________________________________________

# ___________________________________ ChatGPT ___________________________________

CellType <- c(unique(top100$cluster))
top_genes_per_cell <- data.frame()

for (i in 1:length(CellType)) {
  top_genes_per_cell[i,1] <- CellType[i]
  top_genes_per_cell[i,2:101] <- c(top100[top100$cluster==CellType[i],7])
}


# for (j in 1:length(CellType)) {
#   testo <- "Hi, this is a list of the 100 most diffrentially expressed genes of a cluster derived from the single cell RNA sequencing of cells taken from a brain organid. Can you give me your 3 best guesses wich cell types are they? Write giust the 3 guesses in decreasing order of probability. "
#   b <- ""
#   if (b!=testo) {b<-testo}
#   for (i in 2:101) {
#     b <- paste(b,top_genes_per_cell[j,i])
#   }
#   print(CellType[j])
#   print(b)
#   print("----------------------------------------------")
# }

name_new_dir_1 <- paste(name_new_dir_0,"/Annotations_ChatGPT",sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 

setwd(name_new_dir_1)
# Define file name
output_file <- paste("Ann_ChatGPT_",timepoints[y],".txt",sep="")

# Open file connection
file_conn <- file(output_file, "w")

for (j in 1:length(CellType)) {
  testo <- "Hi, this is a list of the 100 most diffrentially expressed genes of a cluster derived from the single cell RNA sequencing of cells taken from a brain organid. Can you give me your 3 best guesses wich cell types are they? Write giust the 3 guesses in decreasing order of probability. "
  b <- ""
  if (b != testo) { b <- testo }
  for (i in 2:101) {
    b <- paste(b, top_genes_per_cell[j, i])
  }
  
  # Write to file
  write(as.character(CellType[j]), file = file_conn, append = TRUE)  # Title
  write(b, file = file_conn, append = TRUE)            # Main text
  write("", file = file_conn, append = TRUE)           # Empty line
}

# Close file connection
close(file_conn)

print(paste("File saved as", output_file))

setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")




# ___________________________________ EnirchR ___________________________________

name_new_dir_1 <- paste(name_new_dir_0,"/Annotations_EnrichR",sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 


# List available databases
dbs <- listEnrichrDbs()
# print(dbs)


gene_list <- top100[top100$cluster==CellType[1],]
gene_list <- gene_list[order(gene_list$avg_log2FC, decreasing = TRUE),] 
gene_list <- c(gene_list[,7])



# Perform enrichment analysis
enriched <- enrichr(gene_list, databases = c("PanglaoDB_Augmented_2021", 
                                             "ARCHS4_Tissues", 
                                             "Allen_Brain_Atlas_10x_scRNA_2021",
                                             "KEGG_2021_Human",
                                             "CellMarker_Augmented_2021",
                                             "Azimuth_Cell_Types_2021"))

# View results
# print(enriched$PanglaoDB_Augmented_2021)


# __________________________________________________________________________________
# ___________________________________ Comparison ___________________________________
# __________________________________________________________________________________



# ___________________________________ Violin ___________________________________


# Fetch the expression data for SRCIN1 and Actin
genes_of_interest <- c("SRCIN1", "ACTB")
expression_data <- FetchData(sc_data, vars = genes_of_interest)



# Plot the expression levels of SRCIN1 and Actin across external clusters
violin <- VlnPlot(sc_data, features = genes_of_interest, group.by = "external_cluster")

name_new_dir_1 <- paste(name_new_dir_0,"/Violin",sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 
jpeg(paste(name_new_dir_1,"/Violin_SRCIN1_ACTB",timepoints[y],".jpg",sep=""),width=480*2, height=270*2)
print(violin)
dev.off()

# ___________________________________ Ratio ___________________________________


# Fetch expression data for SRCIN1 and ACTB
genes_of_interest <- c("SRCIN1", "ACTB")
expression_data <- FetchData(sc_data, vars = genes_of_interest)

# Calculate the ratio of SRCIN1 to ACTB for each cell
expression_data$SRCIN1_ACTB_ratio <- expression_data$SRCIN1 / expression_data$ACTB

# Add the cluster information to the data
expression_data$cluster <- sc_data$external_cluster[colnames(sc_data)]

# Now, we have a data frame where each cell has the SRCIN1/ACTB ratio and its associated cluster



library(ggplot2)

# Plot a histogram of the SRCIN1/ACTB ratio for each cluster
ratio <-ggplot(expression_data, aes(x = SRCIN1_ACTB_ratio, fill = cluster)) +
            geom_histogram(binwidth = 0.2, position = "dodge", alpha = 0.7) +
            labs(title = "Histogram of SRCIN1/ACTB Ratio by Cluster",
                 x = "SRCIN1/ACTB Ratio",
                 y = "Number of Cells") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_manual(values = rainbow(length(unique(expression_data$cluster)))) 

name_new_dir_1 <- paste(name_new_dir_0,"/Ratio",sep="")
if (dir.exists(name_new_dir_1)==FALSE) {
  dir.create(name_new_dir_1)
} 
jpeg(paste(name_new_dir_1,"/Ratio_SRCIN1_ACTB",timepoints[y],".jpg",sep=""),width=480*2, height=270*2)
print(ratio)
dev.off()

# # Calculate the mean and standard deviation of ACTB expression for each cluster
# actb_stats <- aggregate(ACTB ~ cluster, data = expression_data, FUN = function(x) c(mean = mean(x), sd = sd(x)))
# 
# # Flatten the data for better visualization
# actb_stats <- do.call(data.frame, actb_stats)
# 
# # Plot the mean and SD of ACTB expression across clusters
# ggplot(actb_stats, aes(x = cluster, y = ACTB.mean)) +
#   geom_bar(stat = "identity", fill = "lightblue", color = "black") +
#   geom_errorbar(aes(ymin = ACTB.mean - ACTB.sd, ymax = ACTB.mean + ACTB.sd), width = 0.2) +
#   labs(title = "Variation of ACTB Expression Across Clusters",
#        x = "Cluster",
#        y = "ACTB Expression (Mean ± SD)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



# _______ 17 marzo

# gene_count <- FetchData(sc_data, vars = genes_of_interest)
# 
# gene_counts <- GetAssayData(sc_data, slot = "counts")[genes_of_interest, ]
# 
# summary(FetchData(sc_data, vars = genes_of_interest))
# summary(GetAssayData(sc_data, slot = "counts")[genes_of_interest, ])
# 
# cluster_cells <- WhichCells(sc_data, ident = "aRG")
# gene_counts_cluster1 <- GetAssayData(sc_data, slot = "counts")[genes_of_interest, cluster_cells]






