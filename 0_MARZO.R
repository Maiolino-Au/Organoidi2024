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

y <- 1
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_CellType",timepoints[y],".Robj", sep = ""))
load(paste(getwd(),"/Results/Cluster_from_paper/Named_CellType_",timepoints[y],"/Named_",timepoints[y],"_CellType_markers.Robj", sep = ""))

top100 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()

# __________________________________________________________________________________
# ___________________________________ Annotation ___________________________________
# __________________________________________________________________________________

# ___________________________________ ChatGPT ___________________________________

CellType <- unique(top100$cluster)
top_genes_per_cell <- data.frame()

for (i in 1:length(CellType)) {
  top_genes_per_cell[i,1] <- CellType[i]
  top_genes_per_cell[i,2:101] <- c(top100[top100$cluster==CellType[i],7])
}


for (j in 1:length(CellType)) {
  testo <- "Hi, this is a list of the 100 most diffrentially expressed genes of a cluster derived from the single cell RNA sequencing of cells taken from a brain organid. Can you give me your 3 best guesses wich cell types are they? Write giust the 3 guesses in decreasing order of probability. "
  b <- ""
  if (b!=testo) {b<-testo}
  for (i in 2:101) {
    b <- paste(b,top_genes_per_cell[j,i])
  }
  print(CellType[j])
  print(b)
  print("----------------------------------------------")
}


# ___________________________________ EnirchR ___________________________________

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
print(enriched$PanglaoDB_Augmented_2021)

# __________________________________________________________________________________
# ___________________________________ Comparison ___________________________________
# __________________________________________________________________________________



# ______________


# Fetch the expression data for SRCIN1 and Actin
genes_of_interest <- c("SRCIN1", "ACTA1")
expression_data <- FetchData(sc_data, vars = genes_of_interest)



# Plot the expression levels of SRCIN1 and Actin across external clusters
VlnPlot(sc_data, features = genes_of_interest, group.by = "external_cluster")



# Perform ANOVA for each gene (SRCIN1 and Actin) to compare their expression levels across clusters
anova_results <- lapply(genes_of_interest, function(gene) {
  # Perform ANOVA
  aov_result <- aov(expression_data[[gene]] ~ sc_data$external_cluster)
  summary(aov_result)
})

# View the results
anova_results



# Fetch normalized data
expression_data_norm <- FetchData(sc_data, vars = genes_of_interest, layer = "data")

# Plot with normalized data
VlnPlot(sc_data, features = genes_of_interest, group.by = "external_cluster", layer = "data")


# _____________


# Fetch expression data for SRCIN1 and ACTA1
genes_of_interest <- c("SRCIN1", "ACTA1")
expression_data <- FetchData(sc_data, vars = genes_of_interest)

# Calculate the ratio of SRCIN1 to ACTA1 for each cell
expression_data$SRCIN1_ACTA1_ratio <- expression_data$SRCIN1 / expression_data$ACTA1

# Add the cluster information to the data
expression_data$cluster <- sc_data$external_cluster[colnames(sc_data)]

# Now, we have a data frame where each cell has the SRCIN1/ACTA1 ratio and its associated cluster



library(ggplot2)

# Plot a histogram of the SRCIN1/ACTA1 ratio for each cluster
ggplot(expression_data, aes(x = SRCIN1_ACTA1_ratio, fill = cluster)) +
  geom_histogram(binwidth = 0.2, position = "dodge", alpha = 0.7) +
  labs(title = "Histogram of SRCIN1/ACTA1 Ratio by Cluster",
       x = "SRCIN1/ACTA1 Ratio",
       y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = rainbow(length(unique(expression_data$cluster)))) 



# Calculate the mean and standard deviation of ACTA1 expression for each cluster
acta1_stats <- aggregate(ACTA1 ~ cluster, data = expression_data, FUN = function(x) c(mean = mean(x), sd = sd(x)))

# Flatten the data for better visualization
acta1_stats <- do.call(data.frame, acta1_stats)

# Plot the mean and SD of ACTA1 expression across clusters
ggplot(acta1_stats, aes(x = cluster, y = ACTA1.mean)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  geom_errorbar(aes(ymin = ACTA1.mean - ACTA1.sd, ymax = ACTA1.mean + ACTA1.sd), width = 0.2) +
  labs(title = "Variation of ACTA1 Expression Across Clusters",
       x = "Cluster",
       y = "ACTA1 Expression (Mean ± SD)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# _______ 17 marzo

gene_count <- FetchData(sc_data, vars = genes_of_interest)

gene_counts <- GetAssayData(sc_data, slot = "counts")[genes_of_interest, ]

summary(FetchData(sc_data, vars = genes_of_interest))
summary(GetAssayData(sc_data, slot = "counts")[genes_of_interest, ])

cluster_cells <- WhichCells(sc_data, ident = "aRG")
gene_counts_cluster1 <- GetAssayData(sc_data, slot = "counts")[genes_of_interest, cluster_cells]






