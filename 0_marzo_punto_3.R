setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

# setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)

# install.packages('devtools')
# library(devtools)
# devtools::install_github('immunogenomics/presto')

library(presto)

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

name_new_dir_0 <- paste(getwd(),"/Results/Punto3",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________

sc_data <- Read10X(data.dir = path_to_data, gene.column = 1) 
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = timepoints[i], names.delim = "-", names.field = 2)
# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6)
# Identify variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000)


UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[i],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (j in 2:7) {if (j!=4) {UMAP_coordinates[,j]<-as.numeric(UMAP_coordinates[,j])}}

# Create a named vector of clusters
cluster_map <- setNames(UMAP_coordinates$CellType, UMAP_coordinates$NAME)

# Assign external clusters to the Seurat object
sc_data$external_cluster <- cluster_map[colnames(sc_data)]

# Set 'external_cluster' as the active identity
Idents(sc_data) <- sc_data$external_cluster

# Ensure 'external_cluster' is a factor with ordered levels
sc_data$external_cluster <- factor(sc_data$external_cluster, levels = sort(unique(sc_data$external_cluster)))

# Set active identity using the ordered levels
Idents(sc_data) <- sc_data$external_cluster

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
