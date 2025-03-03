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

cores_ram <- function(cores,ram) {
  plan("multisession", workers = cores)
  options(future.globals.maxSize = ram * 1024^3)
}

name_new_dir_0 <- paste(getwd(),"/Results/Cluster&CellTypes",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________



for (y in 1:length(timepoints)) {
  UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[y],".txt",sep=""))
  UMAP_coordinates <- UMAP_coordinates[-1,]
  for (i in 2:7) {if (i!=4) {UMAP_coordinates[,i]<-as.numeric(UMAP_coordinates[,i])}}
  
  CellType <- sort(unique(UMAP_coordinates$CellType))
  N_clusters <- max(as.numeric(UMAP_coordinates$Cluster))
  
  results <- matrix(0, ncol = N_clusters+1, nrow = length(CellType))
  results <- data.frame(results)
  rownames(results) <- CellType
  colnames(results) <- 0:N_clusters
  
  print(timepoints[y])
  for (CT in 1:length(CellType)) {
    Tot <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT],7])
    for (Cl in 0:N_clusters) {
      results[CT,Cl+1] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT] & UMAP_coordinates$Cluster==Cl,7])/Tot
    }
    print(paste(CT,"-",sum(results[CT,])==1))
  }
  write.csv(results, file = paste(name_new_dir_0,"/",timepoints[y],"_Cluster&CellTypes.csv",sep=""))
  
  results <- results %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  
  results$colname <- as.numeric(results$colname)
  
  heatmap_plot <- ggplot(results,aes(colname, rowname, fill=value)) + 
    geom_tile() +
    
    # further customizing the heatmap by
    theme_minimal() + # applying colors and title
    
    scale_fill_gradient(low="white", high="blue") + # setting gradient color as red and white
    
    # setting the title and subtitles using
    # title and subtitle
    labs(title = timepoints[y]) +
    labs(subtitle = NULL) +
    
    labs(x =NULL, y =NULL) + # setting x and y labels using labs
    
    scale_x_continuous(labels = as.character(0:N_clusters), breaks = 0:N_clusters) +
    
    theme(
      title = element_text(size = 30),
      
      axis.text = element_text(size = 20),
      
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      
      legend.key.size = unit(1, 'cm')
    )
  
  jpeg(paste(name_new_dir_0,"/",timepoints[y],"_Cluster&CellTypes.jpg",sep=""),width=1920, height=1080)
  print(heatmap_plot)
  dev.off()
}



# _______________________________________________________________________________



y <- 2
UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",timepoints[y],".txt",sep=""))
UMAP_coordinates <- UMAP_coordinates[-1,]
for (i in 2:7) {if (i!=4) {UMAP_coordinates[,i]<-as.numeric(UMAP_coordinates[,i])}}

CellType <- sort(unique(UMAP_coordinates$CellType))
N_clusters <- max(as.numeric(UMAP_coordinates$Cluster))

results <- matrix(0, ncol = N_clusters+1, nrow = length(CellType))
results <- data.frame(results)
rownames(results) <- CellType
colnames(results) <- 0:N_clusters

print(timepoints[y])
for (CT in 1:length(CellType)) {
  Tot <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT],7])
  for (Cl in 0:N_clusters) {
    results[CT,Cl+1] <- length(UMAP_coordinates[UMAP_coordinates$CellType==CellType[CT] & UMAP_coordinates$Cluster==Cl,7])/Tot
  }
  print(paste(CT,"-",sum(results[CT,])==1))
}


results <- results %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

results$colname <- as.numeric(results$colname)

heatmap_plot <- ggplot(results,aes(colname, rowname, fill=value)) + 
  geom_tile() +
  
  # further customizing the heatmap by
  theme_minimal() + # applying colors and title
  
  scale_fill_gradient(low="white", high="blue", limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) + # setting gradient color as red and white
  
  # setting the title and subtitles using
  # title and subtitle
  labs(title = timepoints[y]) +
  labs(subtitle = NULL) +
  
  labs(x =NULL, y =NULL) + # setting x and y labels using labs
  
  scale_x_continuous(labels = as.character(0:N_clusters), breaks = 0:N_clusters) +
  
  theme(
    title = element_text(size = 30),
    
    axis.text = element_text(size = 20),

    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),

    legend.key.size = unit(1, 'cm')
  )

heatmap_plot

jpeg(paste(name_new_dir_0,"/",timepoints[y],"_Cluster&CellTypes.jpg",sep=""),width=1920, height=1080)
print(heatmap_plot)
dev.off()


