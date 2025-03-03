setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

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

name_new_dir_0 <- paste(getwd(),"/Results/Punto2",sep="")
if (dir.exists(name_new_dir_0)==FALSE) {
  dir.create(name_new_dir_0)
} 

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")
genes_of_interest <- c("SRCIN1")

cores_ram(1,30)


# _______________________________________________________________________________


a <- read.csv("D:/Organoidi_2024/Results/Cluster_from_paper/Named_CellType_23days/Named_23days_top_30_markers_CellType.csv")
b<-c(a[a$cluster=="aRG",8])
save(b,file = paste(getwd(),"/top_genes_23d_aRG.csv",sep=""))
