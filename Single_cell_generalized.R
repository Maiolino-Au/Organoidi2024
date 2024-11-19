library(Seurat)
library(tidyverse)
library(future)

# Ask for initial conditions for program
# []

set_working_directory <- function() {
  cat("Current directory:", getwd())
  work_dir <- readline(prompt="Enter your working directory: ")
  setwd(work_dir)
}
set_working_directory()

# Ask for number of cores, ram

number_of_cores <- as.numeric(readline(prompt="how many cores do you whant to use"))

# Create folder for results

# Ask for the raw data
data_dir <- readline(prompt="Enter directory of the raw data (sparse matrix")
raw_data <- Read10X(data.dir = data_dir)

# Create Seurat object
data_SO <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = 200, project = "data_SO", names.delim = "-", names.field = 2)

# add the percentage of these mitochondrial or ribosomal genes to the meta.data
data_SO[["percent.mito"]] <- PercentageFeatureSet(object = data_SO, pattern = "^mt-")
data_SO[["percent.ribo"]] <- PercentageFeatureSet(object = data_SO, pattern = "Rps|Rpl|Mrpl|Mrps")

data_SO <- subset(x = data_SO, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40 )


data_SO <- NormalizeData(data_SO, normalization.method = "LogNormalize", scale.factor = 10000)
data_SO <- FindVariableFeatures(data_SO, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data_SO)
plan("multisession", workers = 4)
# Increase the maximum allowed size for objects passed to workers
options(future.globals.maxSize = 6 * 1024^3)  # 6 GB
data_SO <- ScaleData(data_SO, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multisession", workers = 1)
# Restore the maximum allowed size for objects passed to workers
options(future.globals.maxSize = 0.5 * 1024^3)  # 512 MB


data_SO <- RunPCA (data_SO, features = VariableFeatures(object = data_SO), ndims.print = 1:2)


ElbowPlot(object = data_SO, ndims = 50)
number_of_PCs <- as.numeric(readline(prompt="How many PCs do you want to use? "))

data_SO <- FindNeighbors(data_SO, dims = 1:number_of_PCs)
data_SO <- FindClusters(data_SO, resolution = c(0.2, 0.4, 0.6))
data_SO <- RunTSNE(data_SO, dims = 1:number_of_PCs)

r02 <- DimPlot(data_SO, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.2")+labs(title="resolution 0.2")
r04 <- DimPlot(data_SO, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.4")+labs(title="resolution 0.4")
r06 <- DimPlot(data_SO, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.6")+labs(title="resolution 0.6")

CombinePlots(plots=list(r02,r04,r06), ncol = 3)