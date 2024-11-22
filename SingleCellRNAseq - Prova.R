# getwd()
# setwd("/home/user/Documents/organoidi_2024/Data - 3 - Organoidi_velasco")

cat("Current directory:", getwd())
i <- readline(prompt="Change directory? [y/n] ")
if (i=="y") {
  work_dir <- readline(prompt="Enter your working directory: ")
  setwd(work_dir)
  cat("Current directory:", getwd())
}


library(Seurat)
library(tidyverse)
library(future)


# Read in data
# For some reasons it wants barcodes.tsv.gz ; features.tsv.gz ; matrix.mtx.gz (in windows that was not necessary)
raw_data_days23 <- Read10X(data.dir = "Data/expression_23days", gene.column = 1)


# Create Seurat object
days23 <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = 200, project = "days23", names.delim = "-", names.field = 2)


# # add the percentage of these mitochondrial or ribosomal genes to the meta.data
# days23[["percent.mito"]] <- PercentageFeatureSet(object = days23, pattern = "^mt-")
# days23[["percent.ribo"]] <- PercentageFeatureSet(object = days23, pattern = "Rps|Rpl|Mrpl|Mrps")
# 
# plot1 <- FeatureScatter(object = days23, feature1 = "nCount_RNA", feature2 = "percent.mito")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
# plot2 <- FeatureScatter(object = days23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
# plot3 <- FeatureScatter(object = days23, feature1 = "nCount_RNA", feature2 = "percent.ribo")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)


days23 <- subset(x = days23, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40 )


days23 <- NormalizeData(days23, normalization.method = "LogNormalize", scale.factor = 10000)
days23 <- FindVariableFeatures(days23, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(days23)
# plan("multisession", workers = 2)
# Increase the maximum allowed size for objects passed to workers
options(future.globals.maxSize = 24 * 1024^3)  # 24 GB
days23 <- ScaleData(days23, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
# plan("multisession", workers = 1)
# Restore the maximum allowed size for objects passed to workers
options(future.globals.maxSize = 0.5 * 1024^3)  # 512 MB


days23 <- RunPCA (days23, features = VariableFeatures(object = days23), ndims.print = 1:2)


ElbowPlot(object = days23, ndims = 50)

n_dims <- as.numeric(readline(ptompt="Number of dimension for clustering: "))

days23 <- FindNeighbors(days23, dims = 1:40)
days23 <- FindClusters(days23, resolution = c(0.2, 0.4, 0.6))
days23 <- RunTSNE(days23, dims = 1:40)

# save the processed object now in case you want to come back to it later
save(days23, file = "days23_all_v3.Robj")


r02 <- DimPlot(days23, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.2")+labs(title="resolution 0.2")
r04 <- DimPlot(days23, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.4")+labs(title="resolution 0.4")
r06 <- DimPlot(days23, label = T, reduction = "tsne", group.by = "RNA_snn_res.0.6")+labs(title="resolution 0.6")

CombinePlots(plots=list(r02,r04,r06), ncol = 3)

# To see specifically cells that express the SRCIN1
FeaturePlot(days23, c("SRCIN1"))



# plan("multisession", workers = 4)
options(future.globals.maxSize = 24 * 1024^3)  # 24 GB
days23.markers <- FindAllMarkers(days23, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plan("multisession", workers = 1)
options(future.globals.maxSize = 0.5 * 1024^3)  # 512 MB

# save this to file in case you want to come back to it later
save(days23.markers, file = "days23_all_markers_res06_v3.Robj")

days23.markers[days23.markers$gene=="SRCIN1"]

# top_markers_all <- days23.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# # write.csv(top_markers_all, file = "days23_all_markers_res06_v3.csv")
# top_markers_all


# find_gene <- function(genes_of_interest) {
#   FeaturePlot(days23, c(genes_of_interest))
#   days23.markers[days23.markers$gene==genes_of_interest]
# }

