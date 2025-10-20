#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SeuratDisk)
    library(Seurat)
})

# Get arguments: input .h5ad and output .rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
    stop("Usage: Rscript convert.r <input.h5seurat>\n Output: <output.rds>")
}

if (grepl("\\.h5ad$", args)) {

}

# Create temporary h5Seurat filename
h5seurat_file <- sub("\\.h5ad$", ".h5seurat", input_h5ad)

# Step 2: Load into Seurat
cat("Loading", h5seurat_file, "into Seurat...\n")
obj <- LoadH5Seurat(h5seurat_file)

# Step 3: Save as .RDS
cat("Saving Seurat object to", output_rds, "...\n")
saveRDS(obj, file = output_rds)

cat("✅ Conversion complete!\n")
