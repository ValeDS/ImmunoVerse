# ============================================================
# 04_export_for_foundation_models.R
# Export Seurat object for Python/scGPT/Geneformer
# ============================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

in_rds <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_seurat_qc_clustered.rds"))
if (!file.exists(in_rds)) stop("Processed Seurat object not found. Run 03_qc_preprocess_clustering.R first.")

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"

export_dir <- file.path(PROCESSED_DIR, "foundation_model_export")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- seu@meta.data %>% as_tibble(rownames = "cell_id")
readr::write_csv(metadata, file.path(export_dir, "cell_metadata.csv"))

# Gene list
readr::write_csv(tibble(gene = rownames(seu)), file.path(export_dir, "genes.csv"))

# Raw counts as sparse Matrix Market
if (EXPORT_RAW_COUNTS) {
  counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
  Matrix::writeMM(counts, file.path(export_dir, "raw_counts.mtx"))
  readr::write_csv(tibble(cell_id = colnames(counts)), file.path(export_dir, "cells.csv"))
}

# Normalized data as sparse Matrix Market
if (EXPORT_NORMALIZED_DATA) {
  norm_data <- GetAssayData(seu, assay = "RNA", slot = "data")
  Matrix::writeMM(norm_data, file.path(export_dir, "normalized_log_data.mtx"))
}

# Optional h5Seurat/h5ad export when SeuratDisk is available.
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  message("SeuratDisk not installed. Skipping h5Seurat/h5ad export.")
  message("To enable it: remotes::install_github('mojaveazure/seurat-disk')")
} else {
  h5seurat_path <- file.path(export_dir, paste0(GEO_ACCESSION, ".h5Seurat"))
  h5ad_path <- file.path(export_dir, paste0(GEO_ACCESSION, ".h5ad"))
  SeuratDisk::SaveH5Seurat(seu, filename = h5seurat_path, overwrite = TRUE)
  SeuratDisk::Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
  message("Saved h5ad export: ", h5ad_path)
}

message("Export completed. Files are in: ", export_dir)
