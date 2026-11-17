# ============================================================
# 02_create_seurat_object.R
# Create Seurat object from GSE144744 RNA counts + per-cell metadata
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

message("Creating Seurat object for ", GEO_ACCESSION)

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

gse_raw_dir <- file.path(RAW_DIR, GEO_ACCESSION)
gse_processed_dir <- file.path(PROCESSED_DIR, GEO_ACCESSION)

dir.create(gse_processed_dir, recursive = TRUE, showWarnings = FALSE)

rna_tar <- file.path(gse_raw_dir, paste0(GEO_ACCESSION, "_RNA_counts.tar.gz"))
meta_cell_file <- file.path(gse_raw_dir, paste0(GEO_ACCESSION, "_metadata_per_cell.csv.gz"))
meta_sample_file <- file.path(gse_raw_dir, paste0(GEO_ACCESSION, "_metadata_per_sample.csv.gz"))

if (!file.exists(rna_tar)) {
  stop("RNA counts archive not found: ", rna_tar)
}

if (!file.exists(meta_cell_file)) {
  stop("Per-cell metadata file not found: ", meta_cell_file)
}

# ------------------------------------------------------------
# 2. Extract RNA counts archive if needed
# ------------------------------------------------------------

rna_counts_dir <- file.path(gse_processed_dir, "RNA_counts")

if (!dir.exists(rna_counts_dir)) {
  message("Extracting RNA counts archive...")
  untar(rna_tar, exdir = gse_processed_dir)
} else {
  message("RNA counts directory already exists. Skipping extraction.")
}

expected_files <- file.path(
  rna_counts_dir,
  c("matrix.mtx", "barcodes.tsv", "genes.tsv")
)

if (!all(file.exists(expected_files))) {
  stop(
    "Expected RNA count files not found. Expected:\n",
    paste(expected_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 3. Read RNA count matrix
# ------------------------------------------------------------

message("Reading RNA count matrix...")

rna_counts <- ReadMtx(
  mtx = file.path(rna_counts_dir, "matrix.mtx"),
  cells = file.path(rna_counts_dir, "barcodes.tsv"),
  features = file.path(rna_counts_dir, "genes.tsv"),
  feature.column = 1
)

message("RNA matrix dimensions: ", nrow(rna_counts), " genes x ", ncol(rna_counts), " cells")

# ------------------------------------------------------------
# 4. Read metadata
# ------------------------------------------------------------

message("Reading per-cell metadata...")

meta_cell <- readr::read_csv(
  meta_cell_file,
  show_col_types = FALSE
)

if (!"cell_names" %in% colnames(meta_cell)) {
  stop("Column 'cell_names' not found in per-cell metadata.")
}

meta_cell <- as.data.frame(meta_cell)
rownames(meta_cell) <- meta_cell$cell_names

# ------------------------------------------------------------
# 5. Match cells between matrix and metadata
# ------------------------------------------------------------

common_cells <- intersect(colnames(rna_counts), rownames(meta_cell))

message("Cells in RNA matrix: ", ncol(rna_counts))
message("Cells in metadata: ", nrow(meta_cell))
message("Common cells: ", length(common_cells))

if (length(common_cells) == 0) {
  stop("No matching cell barcodes between RNA matrix and metadata.")
}

rna_counts <- rna_counts[, common_cells]
meta_cell <- meta_cell[common_cells, , drop = FALSE]

# ------------------------------------------------------------
# 6. Create Seurat object
# ------------------------------------------------------------

message("Creating Seurat object...")

seu <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = meta_cell,
  project = PROJECT_NAME,
  min.cells = 3,
  min.features = MIN_FEATURES
)

seu$dataset <- GEO_ACCESSION

# ------------------------------------------------------------
# 7. Save object and summaries
# ------------------------------------------------------------

out_rds <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_seurat_raw.rds"))
saveRDS(seu, out_rds)

message("Saved raw Seurat object: ", out_rds)

summary_tbl <- tibble(
  n_cells = ncol(seu),
  n_genes = nrow(seu),
  assays = paste(names(seu@assays), collapse = ";"),
  groups = paste(names(table(seu$group)), collapse = ";"),
  cohorts = paste(names(table(seu$cohort)), collapse = ";"),
  basictypes = paste(names(table(seu$basictype)), collapse = ";")
)

readr::write_csv(
  summary_tbl,
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_seurat_raw_summary.csv"))
)

# Useful quick checks
message("Group distribution:")
print(table(seu$group))

message("Cohort distribution:")
print(table(seu$cohort))

message("Basic cell type distribution:")
print(table(seu$basictype))

message("Create Seurat object step completed.")
