# ============================================================
# 00_setup.R
# ImmunoVerse MS Proof-of-Concept
# ============================================================

message("Setting up ImmunoVerse MS PoC environment...")
source("scripts/config.R")

required_cran <- c(
  "tidyverse", "data.table", "Matrix", "patchwork", "remotes"
)

required_bioc <- c(
  "GEOquery", "SingleCellExperiment", "zellkonverter", "scater", "scran"
)

required_github <- character(0)

install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in required_cran) install_if_missing_cran(pkg)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat", repos = "https://cloud.r-project.org")
}

folders <- c(RAW_DIR, PROCESSED_DIR, RESULTS_DIR, FIG_DIR, TAB_DIR)
invisible(lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Setup completed.")
