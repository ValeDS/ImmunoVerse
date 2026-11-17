# ============================================================
# 01_download_GSE144744.R
# Download GEO supplementary files and metadata
# ============================================================

source("scripts/config.R")
suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

message("Downloading GEO series matrix metadata for ", GEO_ACCESSION, "...")
gse <- GEOquery::getGEO(GEO_ACCESSION, GSEMatrix = TRUE, getGPL = FALSE)

if (length(gse) > 1) {
  message("Multiple platforms detected. Using the first object; inspect manually if needed.")
}

pheno <- pData(gse[[1]])
metadata_path <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_sample_metadata_raw.csv"))
readr::write_csv(as_tibble(pheno, rownames = "sample_id"), metadata_path)
message("Saved metadata: ", metadata_path)

message("Downloading supplementary files. This may take some time depending on GEO size...")
GEOquery::getGEOSuppFiles(GEO_ACCESSION, baseDir = RAW_DIR, makeDirectory = TRUE)

message("Supplementary files downloaded to: ", RAW_DIR)

# Decompress common archive formats
files <- list.files(RAW_DIR, recursive = TRUE, full.names = TRUE)

for (f in files) {
  if (grepl("\\.tar$", f, ignore.case = TRUE)) {
    message("Untarring: ", basename(f))
    utils::untar(f, exdir = dirname(f))
  }
  if (grepl("\\.gz$", f, ignore.case = TRUE) && !grepl("(matrix|features|barcodes).*\\.gz$", basename(f), ignore.case = TRUE)) {
    message("Keeping compressed file as-is: ", basename(f))
  }
}

message("Download step completed.")
