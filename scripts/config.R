# ============================================================
# config.R
# ImmunoVerse MS Proof-of-Concept
# ============================================================

# ------------------------------------------------------------
# Project information
# ------------------------------------------------------------

PROJECT_NAME <- "ImmunoVerse_MS_PoC"
GEO_ACCESSION <- "GSE144744"

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

RAW_DIR <- "data/raw"
PROCESSED_DIR <- "data/processed"

RESULTS_DIR <- "results"
TAB_DIR <- file.path(RESULTS_DIR, "tables")
FIG_DIR <- file.path(RESULTS_DIR, "figures")

MODEL_DIR <- "models"

# ------------------------------------------------------------
# QC parameters
# ------------------------------------------------------------

MIN_FEATURES <- 200
MAX_FEATURES <- 6000

MIN_COUNTS <- 500
MAX_COUNTS <- 50000

MAX_MITO <- 15

# ------------------------------------------------------------
# Clustering parameters
# ------------------------------------------------------------

N_VARIABLE_FEATURES <- 3000

N_PCS <- 30

UMAP_DIMS <- 1:30

CLUSTER_RESOLUTION <- 0.5

# ------------------------------------------------------------
# Foundation model export
# ------------------------------------------------------------

EXPORT_H5AD <- TRUE

# ------------------------------------------------------------
# Create directories if missing
# ------------------------------------------------------------

dirs_to_create <- c(
  RAW_DIR,
  PROCESSED_DIR,
  RESULTS_DIR,
  TAB_DIR,
  FIG_DIR,
  MODEL_DIR
)

for (d in dirs_to_create) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
