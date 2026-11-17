# ============================================================
# 07_transfer_states_to_full_tcells.R
# Transfer the 36 reference states to the complete T-cell set
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(123)

reference_file <- file.path(
  "data", "processed",
  "GSE144744_tcells_reference_slim.rds"
)

full_tcells_file <- file.path(
  "data", "processed",
  "GSE144744_tcells_raw.rds"
)

chunk_dir <- file.path(
  "data", "processed",
  "state_transfer_chunks"
)

dir.create(
  chunk_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load objects only if not already present
# ------------------------------------------------------------

if (!exists("reference")) {
  reference <- readRDS(reference_file)
}

if (!exists("full_tcells")) {
  full_tcells <- readRDS(full_tcells_file)
}

DefaultAssay(reference) <- "RNA"
DefaultAssay(full_tcells) <- "RNA"

stopifnot(
  all(colnames(reference) %in% colnames(full_tcells)),
  "immune_state_id" %in% colnames(reference@meta.data)
)

# ------------------------------------------------------------
# Mapping features
# ------------------------------------------------------------

features_use <- VariableFeatures(reference)

if (length(features_use) < 500) {
  features_use <- rownames(
    Loadings(reference[["pca"]])
  )
}

features_use <- intersect(
  features_use,
  rownames(full_tcells)
)

stopifnot(length(features_use) >= 500)

message(
  "Mapping with ",
  length(features_use),
  " reference features."
)

# ------------------------------------------------------------
# Identify cells requiring label transfer
# ------------------------------------------------------------

query_cells <- setdiff(
  colnames(full_tcells),
  colnames(reference)
)

query_metadata <- full_tcells@meta.data[
  query_cells,
  ,
  drop = FALSE
]

if (anyNA(query_metadata$sample)) {
  stop("Missing sample identifiers among query cells.")
}

sample_ids <- sort(
  unique(as.character(query_metadata$sample))
)

ref_labels <- setNames(
  as.character(reference$immune_state_id),
  colnames(reference)
)

message(
  length(query_cells),
  " cells will be transferred across ",
  length(sample_ids),
  " samples."
)

# ------------------------------------------------------------
# Transfer sample by sample
# ------------------------------------------------------------

for (i in seq_along(sample_ids)) {
  
  sample_id <- sample_ids[i]
  
  safe_sample_id <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    sample_id
  )
  
  output_file <- file.path(
    chunk_dir,
    paste0(
      sprintf("%03d", i),
      "_",
      safe_sample_id,
      "_state_predictions.rds"
    )
  )
  
  if (file.exists(output_file)) {
    message(
      "[", i, "/", length(sample_ids), "] ",
      sample_id,
      " already completed; skipping."
    )
    next
  }
  
  cells_i <- query_cells[
    query_metadata$sample == sample_id
  ]
  
  message(
    "[", i, "/", length(sample_ids), "] ",
    "Mapping ", sample_id,
    " (", length(cells_i), " cells)"
  )
  
  set.seed(123 + i)
  
  query_i <- subset(
    full_tcells,
    cells = cells_i
  )
  
  DefaultAssay(query_i) <- "RNA"
  
  if (!"data" %in% Layers(query_i[["RNA"]])) {
    query_i <- NormalizeData(
      query_i,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
  
  query_i <- DietSeurat(
    object = query_i,
    assays = "RNA",
    layers = "data",
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )
  
  anchors_i <- FindTransferAnchors(
    reference = reference,
    query = query_i,
    reference.assay = "RNA",
    query.assay = "RNA",
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    reference.reduction = "pca",
    features = features_use,
    dims = 1:30,
    k.filter = NA,
    verbose = FALSE
  )
  
  predictions_i <- TransferData(
    anchorset = anchors_i,
    refdata = ref_labels,
    dims = 1:30,
    verbose = FALSE
  )
  
  predictions_i$cell_id <- rownames(
    predictions_i
  )
  
  predictions_i$sample <- sample_id
  
  saveRDS(
    predictions_i,
    output_file
  )
  
  message(
    "Completed: ",
    sample_id,
    "; median score = ",
    round(
      median(
        predictions_i$prediction.score.max,
        na.rm = TRUE
      ),
      3
    )
  )
  
  rm(
    query_i,
    anchors_i,
    predictions_i
  )
  
  gc(verbose = FALSE)
}

message("All available samples have been processed.")