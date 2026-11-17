# ============================================================
# 33_transfer_curated_states_to_GSE138266.R
# External transfer of the 11 curated ImmunoVerse T-cell states
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

set.seed(123)

N_TRANSFER_DIMS <- 30L
DEFAULT_K_ANCHOR <- 5L
DEFAULT_K_SCORE <- 30L
DEFAULT_K_WEIGHT <- 50L
MIN_QUERY_CELLS <- 3L

REFERENCE_FILE <- paste0(
  "data/processed/",
  "GSE144744_Tcells_V3_Harmony_reference_curated_atlas.rds"
)

QUERY_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_Tcells_external_validation.rds"
)

OUTPUT_DIR <- paste0(
  "data/processed/external_validation/GSE138266"
)
TABLE_DIR <- paste0(
  "results/tables/external_validation/GSE138266/",
  "state_transfer"
)

CHUNK_DIR <- file.path(
  OUTPUT_DIR,
  "curated_state_transfer_chunks"
)

PREDICTION_FILE <- file.path(
  OUTPUT_DIR,
  "GSE138266_curated_state_predictions.rds"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(REFERENCE_FILE, QUERY_FILE)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Loading the curated GSE144744 V3 reference...")
reference <- readRDS(REFERENCE_FILE)

message("Loading GSE138266 external-validation T cells...")
query <- readRDS(QUERY_FILE)

if (!inherits(reference, "Seurat") || !inherits(query, "Seurat")) {
  stop("Both input files must contain Seurat objects.")
}

if (!"RNA" %in% Assays(reference) || !"RNA" %in% Assays(query)) {
  stop("The RNA assay is missing from one or both objects.")
}

DefaultAssay(reference) <- "RNA"
DefaultAssay(query) <- "RNA"

required_reference_metadata <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family"
)

missing_reference_metadata <- setdiff(
  required_reference_metadata,
  colnames(reference[[]])
)

if (length(missing_reference_metadata) > 0L) {
  stop(
    "The reference is missing metadata column(s): ",
    paste(missing_reference_metadata, collapse = ", ")
  )
}

required_query_metadata <- c(
  "sample",
  "donor",
  "group",
  "compartment",
  "dataset",
  "chemistry"
)

missing_query_metadata <- setdiff(
  required_query_metadata,
  colnames(query[[]])
)

if (length(missing_query_metadata) > 0L) {
  stop(
    "The query is missing metadata column(s): ",
    paste(missing_query_metadata, collapse = ", ")
  )
}

if (!"pca" %in% Reductions(reference)) {
  stop("The curated reference does not contain a PCA reduction.")
}

if (ncol(Embeddings(reference, reduction = "pca")) < N_TRANSFER_DIMS) {
  stop("The curated reference contains too few PCA dimensions.")
}

# ------------------------------------------------------------
# Reference labels, lookup table, and mapping features
# ------------------------------------------------------------

expected_state_ids <- sort(
  unique(as.character(reference$immune_state_id))
)

if (length(expected_state_ids) != 11L || anyNA(expected_state_ids)) {
  stop(
    "Expected 11 non-missing curated immune states but found ",
    length(expected_state_ids),
    "."
  )
}

reference_state_labels <- setNames(
  as.character(reference$immune_state_id),
  colnames(reference)
)

state_lookup <- reference[[]] %>%
  transmute(
    immune_state_id = as.character(immune_state_id),
    immune_state_name = as.character(immune_state_name),
    state_family = as.character(state_family)
  ) %>%
  distinct() %>%
  arrange(immune_state_id)

if (
  nrow(state_lookup) != 11L ||
    anyDuplicated(state_lookup$immune_state_id)
) {
  stop("The reference state lookup table is not one-to-one.")
}

features_use <- VariableFeatures(reference)

if (length(features_use) < 500L) {
  features_use <- rownames(Loadings(reference[["pca"]]))
}

features_use <- Reduce(
  intersect,
  list(
    features_use,
    rownames(reference),
    rownames(query)
  )
)

features_use <- unique(
  features_use[
    !grepl("^(MT-|RPL|RPS|HBA|HBB|HBD)", features_use)
  ]
)

if (length(features_use) < 500L) {
  stop(
    "Fewer than 500 suitable cross-study mapping features remain."
  )
}

message(
  "External mapping will use ",
  length(features_use),
  " shared reference features."
)

# ------------------------------------------------------------
# Query design and metadata audit
# ------------------------------------------------------------

query_cells <- colnames(query)
query_metadata <- query[[]][
  query_cells,
  required_query_metadata,
  drop = FALSE
]

if (
  anyNA(query_metadata$sample) ||
    any(as.character(query_metadata$sample) == "")
) {
  stop("Missing sample identifiers were found in the query.")
}

sample_ids <- sort(unique(as.character(query_metadata$sample)))

sample_metadata_check <- query_metadata %>%
  mutate(sample_check = as.character(sample)) %>%
  group_by(sample_check) %>%
  summarise(
    across(
      all_of(c("donor", "group", "compartment", "dataset", "chemistry")),
      ~ n_distinct(.x, na.rm = FALSE)
    ),
    .groups = "drop"
  )

inconsistent_samples <- sample_metadata_check %>%
  filter(if_any(-sample_check, ~ .x != 1L))

if (nrow(inconsistent_samples) > 0L) {
  stop(
    "Inconsistent metadata were found within sample(s): ",
    paste(inconsistent_samples$sample_check, collapse = ", ")
  )
}

query_manifest <- query_metadata %>%
  mutate(sample_manifest = as.character(sample)) %>%
  count(
    sample_manifest,
    donor,
    group,
    compartment,
    dataset,
    chemistry,
    name = "query_cells"
  ) %>%
  rename(sample = sample_manifest) %>%
  arrange(compartment, group, sample)

if (length(sample_ids) != 22L || nrow(query_manifest) != 22L) {
  stop(
    "Expected 22 external-validation samples but found ",
    length(sample_ids),
    "."
  )
}

design_check <- query_manifest %>%
  count(compartment, group, name = "n_samples") %>%
  arrange(compartment, group)

expected_design <- data.frame(
  compartment = c("CSF", "CSF", "PBMC", "PBMC"),
  group = c("Control", "MS", "Control", "MS"),
  n_samples = c(6L, 6L, 5L, 5L),
  stringsAsFactors = FALSE
)

if (!identical(
  as.data.frame(design_check),
  as.data.frame(expected_design)
)) {
  stop("The external-validation design failed validation.")
}

write.csv(
  query_manifest,
  file.path(TABLE_DIR, "GSE138266_state_transfer_manifest.csv"),
  row.names = FALSE
)

write.csv(
  state_lookup,
  file.path(TABLE_DIR, "GSE138266_reference_state_lookup.csv"),
  row.names = FALSE
)

transfer_configuration <- data.frame(
  reference_dataset = "GSE144744",
  query_dataset = "GSE138266",
  reference_cells = ncol(reference),
  query_cells = ncol(query),
  query_samples = length(sample_ids),
  query_donors = n_distinct(query_metadata$donor),
  mapping_features = length(features_use),
  transfer_dimensions = N_TRANSFER_DIMS,
  reference_states = length(expected_state_ids),
  default_k_anchor = DEFAULT_K_ANCHOR,
  default_k_score = DEFAULT_K_SCORE,
  default_k_weight = DEFAULT_K_WEIGHT,
  stringsAsFactors = FALSE
)

write.csv(
  transfer_configuration,
  file.path(TABLE_DIR, "GSE138266_state_transfer_configuration.csv"),
  row.names = FALSE
)

message(
  length(query_cells),
  " external T cells will be mapped across ",
  length(sample_ids),
  " samples."
)

# ------------------------------------------------------------
# Reusable chunk validation
# ------------------------------------------------------------

validate_prediction_chunk <- function(
    prediction_object,
    expected_cells,
    expected_states) {

  required_prediction_columns <- c(
    "cell_id",
    "predicted.id",
    "prediction.score.max",
    "sample",
    "donor",
    "group",
    "compartment",
    "dataset",
    "chemistry"
  )

  if (
    is.null(prediction_object) ||
      !is.data.frame(prediction_object) ||
      !all(required_prediction_columns %in% colnames(prediction_object))
  ) {
    return(FALSE)
  }

  if (
    nrow(prediction_object) != length(expected_cells) ||
      anyDuplicated(prediction_object$cell_id) ||
      !setequal(prediction_object$cell_id, expected_cells)
  ) {
    return(FALSE)
  }

  if (
    anyNA(prediction_object$predicted.id) ||
      anyNA(prediction_object$prediction.score.max) ||
      any(!is.finite(prediction_object$prediction.score.max))
  ) {
    return(FALSE)
  }

  score_tolerance <- sqrt(.Machine$double.eps)

  if (
    any(
      prediction_object$prediction.score.max < -score_tolerance |
        prediction_object$prediction.score.max > 1 + score_tolerance
    )
  ) {
    return(FALSE)
  }

  if (!all(
    as.character(prediction_object$predicted.id) %in% expected_states
  )) {
    return(FALSE)
  }

  TRUE
}

# ------------------------------------------------------------
# Transfer labels one external sample at a time
# ------------------------------------------------------------

for (i in seq_along(sample_ids)) {

  sample_id <- sample_ids[[i]]

  cells_i <- query_cells[
    as.character(query_metadata$sample) == sample_id
  ]

  if (length(cells_i) < MIN_QUERY_CELLS) {
    stop(
      "Sample ", sample_id, " contains only ", length(cells_i),
      " cells; at least ", MIN_QUERY_CELLS, " are required."
    )
  }

  safe_sample_id <- gsub("[^A-Za-z0-9_-]", "_", sample_id)

  output_file <- file.path(
    CHUNK_DIR,
    paste0(
      sprintf("%03d", i), "_", safe_sample_id,
      "_external_state_predictions.rds"
    )
  )

  temporary_output_file <- paste0(output_file, ".tmp")

  if (file.exists(output_file)) {
    existing_prediction <- tryCatch(
      readRDS(output_file),
      error = function(e) NULL
    )

    existing_valid <- validate_prediction_chunk(
      existing_prediction,
      expected_cells = cells_i,
      expected_states = expected_state_ids
    )

    rm(existing_prediction)

    if (existing_valid) {
      message(
        "[", i, "/", length(sample_ids), "] ", sample_id,
        " already completed and validated; skipping."
      )
      next
    }

    stop(
      "An invalid or incomplete chunk already exists: ",
      output_file,
      ". Move or inspect it before restarting."
    )
  }

  if (file.exists(temporary_output_file)) {
    unlink(temporary_output_file)
  }

  message(
    "[", i, "/", length(sample_ids), "] Mapping ",
    sample_id, " (", length(cells_i), " cells)"
  )

  set.seed(2000 + i)

  query_i <- subset(query, cells = cells_i)
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
    query_i,
    assays = "RNA",
    layers = "data",
    dimreducs = character(0),
    graphs = character(0),
    misc = FALSE
  )

  k_anchor_i <- as.integer(min(
    DEFAULT_K_ANCHOR,
    ncol(query_i) - 1L,
    ncol(reference) - 1L
  ))

  k_score_i <- as.integer(min(
    DEFAULT_K_SCORE,
    ncol(query_i) - 1L,
    ncol(reference) - 1L
  ))

  if (k_anchor_i < 1L || k_score_i < 1L) {
    stop("Unable to define valid anchor parameters for ", sample_id, ".")
  }

  message(
    "Using k.anchor = ", k_anchor_i,
    " and k.score = ", k_score_i, "."
  )

  anchors_i <- tryCatch(
    FindTransferAnchors(
      reference = reference,
      query = query_i,
      reference.assay = "RNA",
      query.assay = "RNA",
      normalization.method = "LogNormalize",
      reduction = "pcaproject",
      reference.reduction = "pca",
      features = features_use,
      dims = seq_len(N_TRANSFER_DIMS),
      k.anchor = k_anchor_i,
      k.filter = NA,
      k.score = k_score_i,
      verbose = FALSE
    ),
    error = function(e) {
      stop(
        "FindTransferAnchors failed for sample ", sample_id,
        ": ", conditionMessage(e), call. = FALSE
      )
    }
  )

  anchor_matrix_i <- slot(anchors_i, "anchors")

  if (
    nrow(anchor_matrix_i) == 0L ||
      !"cell2" %in% colnames(anchor_matrix_i)
  ) {
    stop("No usable query anchors were found for ", sample_id, ".")
  }

  available_query_anchor_cells_i <- length(
    unique(anchor_matrix_i[, "cell2"])
  )

  if (available_query_anchor_cells_i < 2L) {
    stop("Fewer than two query anchor cells were found for ", sample_id, ".")
  }

  maximum_k_weight <- as.integer(min(
    DEFAULT_K_WEIGHT,
    ncol(query_i) - 1L,
    available_query_anchor_cells_i - 1L
  ))

  candidate_k_weights <- unique(as.integer(c(
    maximum_k_weight, 30L, 20L, 10L, 5L, 2L, 1L
  )))

  candidate_k_weights <- candidate_k_weights[
    candidate_k_weights >= 1L &
      candidate_k_weights <= maximum_k_weight
  ]

  predictions_i <- NULL
  transfer_errors_i <- character(0)
  selected_k_weight_i <- NA_integer_

  for (candidate_k_weight in candidate_k_weights) {
    message(
      "Trying k.weight = ", candidate_k_weight,
      " with ", available_query_anchor_cells_i,
      " available query anchor cells."
    )

    transfer_attempt <- tryCatch(
      TransferData(
        anchorset = anchors_i,
        refdata = reference_state_labels,
        dims = seq_len(N_TRANSFER_DIMS),
        k.weight = candidate_k_weight,
        verbose = FALSE
      ),
      error = function(e) e
    )

    if (!inherits(transfer_attempt, "error")) {
      predictions_i <- transfer_attempt
      selected_k_weight_i <- candidate_k_weight
      break
    }

    transfer_errors_i <- c(
      transfer_errors_i,
      paste0(
        "k.weight=", candidate_k_weight, ": ",
        conditionMessage(transfer_attempt)
      )
    )
  }

  if (is.null(predictions_i)) {
    stop(
      "TransferData failed for sample ", sample_id,
      " after all valid k.weight attempts:\n",
      paste(transfer_errors_i, collapse = "\n"),
      call. = FALSE
    )
  }

  predictions_i$cell_id <- rownames(predictions_i)
  predictions_i$transfer_k_anchor <- k_anchor_i
  predictions_i$transfer_k_score <- k_score_i
  predictions_i$transfer_k_weight <- selected_k_weight_i
  predictions_i$available_query_anchor_cells <-
    available_query_anchor_cells_i

  sample_metadata_i <- query_metadata[
    predictions_i$cell_id,
    required_query_metadata,
    drop = FALSE
  ]
  rownames(sample_metadata_i) <- NULL

  predictions_i <- bind_cols(predictions_i, sample_metadata_i)

  if (!validate_prediction_chunk(
    predictions_i,
    expected_cells = cells_i,
    expected_states = expected_state_ids
  )) {
    stop("The completed prediction chunk failed validation for ", sample_id, ".")
  }

  saveRDS(predictions_i, temporary_output_file, compress = TRUE)

  saved_check <- tryCatch(
    readRDS(temporary_output_file),
    error = function(e) NULL
  )

  saved_valid <- validate_prediction_chunk(
    saved_check,
    expected_cells = cells_i,
    expected_states = expected_state_ids
  )

  rm(saved_check)

  if (!saved_valid) {
    stop("The temporary prediction chunk failed validation for ", sample_id, ".")
  }

  if (!file.rename(temporary_output_file, output_file)) {
    stop("Unable to finalize the validated chunk for ", sample_id, ".")
  }

  message(
    "Completed: ", sample_id,
    "; median score = ",
    round(median(predictions_i$prediction.score.max, na.rm = TRUE), 3)
  )

  rm(
    query_i,
    anchors_i,
    anchor_matrix_i,
    predictions_i,
    sample_metadata_i,
    transfer_errors_i
  )
  invisible(gc())
}

# ------------------------------------------------------------
# Assemble and validate all external predictions
# ------------------------------------------------------------

expected_chunk_files <- file.path(
  CHUNK_DIR,
  paste0(
    sprintf("%03d", seq_along(sample_ids)), "_",
    gsub("[^A-Za-z0-9_-]", "_", sample_ids),
    "_external_state_predictions.rds"
  )
)

missing_chunk_files <- expected_chunk_files[
  !file.exists(expected_chunk_files)
]

if (length(missing_chunk_files) > 0L) {
  stop(
    "External mapping ended with missing chunk(s):\n",
    paste(missing_chunk_files, collapse = "\n")
  )
}

message("Loading and validating all external prediction chunks...")

prediction_chunks <- lapply(
  seq_along(expected_chunk_files),
  function(chunk_index) {
    prediction_chunk <- readRDS(expected_chunk_files[[chunk_index]])
    sample_id <- sample_ids[[chunk_index]]
    cells_i <- query_cells[
      as.character(query_metadata$sample) == sample_id
    ]

    if (!validate_prediction_chunk(
      prediction_chunk,
      expected_cells = cells_i,
      expected_states = expected_state_ids
    )) {
      stop(
        "Final validation failed for chunk ",
        expected_chunk_files[[chunk_index]], "."
      )
    }

    prediction_chunk
  }
)

external_predictions <- bind_rows(prediction_chunks)
rm(prediction_chunks)
invisible(gc())

exact_cell_set_match <- setequal(
  external_predictions$cell_id,
  query_cells
)

metadata_match <- all(
  as.character(external_predictions$sample) ==
    as.character(query_metadata[external_predictions$cell_id, "sample"])
)

if (
  nrow(external_predictions) != length(query_cells) ||
    anyDuplicated(external_predictions$cell_id) ||
    anyNA(external_predictions$predicted.id) ||
    anyNA(external_predictions$prediction.score.max) ||
    !exact_cell_set_match ||
    !metadata_match
) {
  stop("The assembled external predictions failed final validation.")
}

external_predictions <- external_predictions %>%
  mutate(
    predicted.id = as.character(predicted.id)
  ) %>%
  left_join(
    state_lookup,
    by = c("predicted.id" = "immune_state_id")
  )

if (
  anyNA(external_predictions$immune_state_name) ||
    anyNA(external_predictions$state_family)
) {
  stop("State names or families are missing after lookup.")
}

external_predictions$predicted.id <- factor(
  external_predictions$predicted.id,
  levels = expected_state_ids
)

saveRDS(
  external_predictions,
  PREDICTION_FILE,
  compress = "xz"
)

# ------------------------------------------------------------
# External-transfer QC summaries
# ------------------------------------------------------------

sample_qc <- external_predictions %>%
  group_by(sample, donor, group, compartment) %>%
  summarise(
    predicted_cells = n(),
    predicted_states = n_distinct(predicted.id),
    mean_score = mean(prediction.score.max, na.rm = TRUE),
    median_score = median(prediction.score.max, na.rm = TRUE),
    q25_score = quantile(
      prediction.score.max,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    pct_below_050 = 100 * mean(prediction.score.max < 0.50),
    pct_below_075 = 100 * mean(prediction.score.max < 0.75),
    pct_below_090 = 100 * mean(prediction.score.max < 0.90),
    .groups = "drop"
  ) %>%
  arrange(compartment, group, sample)

state_qc <- external_predictions %>%
  group_by(predicted.id, immune_state_name, state_family) %>%
  summarise(
    predicted_cells = n(),
    samples_with_cells = n_distinct(sample),
    donors_with_cells = n_distinct(donor),
    mean_score = mean(prediction.score.max, na.rm = TRUE),
    median_score = median(prediction.score.max, na.rm = TRUE),
    q25_score = quantile(
      prediction.score.max,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    pct_below_050 = 100 * mean(prediction.score.max < 0.50),
    pct_below_075 = 100 * mean(prediction.score.max < 0.75),
    pct_below_090 = 100 * mean(prediction.score.max < 0.90),
    .groups = "drop"
  ) %>%
  arrange(predicted.id)

score_summary <- data.frame(
  predicted_cells = nrow(external_predictions),
  minimum = min(external_predictions$prediction.score.max),
  q25 = quantile(
    external_predictions$prediction.score.max,
    0.25,
    names = FALSE
  ),
  median = median(external_predictions$prediction.score.max),
  mean = mean(external_predictions$prediction.score.max),
  q75 = quantile(
    external_predictions$prediction.score.max,
    0.75,
    names = FALSE
  ),
  maximum = max(external_predictions$prediction.score.max),
  pct_below_050 = 100 * mean(
    external_predictions$prediction.score.max < 0.50
  ),
  pct_below_075 = 100 * mean(
    external_predictions$prediction.score.max < 0.75
  ),
  pct_below_090 = 100 * mean(
    external_predictions$prediction.score.max < 0.90
  )
)

final_validation <- data.frame(
  metric = c(
    "chunk_files",
    "expected_samples",
    "predicted_cells",
    "expected_query_cells",
    "unique_predicted_cells",
    "duplicated_cell_ids",
    "missing_predictions",
    "missing_scores",
    "predicted_states",
    "exact_cell_set_match",
    "metadata_match"
  ),
  value = c(
    length(expected_chunk_files),
    length(sample_ids),
    nrow(external_predictions),
    length(query_cells),
    n_distinct(external_predictions$cell_id),
    sum(duplicated(external_predictions$cell_id)),
    sum(is.na(external_predictions$predicted.id)),
    sum(is.na(external_predictions$prediction.score.max)),
    n_distinct(external_predictions$predicted.id),
    as.integer(exact_cell_set_match),
    as.integer(metadata_match)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_qc,
  file.path(TABLE_DIR, "GSE138266_state_transfer_sample_QC.csv"),
  row.names = FALSE
)

write.csv(
  state_qc,
  file.path(TABLE_DIR, "GSE138266_state_transfer_state_QC.csv"),
  row.names = FALSE
)

write.csv(
  score_summary,
  file.path(TABLE_DIR, "GSE138266_state_transfer_score_summary.csv"),
  row.names = FALSE
)

write.csv(
  final_validation,
  file.path(TABLE_DIR, "GSE138266_state_transfer_validation.csv"),
  row.names = FALSE
)

message("")
message("GSE138266 external state transfer completed.")
message("")
message("Final validation:")
print(final_validation, row.names = FALSE)
message("")
message("Global score summary:")
print(score_summary, row.names = FALSE)
message("")
message("State-level QC:")
print(state_qc, row.names = FALSE)
message("")
message("Saved predictions: ", PREDICTION_FILE)

