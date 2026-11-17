# ============================================================
# 21_transfer_curated_states_to_full_V3.R
# Transfer the 11 curated states to all remaining V3 T cells
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

REFERENCE_CHEMISTRY <- "V3"
N_TRANSFER_DIMS <- 30L
DEFAULT_K_ANCHOR <- 5L
DEFAULT_K_SCORE <- 30L
DEFAULT_K_WEIGHT <- 50L
MIN_QUERY_CELLS <- 3L

reference_file <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas.rds"
  )
)

full_tcells_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_tcells_raw.rds"
  )
)

chunk_dir <- file.path(
  PROCESSED_DIR,
  "curated_V3_state_transfer_chunks"
)

dir.create(
  chunk_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

required_files <- c(
  reference_file,
  full_tcells_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# Load reference and complete T-cell object
# ------------------------------------------------------------

message("Loading curated V3 reference...")

reference <- readRDS(
  reference_file
)

message("Loading complete T-cell object...")

full_tcells <- readRDS(
  full_tcells_file
)

if (
  !inherits(reference, "Seurat") ||
    !inherits(full_tcells, "Seurat")
) {
  stop(
    "Both input files must contain Seurat objects."
  )
}

if (
  !"RNA" %in% Assays(reference) ||
    !"RNA" %in% Assays(full_tcells)
) {
  stop(
    "The RNA assay is missing from one or both input objects."
  )
}

DefaultAssay(reference) <- "RNA"
DefaultAssay(full_tcells) <- "RNA"

required_reference_metadata <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family",
  "V_10X"
)

missing_reference_metadata <- setdiff(
  required_reference_metadata,
  colnames(reference@meta.data)
)

if (length(missing_reference_metadata) > 0) {
  stop(
    "Missing reference metadata: ",
    paste(
      missing_reference_metadata,
      collapse = ", "
    )
  )
}

required_query_metadata <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "V_10X",
  "exp_name"
)

missing_query_metadata <- setdiff(
  required_query_metadata,
  colnames(full_tcells@meta.data)
)

if (length(missing_query_metadata) > 0) {
  stop(
    "Missing full-object metadata: ",
    paste(
      missing_query_metadata,
      collapse = ", "
    )
  )
}

if (!"pca" %in% Reductions(reference)) {
  stop(
    "The curated reference does not contain a PCA reduction."
  )
}

reference_pca_dimensions <- ncol(
  Embeddings(
    reference,
    reduction = "pca"
  )
)

if (reference_pca_dimensions < N_TRANSFER_DIMS) {
  stop(
    "The reference PCA contains only ",
    reference_pca_dimensions,
    " dimensions; ",
    N_TRANSFER_DIMS,
    " are required."
  )
}

if (
  !all(
    colnames(reference) %in%
      colnames(full_tcells)
  )
) {
  stop(
    "Some curated reference cells are absent ",
    "from the complete T-cell object."
  )
}

reference_chemistry <- unique(
  as.character(
    reference$V_10X
  )
)

reference_chemistry <- reference_chemistry[
  !is.na(reference_chemistry)
]

if (
  length(reference_chemistry) != 1L ||
    reference_chemistry != REFERENCE_CHEMISTRY
) {
  stop(
    "The discovery reference is not chemistry-homogeneous V3."
  )
}

# ------------------------------------------------------------
# Reference labels
# ------------------------------------------------------------

if (
  anyNA(reference$immune_state_id) ||
    any(as.character(reference$immune_state_id) == "")
) {
  stop(
    "The curated reference contains missing immune-state IDs."
  )
}

expected_state_ids <- sort(
  unique(
    as.character(
      reference$immune_state_id
    )
  )
)

if (length(expected_state_ids) != 11L) {
  stop(
    "Expected 11 curated states, but found ",
    length(expected_state_ids),
    "."
  )
}

reference_state_labels <- setNames(
  as.character(
    reference$immune_state_id
  ),
  colnames(reference)
)

# ------------------------------------------------------------
# Mapping features
# ------------------------------------------------------------

features_use <- VariableFeatures(
  reference
)

if (length(features_use) < 500L) {
  features_use <- rownames(
    Loadings(
      reference[["pca"]]
    )
  )
}

features_use <- Reduce(
  intersect,
  list(
    features_use,
    rownames(reference),
    rownames(full_tcells)
  )
)

features_use <- features_use[
  !grepl(
    "^(MT-|RPL|RPS|HBA|HBB|HBD)",
    features_use
  )
]

features_use <- unique(
  features_use
)

if (length(features_use) < 500L) {
  stop(
    "Fewer than 500 suitable mapping features remain."
  )
}

message(
  "Mapping with ",
  length(features_use),
  " reference features."
)

# ------------------------------------------------------------
# Identify the remaining V3 query cells
# ------------------------------------------------------------

all_v3_cells <- colnames(
  full_tcells
)[
  !is.na(full_tcells$V_10X) &
    as.character(
      full_tcells$V_10X
    ) == REFERENCE_CHEMISTRY
]

query_cells <- setdiff(
  all_v3_cells,
  colnames(reference)
)

if (length(query_cells) == 0L) {
  stop(
    "No remaining V3 query cells were found."
  )
}

query_metadata <- full_tcells@meta.data[
  query_cells,
  required_query_metadata,
  drop = FALSE
]

if (
  anyNA(query_metadata$sample) ||
    any(as.character(query_metadata$sample) == "")
) {
  stop(
    "Missing sample identifiers among V3 query cells."
  )
}

sample_ids <- sort(
  unique(
    as.character(
      query_metadata$sample
    )
  )
)

if (length(sample_ids) == 0L) {
  stop(
    "No biological samples were found among V3 query cells."
  )
}

# Verify that each sample has internally consistent metadata.
sample_metadata_check <- query_metadata %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  group_by(
    sample
  ) %>%
  summarise(
    across(
      all_of(
        c(
          "donor",
          "group",
          "cohort",
          "V_10X",
          "exp_name"
        )
      ),
      ~ n_distinct(
        .x,
        na.rm = FALSE
      )
    ),
    .groups = "drop"
  )

inconsistent_samples <- sample_metadata_check %>%
  filter(
    if_any(
      -sample,
      ~ .x != 1L
    )
  )

if (nrow(inconsistent_samples) > 0L) {
  stop(
    "Inconsistent metadata were found within samples: ",
    paste(
      inconsistent_samples$sample,
      collapse = ", "
    ),
    "."
  )
}

# ------------------------------------------------------------
# Save transfer manifest and configuration
# ------------------------------------------------------------

transfer_manifest <- query_metadata %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  count(
    sample,
    donor,
    group,
    cohort,
    exp_name,
    V_10X,
    name = "query_cells"
  ) %>%
  arrange(
    sample
  )

if (nrow(transfer_manifest) != length(sample_ids)) {
  stop(
    "The transfer manifest does not contain exactly one row ",
    "per biological sample."
  )
}

write.csv(
  transfer_manifest,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_manifest.csv"
    )
  ),
  row.names = FALSE
)

transfer_configuration <- data.frame(
  reference_cells = ncol(reference),
  complete_tcells = ncol(full_tcells),
  complete_V3_cells = length(all_v3_cells),
  V3_query_cells = length(query_cells),
  V3_query_samples = length(sample_ids),
  mapping_features = length(features_use),
  transfer_dimensions = N_TRANSFER_DIMS,
  reference_states = length(expected_state_ids),
  default_k_anchor = DEFAULT_K_ANCHOR,
  default_k_score = DEFAULT_K_SCORE,
  default_k_weight = DEFAULT_K_WEIGHT,
  minimum_query_cells = MIN_QUERY_CELLS
)

write.csv(
  transfer_configuration,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_configuration.csv"
    )
  ),
  row.names = FALSE
)

message(
  length(query_cells),
  " V3 cells will be mapped across ",
  length(sample_ids),
  " biological samples."
)

# ------------------------------------------------------------
# Validate a saved prediction chunk
# ------------------------------------------------------------

validate_prediction_chunk <- function(
  prediction_object,
  expected_cells,
  expected_states
) {

  required_prediction_columns <- c(
    "cell_id",
    "predicted.id",
    "prediction.score.max"
  )

  if (
    is.null(prediction_object) ||
      !is.data.frame(prediction_object) ||
      !all(
        required_prediction_columns %in%
          colnames(prediction_object)
      )
  ) {
    return(FALSE)
  }

  if (
    nrow(prediction_object) != length(expected_cells) ||
      anyDuplicated(prediction_object$cell_id) ||
      !setequal(
        prediction_object$cell_id,
        expected_cells
      )
  ) {
    return(FALSE)
  }

  if (
    anyNA(prediction_object$predicted.id) ||
      anyNA(prediction_object$prediction.score.max) ||
      any(
        !is.finite(
          prediction_object$prediction.score.max
        )
      )
  ) {
    return(FALSE)
  }

  score_tolerance <- sqrt(
    .Machine$double.eps
  )

  if (
    any(
      prediction_object$prediction.score.max <
        -score_tolerance |
        prediction_object$prediction.score.max >
        1 + score_tolerance
    )
  ) {
    return(FALSE)
  }

  if (
    !all(
      as.character(
        prediction_object$predicted.id
      ) %in% expected_states
    )
  ) {
    return(FALSE)
  }

  TRUE
}

# ------------------------------------------------------------
# Recover label predictions when Seurat returns non-finite
# weights for an exact-distance tie in a very small query
# ------------------------------------------------------------

transfer_labels_with_tie_repair <- function(
  anchorset,
  reference_labels,
  expected_states,
  transfer_dims,
  k_weight
) {

  weights <- TransferData(
    anchorset = anchorset,
    refdata = reference_labels,
    dims = transfer_dims,
    k.weight = k_weight,
    only.weights = TRUE,
    verbose = FALSE
  )

  weights_dense <- as.matrix(
    weights
  )

  non_finite_columns <- which(
    apply(
      weights_dense,
      2,
      function(x) {
        any(
          !is.finite(x)
        )
      }
    )
  )

  if (length(non_finite_columns) == 0L) {
    stop(
      "The tie-repair fallback was invoked, but no ",
      "non-finite weight columns were found."
    )
  }

  repaired_query_cells <- character(
    length(non_finite_columns)
  )

  query_cell_names <- sub(
    "_query$",
    "",
    slot(
      anchorset,
      "query.cells"
    )
  )

  for (repair_index in seq_along(non_finite_columns)) {

    column_index <- non_finite_columns[
      repair_index
    ]

    tied_anchor_rows <- which(
      !is.finite(
        weights_dense[, column_index]
      )
    )

    if (length(tied_anchor_rows) == 0L) {
      stop(
        "Unable to identify tied anchors for query cell ",
        query_cell_names[column_index],
        "."
      )
    }

    weights_dense[, column_index] <- 0

    weights_dense[
      tied_anchor_rows,
      column_index
    ] <- 1 / length(tied_anchor_rows)

    repaired_query_cells[
      repair_index
    ] <- query_cell_names[
      column_index
    ]
  }

  weight_sums <- colSums(
    weights_dense
  )

  if (
    any(
      !is.finite(weight_sums)
    ) ||
      any(weight_sums <= 0)
  ) {
    stop(
      "Invalid weight sums remained after exact-tie repair."
    )
  }

  weights_dense <- sweep(
    weights_dense,
    2,
    weight_sums,
    "/"
  )

  anchor_table <- as.data.frame(
    slot(
      anchorset,
      "anchors"
    )
  )

  anchor_state_labels <- as.character(
    reference_labels[
      anchor_table$cell1
    ]
  )

  if (
    anyNA(anchor_state_labels) ||
      !all(
        anchor_state_labels %in%
          expected_states
      )
  ) {
    stop(
      "Unable to match one or more anchors to curated ",
      "reference-state labels."
    )
  }

  prediction_indicator <- matrix(
    0,
    nrow = nrow(anchor_table),
    ncol = length(expected_states),
    dimnames = list(
      NULL,
      expected_states
    )
  )

  state_column_indices <- match(
    anchor_state_labels,
    expected_states
  )

  prediction_indicator[
    cbind(
      seq_len(
        nrow(anchor_table)
      ),
      state_column_indices
    )
  ] <- 1

  prediction_scores <- t(
    weights_dense
  ) %*%
    prediction_indicator

  prediction_scores <- as.matrix(
    prediction_scores
  )

  rownames(prediction_scores) <- query_cell_names
  colnames(prediction_scores) <- expected_states

  if (
    any(
      !is.finite(
        prediction_scores
      )
    )
  ) {
    stop(
      "Non-finite prediction scores remained after ",
      "exact-tie repair."
    )
  }

  prediction_indices <- max.col(
    prediction_scores,
    ties.method = "first"
  )

  prediction_ids <- expected_states[
    prediction_indices
  ]

  prediction_maxima <- prediction_scores[
    cbind(
      seq_len(
        nrow(prediction_scores)
      ),
      prediction_indices
    )
  ]

  predictions <- data.frame(
    predicted.id = prediction_ids,
    prediction.score = prediction_scores,
    row.names = query_cell_names,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  predictions$prediction.score.max <-
    prediction_maxima

  list(
    predictions = predictions,
    repaired_query_cells = repaired_query_cells
  )
}

# ------------------------------------------------------------
# Transfer labels sample by sample
# ------------------------------------------------------------

for (i in seq_along(sample_ids)) {

  sample_id <- sample_ids[i]

  cells_i <- query_cells[
    as.character(
      query_metadata$sample
    ) == sample_id
  ]

  if (length(cells_i) == 0L) {
    message(
      "[",
      i,
      "/",
      length(sample_ids),
      "] ",
      sample_id,
      " has no query cells; skipping."
    )
    next
  }

  safe_sample_id <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    sample_id
  )

  output_file <- file.path(
    chunk_dir,
    paste0(
      sprintf(
        "%03d",
        i
      ),
      "_",
      safe_sample_id,
      "_curated_V3_predictions.rds"
    )
  )

  temporary_output_file <- paste0(
    output_file,
    ".tmp"
  )

  if (file.exists(output_file)) {

    existing_prediction <- tryCatch(
      readRDS(output_file),
      error = function(e) NULL
    )

    existing_valid <- validate_prediction_chunk(
      prediction_object = existing_prediction,
      expected_cells = cells_i,
      expected_states = expected_state_ids
    )

    rm(
      existing_prediction
    )

    if (existing_valid) {
      message(
        "[",
        i,
        "/",
        length(sample_ids),
        "] ",
        sample_id,
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
    unlink(
      temporary_output_file
    )
  }

  if (length(cells_i) < MIN_QUERY_CELLS) {
    stop(
      "Sample ",
      sample_id,
      " contains only ",
      length(cells_i),
      " query cells; at least ",
      MIN_QUERY_CELLS,
      " are required."
    )
  }

  message(
    "[",
    i,
    "/",
    length(sample_ids),
    "] Mapping ",
    sample_id,
    " (",
    length(cells_i),
    " cells)"
  )

  set.seed(
    123 + i
  )

  query_i <- subset(
    full_tcells,
    cells = cells_i
  )

  DefaultAssay(query_i) <- "RNA"

  if (
    !"data" %in%
      Layers(
        query_i[["RNA"]]
      )
  ) {
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

  # Seurat requires k.anchor and k.score to be smaller than
  # both the query and reference cell counts.
  k_anchor_i <- min(
    DEFAULT_K_ANCHOR,
    ncol(query_i) - 1L,
    ncol(reference) - 1L
  )

  k_score_i <- min(
    DEFAULT_K_SCORE,
    ncol(query_i) - 1L,
    ncol(reference) - 1L
  )

  k_anchor_i <- as.integer(
    k_anchor_i
  )

  k_score_i <- as.integer(
    k_score_i
  )

  if (
    k_anchor_i < 1L ||
      k_score_i < 1L
  ) {
    stop(
      "Unable to define valid anchor parameters for sample ",
      sample_id,
      "."
    )
  }

  message(
    "Using k.anchor = ",
    k_anchor_i,
    " and k.score = ",
    k_score_i,
    "."
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
      dims = seq_len(
        N_TRANSFER_DIMS
      ),
      k.anchor = k_anchor_i,
      k.filter = NA,
      k.score = k_score_i,
      verbose = FALSE
    ),
    error = function(e) {
      stop(
        "FindTransferAnchors failed for sample ",
        sample_id,
        " (",
        length(cells_i),
        " cells; k.anchor = ",
        k_anchor_i,
        "; k.score = ",
        k_score_i,
        "): ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  anchor_matrix_i <- slot(
    anchors_i,
    "anchors"
  )

  if (
    !is.matrix(anchor_matrix_i) &&
      !is.data.frame(anchor_matrix_i)
  ) {
    stop(
      "The anchor matrix has an unexpected format for sample ",
      sample_id,
      "."
    )
  }

  if (
    nrow(anchor_matrix_i) == 0L ||
      !"cell2" %in% colnames(anchor_matrix_i)
  ) {
    stop(
      "No usable query anchors were found for sample ",
      sample_id,
      "."
    )
  }

  available_query_anchor_cells_i <- length(
    unique(
      anchor_matrix_i[, "cell2"]
    )
  )

  if (available_query_anchor_cells_i < 2L) {
    stop(
      "Fewer than two query anchor cells were found for sample ",
      sample_id,
      "; label transfer cannot be performed reliably."
    )
  }

  # Seurat requires k.weight to be strictly smaller than the
  # number of available query anchor cells.
  k_weight_i <- min(
    DEFAULT_K_WEIGHT,
    ncol(query_i) - 1L,
    available_query_anchor_cells_i - 1L
  )

  k_weight_i <- as.integer(
    k_weight_i
  )

  if (k_weight_i < 1L) {
    stop(
      "Unable to define a valid k.weight for sample ",
      sample_id,
      "."
    )
  }

  message(
    "Using k.weight = ",
    k_weight_i,
    " with ",
    available_query_anchor_cells_i,
    " available query anchor cells."
  )

  transfer_error_i <- NULL
  tie_repair_result_i <- NULL
  used_exact_tie_repair_i <- FALSE

  predictions_i <- tryCatch(
    TransferData(
      anchorset = anchors_i,
      refdata = reference_state_labels,
      dims = seq_len(
        N_TRANSFER_DIMS
      ),
      k.weight = k_weight_i,
      verbose = FALSE
    ),
    error = function(e) {
      transfer_error_i <<- e
      NULL
    }
  )

  if (is.null(predictions_i)) {

    transfer_error_message_i <- conditionMessage(
      transfer_error_i
    )

    is_exact_tie_failure_i <- grepl(
      "invalid subscript type ['\"]list['\"]",
      transfer_error_message_i
    )

    if (!is_exact_tie_failure_i) {
      stop(
        "TransferData failed for sample ",
        sample_id,
        " (",
        length(cells_i),
        " cells; k.weight = ",
        k_weight_i,
        "; available query anchor cells = ",
        available_query_anchor_cells_i,
        "): ",
        transfer_error_message_i,
        call. = FALSE
      )
    }

    message(
      "Seurat encountered an exact-distance tie; ",
      "applying equal-anchor-weight repair."
    )

    tie_repair_result_i <- tryCatch(
      transfer_labels_with_tie_repair(
        anchorset = anchors_i,
        reference_labels = reference_state_labels,
        expected_states = expected_state_ids,
        transfer_dims = seq_len(
          N_TRANSFER_DIMS
        ),
        k_weight = k_weight_i
      ),
      error = function(e) {
        stop(
          "Exact-tie repair failed for sample ",
          sample_id,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    predictions_i <-
      tie_repair_result_i$predictions

    repaired_query_cells_i <-
      tie_repair_result_i$repaired_query_cells

    used_exact_tie_repair_i <- TRUE

  } else {

    repaired_query_cells_i <- character(0)
  }

  predictions_i$cell_id <- rownames(
    predictions_i
  )

  predictions_i$transfer_k_anchor <- k_anchor_i
  predictions_i$transfer_k_score <- k_score_i
  predictions_i$transfer_k_weight <- k_weight_i
  predictions_i$available_query_anchor_cells <-
    available_query_anchor_cells_i
  predictions_i$used_exact_tie_repair <-
    used_exact_tie_repair_i
  predictions_i$cell_received_tie_repair <-
    predictions_i$cell_id %in%
      repaired_query_cells_i

  sample_metadata_i <- query_metadata[
    predictions_i$cell_id,
    required_query_metadata,
    drop = FALSE
  ]

  rownames(sample_metadata_i) <- NULL

  predictions_i <- bind_cols(
    predictions_i,
    sample_metadata_i
  )

  if (
    !validate_prediction_chunk(
      prediction_object = predictions_i,
      expected_cells = cells_i,
      expected_states = expected_state_ids
    )
  ) {
    stop(
      "The completed prediction chunk failed validation for sample ",
      sample_id,
      "."
    )
  }

  saveRDS(
    predictions_i,
    temporary_output_file,
    compress = TRUE
  )

  saved_prediction_check <- tryCatch(
    readRDS(
      temporary_output_file
    ),
    error = function(e) NULL
  )

  saved_prediction_valid <- validate_prediction_chunk(
    prediction_object = saved_prediction_check,
    expected_cells = cells_i,
    expected_states = expected_state_ids
  )

  rm(
    saved_prediction_check
  )

  if (!saved_prediction_valid) {
    stop(
      "The temporary saved chunk failed validation for sample ",
      sample_id,
      "."
    )
  }

  if (
    !file.rename(
      temporary_output_file,
      output_file
    )
  ) {
    stop(
      "Unable to move the validated temporary chunk to ",
      output_file,
      "."
    )
  }

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
    anchor_matrix_i,
    predictions_i,
    sample_metadata_i,
    transfer_error_i,
    tie_repair_result_i,
    repaired_query_cells_i
  )

  gc(
    verbose = FALSE
  )
}

# ------------------------------------------------------------
# Final chunk audit
# ------------------------------------------------------------

expected_chunk_files <- file.path(
  chunk_dir,
  paste0(
    sprintf(
      "%03d",
      seq_along(sample_ids)
    ),
    "_",
    gsub(
      "[^A-Za-z0-9_-]",
      "_",
      sample_ids
    ),
    "_curated_V3_predictions.rds"
  )
)

missing_chunk_files <- expected_chunk_files[
  !file.exists(expected_chunk_files)
]

if (length(missing_chunk_files) > 0L) {
  stop(
    "Mapping finished with ",
    length(missing_chunk_files),
    " missing prediction chunks:\n",
    paste(
      missing_chunk_files,
      collapse = "\n"
    )
  )
}

message(
  "All ",
  length(expected_chunk_files),
  " V3 sample chunks are present."
)

message(
  "All available V3 samples have been processed."
)
