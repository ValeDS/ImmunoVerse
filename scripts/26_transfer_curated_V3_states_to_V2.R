# ============================================================
# 26_transfer_curated_V3_states_to_V2.R
# Cross-chemistry transfer of the 11 curated V3 immune states
# to all V2 T cells in GSE144744
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

set.seed(123)

REFERENCE_CHEMISTRY <- "V3"
QUERY_CHEMISTRY <- "V2"
N_TRANSFER_DIMS <- 30L
DEFAULT_K_ANCHOR <- 5L
DEFAULT_K_SCORE <- 30L
DEFAULT_K_WEIGHT <- 50L
MIN_QUERY_CELLS <- 3L

reference_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_Harmony_reference_curated_atlas.rds"
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
  "curated_V3_to_V2_state_transfer_chunks"
)

prediction_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_to_V2_curated_state_predictions.rds"
  )
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
# Validate and load inputs
# ------------------------------------------------------------

required_files <- c(
  reference_file,
  full_tcells_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Loading the curated V3 reference...")
reference <- readRDS(reference_file)

message("Loading the complete T-cell object...")
full_tcells <- readRDS(full_tcells_file)

if (
  !inherits(reference, "Seurat") ||
    !inherits(full_tcells, "Seurat")
) {
  stop("Both input files must contain Seurat objects.")
}

if (
  !"RNA" %in% Assays(reference) ||
    !"RNA" %in% Assays(full_tcells)
) {
  stop("The RNA assay is missing from one or both objects.")
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
  "cohort",
  "V_10X",
  "exp_name"
)

missing_query_metadata <- setdiff(
  required_query_metadata,
  colnames(full_tcells[[]])
)

if (length(missing_query_metadata) > 0L) {
  stop(
    "The complete T-cell object is missing metadata column(s): ",
    paste(missing_query_metadata, collapse = ", ")
  )
}

if (!"pca" %in% Reductions(reference)) {
  stop("The curated V3 reference does not contain a PCA reduction.")
}

reference_pca_dimensions <- ncol(
  Embeddings(reference, reduction = "pca")
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

reference_chemistry <- unique(
  as.character(reference$V_10X)
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
# Reference labels and mapping features
# ------------------------------------------------------------

if (
  anyNA(reference$immune_state_id) ||
    any(as.character(reference$immune_state_id) == "")
) {
  stop("The curated reference contains missing immune-state IDs.")
}

expected_state_ids <- sort(
  unique(as.character(reference$immune_state_id))
)

if (length(expected_state_ids) != 11L) {
  stop(
    "Expected 11 curated immune states but found ",
    length(expected_state_ids),
    "."
  )
}

reference_state_labels <- setNames(
  as.character(reference$immune_state_id),
  colnames(reference)
)

features_use <- VariableFeatures(reference)

if (length(features_use) < 500L) {
  features_use <- rownames(
    Loadings(reference[["pca"]])
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

features_use <- unique(features_use)

if (length(features_use) < 500L) {
  stop(
    "Fewer than 500 suitable cross-chemistry mapping ",
    "features remain."
  )
}

message(
  "Cross-chemistry mapping will use ",
  length(features_use),
  " reference features."
)

# ------------------------------------------------------------
# Identify and audit the V2 query cells
# ------------------------------------------------------------

query_cells <- colnames(full_tcells)[
  !is.na(full_tcells$V_10X) &
    as.character(full_tcells$V_10X) ==
      QUERY_CHEMISTRY
]

if (length(query_cells) == 0L) {
  stop("No V2 query cells were found.")
}

if (any(query_cells %in% colnames(reference))) {
  stop(
    "V2 query cells unexpectedly overlap the V3 reference."
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
  stop("Missing sample identifiers were found among V2 cells.")
}

sample_ids <- sort(
  unique(as.character(query_metadata$sample))
)

sample_metadata_check <- query_metadata %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(sample) %>%
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
      ~ n_distinct(.x, na.rm = FALSE)
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
    "Inconsistent metadata were found within V2 sample(s): ",
    paste(inconsistent_samples$sample, collapse = ", ")
  )
}

transfer_manifest <- query_metadata %>%
  tibble::rownames_to_column("cell_id") %>%
  count(
    sample,
    donor,
    group,
    cohort,
    exp_name,
    V_10X,
    name = "query_cells"
  ) %>%
  arrange(sample)

if (
  nrow(transfer_manifest) != 22L ||
    length(sample_ids) != 22L
) {
  stop(
    "Expected 22 V2 samples but found ",
    length(sample_ids),
    "."
  )
}

group_audit <- transfer_manifest %>%
  count(cohort, group, name = "n_samples")

expected_group_audit <- data.frame(
  cohort = c("RRMS_HI", "RRMS_HI"),
  group = c("HI2", "MS2"),
  n_samples = c(11L, 11L),
  stringsAsFactors = FALSE
)

group_audit_comparison <- group_audit %>%
  arrange(cohort, group)

expected_group_audit <- expected_group_audit %>%
  arrange(cohort, group)

if (!identical(
  as.data.frame(group_audit_comparison),
  as.data.frame(expected_group_audit)
)) {
  stop(
    "The V2 sample composition is not the expected ",
    "11 HI2 plus 11 MS2 design."
  )
}

write.csv(
  transfer_manifest,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_transfer_manifest.csv"
    )
  ),
  row.names = FALSE
)

transfer_configuration <- data.frame(
  reference_chemistry = REFERENCE_CHEMISTRY,
  query_chemistry = QUERY_CHEMISTRY,
  reference_cells = ncol(reference),
  complete_tcells = ncol(full_tcells),
  V2_query_cells = length(query_cells),
  V2_query_samples = length(sample_ids),
  mapping_features = length(features_use),
  transfer_dimensions = N_TRANSFER_DIMS,
  reference_states = length(expected_state_ids),
  default_k_anchor = DEFAULT_K_ANCHOR,
  default_k_score = DEFAULT_K_SCORE,
  default_k_weight = DEFAULT_K_WEIGHT
)

write.csv(
  transfer_configuration,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_transfer_configuration.csv"
    )
  ),
  row.names = FALSE
)

message(
  length(query_cells),
  " V2 cells will be mapped across ",
  length(sample_ids),
  " biological samples."
)

# ------------------------------------------------------------
# Prediction-chunk validation
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
    "cohort",
    "V_10X"
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
      any(!is.finite(
        prediction_object$prediction.score.max
      ))
  ) {
    return(FALSE)
  }

  score_tolerance <- sqrt(.Machine$double.eps)

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
      as.character(prediction_object$predicted.id) %in%
        expected_states
    ) ||
      any(as.character(prediction_object$V_10X) !=
        QUERY_CHEMISTRY)
  ) {
    return(FALSE)
  }

  TRUE
}

# ------------------------------------------------------------
# Transfer labels one V2 sample at a time
# ------------------------------------------------------------

for (i in seq_along(sample_ids)) {

  sample_id <- sample_ids[[i]]

  cells_i <- query_cells[
    as.character(query_metadata$sample) == sample_id
  ]

  if (length(cells_i) < MIN_QUERY_CELLS) {
    stop(
      "V2 sample ",
      sample_id,
      " contains only ",
      length(cells_i),
      " cells; at least ",
      MIN_QUERY_CELLS,
      " are required."
    )
  }

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
      "_V3_to_V2_predictions.rds"
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
      existing_prediction,
      expected_cells = cells_i,
      expected_states = expected_state_ids
    )

    rm(existing_prediction)

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
      "An invalid or incomplete V2 chunk already exists: ",
      output_file,
      ". Move or inspect it before restarting."
    )
  }

  if (file.exists(temporary_output_file)) {
    unlink(temporary_output_file)
  }

  message(
    "[",
    i,
    "/",
    length(sample_ids),
    "] Mapping V2 sample ",
    sample_id,
    " (",
    length(cells_i),
    " cells)"
  )

  set.seed(1000 + i)

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

  k_anchor_i <- as.integer(
    min(
      DEFAULT_K_ANCHOR,
      ncol(query_i) - 1L,
      ncol(reference) - 1L
    )
  )

  k_score_i <- as.integer(
    min(
      DEFAULT_K_SCORE,
      ncol(query_i) - 1L,
      ncol(reference) - 1L
    )
  )

  if (
    k_anchor_i < 1L ||
      k_score_i < 1L
  ) {
    stop(
      "Unable to define valid anchor parameters for ",
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
      dims = seq_len(N_TRANSFER_DIMS),
      k.anchor = k_anchor_i,
      k.filter = NA,
      k.score = k_score_i,
      verbose = FALSE
    ),
    error = function(e) {
      stop(
        "FindTransferAnchors failed for V2 sample ",
        sample_id,
        ": ",
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
    nrow(anchor_matrix_i) == 0L ||
      !"cell2" %in% colnames(anchor_matrix_i)
  ) {
    stop(
      "No usable V2 query anchors were found for ",
      sample_id,
      "."
    )
  }

  available_query_anchor_cells_i <- length(
    unique(anchor_matrix_i[, "cell2"])
  )

  if (available_query_anchor_cells_i < 2L) {
    stop(
      "Fewer than two query anchor cells were found for ",
      sample_id,
      "."
    )
  }

  maximum_k_weight <- as.integer(
    min(
      DEFAULT_K_WEIGHT,
      ncol(query_i) - 1L,
      available_query_anchor_cells_i - 1L
    )
  )

  candidate_k_weights <- unique(
    as.integer(
      c(
        maximum_k_weight,
        30L,
        20L,
        10L,
        5L,
        2L,
        1L
      )
    )
  )

  candidate_k_weights <- candidate_k_weights[
    candidate_k_weights >= 1L &
      candidate_k_weights <= maximum_k_weight
  ]

  predictions_i <- NULL
  transfer_errors_i <- character(0)
  selected_k_weight_i <- NA_integer_

  for (candidate_k_weight in candidate_k_weights) {

    message(
      "Trying k.weight = ",
      candidate_k_weight,
      " with ",
      available_query_anchor_cells_i,
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
        "k.weight=",
        candidate_k_weight,
        ": ",
        conditionMessage(transfer_attempt)
      )
    )
  }

  if (is.null(predictions_i)) {
    stop(
      "TransferData failed for V2 sample ",
      sample_id,
      " after all valid k.weight attempts:\n",
      paste(transfer_errors_i, collapse = "\n"),
      call. = FALSE
    )
  }

  predictions_i$cell_id <- rownames(predictions_i)
  predictions_i$transfer_k_anchor <- k_anchor_i
  predictions_i$transfer_k_score <- k_score_i
  predictions_i$transfer_k_weight <-
    selected_k_weight_i
  predictions_i$available_query_anchor_cells <-
    available_query_anchor_cells_i

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

  if (!validate_prediction_chunk(
    predictions_i,
    expected_cells = cells_i,
    expected_states = expected_state_ids
  )) {
    stop(
      "The completed V2 prediction chunk failed validation ",
      "for ",
      sample_id,
      "."
    )
  }

  saveRDS(
    predictions_i,
    temporary_output_file,
    compress = TRUE
  )

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
    stop(
      "The temporary V2 chunk failed validation for ",
      sample_id,
      "."
    )
  }

  if (!file.rename(
    temporary_output_file,
    output_file
  )) {
    stop(
      "Unable to finalize the validated chunk for ",
      sample_id,
      "."
    )
  }

  message(
    "Completed: ",
    sample_id,
    "; median cross-chemistry score = ",
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
    transfer_errors_i
  )

  invisible(gc())
}

# ------------------------------------------------------------
# Assemble and validate all V2 predictions
# ------------------------------------------------------------

expected_chunk_files <- file.path(
  chunk_dir,
  paste0(
    sprintf("%03d", seq_along(sample_ids)),
    "_",
    gsub(
      "[^A-Za-z0-9_-]",
      "_",
      sample_ids
    ),
    "_V3_to_V2_predictions.rds"
  )
)

missing_chunk_files <- expected_chunk_files[
  !file.exists(expected_chunk_files)
]

if (length(missing_chunk_files) > 0L) {
  stop(
    "Cross-chemistry mapping ended with missing chunk(s):\n",
    paste(missing_chunk_files, collapse = "\n")
  )
}

message("Loading and validating all 22 V2 chunks...")

prediction_chunks <- lapply(
  seq_along(expected_chunk_files),
  function(chunk_index) {
    prediction_chunk <- readRDS(
      expected_chunk_files[[chunk_index]]
    )

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
        expected_chunk_files[[chunk_index]],
        "."
      )
    }

    prediction_chunk
  }
)

v2_predictions <- bind_rows(prediction_chunks)

rm(prediction_chunks)
invisible(gc())

exact_cell_set_match <- setequal(
  v2_predictions$cell_id,
  query_cells
)

metadata_match <- all(
  as.character(
    v2_predictions$sample
  ) ==
    as.character(
      query_metadata[
        v2_predictions$cell_id,
        "sample"
      ]
    )
) &&
  all(
    as.character(
      v2_predictions$group
    ) ==
      as.character(
        query_metadata[
          v2_predictions$cell_id,
          "group"
        ]
      )
  )

final_validation <- data.frame(
  metric = c(
    "chunk_files",
    "expected_V2_samples",
    "predicted_cells",
    "expected_V2_cells",
    "unique_predicted_cells",
    "duplicated_cell_ids",
    "missing_predictions",
    "missing_scores",
    "exact_cell_set_match",
    "metadata_match",
    "predicted_states"
  ),
  value = c(
    length(expected_chunk_files),
    length(sample_ids),
    nrow(v2_predictions),
    length(query_cells),
    n_distinct(v2_predictions$cell_id),
    sum(duplicated(v2_predictions$cell_id)),
    sum(is.na(v2_predictions$predicted.id)),
    sum(is.na(v2_predictions$prediction.score.max)),
    as.integer(exact_cell_set_match),
    as.integer(metadata_match),
    n_distinct(v2_predictions$predicted.id)
  ),
  stringsAsFactors = FALSE
)

if (
  nrow(v2_predictions) != length(query_cells) ||
    anyDuplicated(v2_predictions$cell_id) ||
    anyNA(v2_predictions$predicted.id) ||
    anyNA(v2_predictions$prediction.score.max) ||
    !exact_cell_set_match ||
    !metadata_match
) {
  stop(
    "The assembled V2 cross-chemistry predictions failed ",
    "final validation."
  )
}

v2_predictions$predicted.id <- factor(
  as.character(v2_predictions$predicted.id),
  levels = expected_state_ids
)

saveRDS(
  v2_predictions,
  prediction_file,
  compress = "xz"
)

# ------------------------------------------------------------
# Cross-chemistry QC summaries
# ------------------------------------------------------------

sample_qc <- v2_predictions %>%
  group_by(
    sample,
    donor,
    group,
    cohort,
    V_10X
  ) %>%
  summarise(
    predicted_cells = n(),
    predicted_states = n_distinct(predicted.id),
    mean_score = mean(
      prediction.score.max,
      na.rm = TRUE
    ),
    median_score = median(
      prediction.score.max,
      na.rm = TRUE
    ),
    q25_score = quantile(
      prediction.score.max,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    pct_below_050 = 100 * mean(
      prediction.score.max < 0.50,
      na.rm = TRUE
    ),
    pct_below_075 = 100 * mean(
      prediction.score.max < 0.75,
      na.rm = TRUE
    ),
    pct_below_090 = 100 * mean(
      prediction.score.max < 0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(sample)

state_qc <- v2_predictions %>%
  group_by(predicted.id) %>%
  summarise(
    predicted_cells = n(),
    samples_with_cells = n_distinct(sample),
    donors_with_cells = n_distinct(donor),
    mean_score = mean(
      prediction.score.max,
      na.rm = TRUE
    ),
    median_score = median(
      prediction.score.max,
      na.rm = TRUE
    ),
    q25_score = quantile(
      prediction.score.max,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    pct_below_050 = 100 * mean(
      prediction.score.max < 0.50,
      na.rm = TRUE
    ),
    pct_below_075 = 100 * mean(
      prediction.score.max < 0.75,
      na.rm = TRUE
    ),
    pct_below_090 = 100 * mean(
      prediction.score.max < 0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    predicted.id = factor(
      expected_state_ids,
      levels = expected_state_ids
    ),
    fill = list(
      predicted_cells = 0L,
      samples_with_cells = 0L,
      donors_with_cells = 0L
    )
  ) %>%
  arrange(predicted.id)

score_summary <- data.frame(
  minimum = min(
    v2_predictions$prediction.score.max,
    na.rm = TRUE
  ),
  q25 = quantile(
    v2_predictions$prediction.score.max,
    0.25,
    na.rm = TRUE,
    names = FALSE
  ),
  median = median(
    v2_predictions$prediction.score.max,
    na.rm = TRUE
  ),
  mean = mean(
    v2_predictions$prediction.score.max,
    na.rm = TRUE
  ),
  q75 = quantile(
    v2_predictions$prediction.score.max,
    0.75,
    na.rm = TRUE,
    names = FALSE
  ),
  maximum = max(
    v2_predictions$prediction.score.max,
    na.rm = TRUE
  ),
  pct_below_050 = 100 * mean(
    v2_predictions$prediction.score.max < 0.50,
    na.rm = TRUE
  ),
  pct_below_075 = 100 * mean(
    v2_predictions$prediction.score.max < 0.75,
    na.rm = TRUE
  ),
  pct_below_090 = 100 * mean(
    v2_predictions$prediction.score.max < 0.90,
    na.rm = TRUE
  )
)

write.csv(
  final_validation,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_transfer_validation.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  sample_qc,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_sample_QC.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  state_qc,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_state_QC.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  score_summary,
  file.path(
    TAB_DIR,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_to_V2_score_summary.csv"
    )
  ),
  row.names = FALSE
)

message("\nFinal V3-to-V2 transfer validation:")
print(as.data.frame(final_validation), row.names = FALSE)

message("\nCross-chemistry score summary:")
print(as.data.frame(score_summary), row.names = FALSE)

message("\nV2 sample-level QC:")
print(as.data.frame(sample_qc), row.names = FALSE)

message("\nV2 state-level QC:")
print(as.data.frame(state_qc), row.names = FALSE)

message("\nSaved V2 predictions: ", prediction_file)
message(
  "Cross-chemistry V3-to-V2 immune-state transfer completed."
)
