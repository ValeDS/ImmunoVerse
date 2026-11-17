# ============================================================
# 22_assemble_and_validate_full_V3_atlas.R
# Assemble the curated reference and transferred V3 cells,
# preserve curation exclusions, and generate final QC outputs
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

REFERENCE_CHEMISTRY <- "V3"
EXPECTED_STATES <- 11L
EXPECTED_QUERY_SAMPLES <- 53L

curated_reference_file <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas.rds"
  )
)

original_reference_file <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    ".rds"
  )
)

full_tcells_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_tcells_raw.rds"
  )
)

atlas_definition_file <- file.path(
  TAB_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_immune_state_atlas.csv"
  )
)

excluded_definition_file <- file.path(
  TAB_DIR,
  paste0(
    OUT_PREFIX,
    "_excluded_state_definitions.csv"
  )
)

chunk_dir <- file.path(
  PROCESSED_DIR,
  "curated_V3_state_transfer_chunks"
)

qc_dir <- file.path(
  TAB_DIR,
  "curated_V3_mapping_QC"
)

output_object_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_curated_atlas_full_mapped.rds"
  )
)

output_assignments_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_curated_atlas_full_assignments.rds"
  )
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Validate required inputs
# ------------------------------------------------------------

required_files <- c(
  curated_reference_file,
  original_reference_file,
  full_tcells_file,
  atlas_definition_file,
  excluded_definition_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

chunk_files <- sort(
  list.files(
    chunk_dir,
    pattern =
      "_curated_V3_predictions\\.rds$",
    full.names = TRUE
  )
)

if (length(chunk_files) != EXPECTED_QUERY_SAMPLES) {
  stop(
    "Expected ",
    EXPECTED_QUERY_SAMPLES,
    " prediction chunks, but found ",
    length(chunk_files),
    "."
  )
}

# ------------------------------------------------------------
# Load state definitions
# ------------------------------------------------------------

atlas_definition <- read.csv(
  atlas_definition_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_atlas_columns <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family",
  "confidence",
  "interpretation"
)

if (
  !all(
    required_atlas_columns %in%
      colnames(atlas_definition)
  )
) {
  stop(
    "The curated atlas definition is missing required columns."
  )
}

atlas_lookup <- atlas_definition %>%
  select(
    all_of(
      required_atlas_columns
    )
  ) %>%
  distinct()

if (
  nrow(atlas_lookup) != EXPECTED_STATES ||
    anyDuplicated(atlas_lookup$immune_state_id)
) {
  stop(
    "The atlas lookup must contain exactly ",
    EXPECTED_STATES,
    " unique immune-state IDs."
  )
}

expected_state_ids <- sort(
  atlas_lookup$immune_state_id
)

excluded_definition <- read.csv(
  excluded_definition_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_excluded_columns <- c(
  "V3_state_cluster",
  "exclusion_reason"
)

if (
  !all(
    required_excluded_columns %in%
      colnames(excluded_definition)
  )
) {
  stop(
    "The excluded-state definition is missing required columns."
  )
}

excluded_definition <- excluded_definition %>%
  mutate(
    V3_state_cluster =
      as.character(V3_state_cluster)
  )

# ------------------------------------------------------------
# Load curated reference metadata
# ------------------------------------------------------------

message("Loading curated V3 reference metadata...")

curated_reference <- readRDS(
  curated_reference_file
)

if (!inherits(curated_reference, "Seurat")) {
  stop(
    "The curated reference file does not contain a Seurat object."
  )
}

required_reference_metadata <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family"
)

if (
  !all(
    required_reference_metadata %in%
      colnames(curated_reference@meta.data)
  )
) {
  stop(
    "The curated reference is missing immune-state metadata."
  )
}

reference_cell_ids <- colnames(
  curated_reference
)

reference_assignments <- curated_reference@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  transmute(
    cell_id,
    curated_immune_state_id =
      as.character(immune_state_id),
    transfer_candidate_state_id =
      as.character(immune_state_id),
    transfer_prediction_score =
      NA_real_,
    state_assignment_origin =
      "curated_reference",
    atlas_inclusion_status =
      "included",
    curation_exclusion_reason =
      NA_character_,
    transfer_k_anchor =
      NA_integer_,
    transfer_k_score =
      NA_integer_,
    transfer_k_weight =
      NA_integer_,
    available_query_anchor_cells =
      NA_integer_,
    used_exact_tie_repair =
      FALSE,
    cell_received_tie_repair =
      FALSE
  )

if (
  nrow(reference_assignments) !=
    length(reference_cell_ids) ||
    anyDuplicated(reference_assignments$cell_id)
) {
  stop(
    "Invalid cell identifiers in the curated reference."
  )
}

if (
  !all(
    reference_assignments$curated_immune_state_id %in%
      expected_state_ids
  )
) {
  stop(
    "Unexpected immune-state IDs in the curated reference."
  )
}

rm(
  curated_reference
)

gc(
  verbose = FALSE
)

# ------------------------------------------------------------
# Recover cells intentionally excluded during V3 curation
# ------------------------------------------------------------

message("Recovering the original V3 curation exclusions...")

original_reference <- readRDS(
  original_reference_file
)

if (
  !inherits(original_reference, "Seurat") ||
    !"V3_state_cluster" %in%
      colnames(original_reference@meta.data)
) {
  stop(
    "The original V3 reference lacks V3_state_cluster metadata."
  )
}

original_reference_metadata <-
  original_reference@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  transmute(
    cell_id,
    V3_state_cluster =
      as.character(V3_state_cluster)
  )

excluded_cells <- original_reference_metadata %>%
  filter(
    V3_state_cluster %in%
      excluded_definition$V3_state_cluster
  ) %>%
  left_join(
    excluded_definition,
    by = "V3_state_cluster"
  )

unexpected_removed_cells <- setdiff(
  setdiff(
    colnames(original_reference),
    reference_cell_ids
  ),
  excluded_cells$cell_id
)

if (length(unexpected_removed_cells) > 0L) {
  stop(
    length(unexpected_removed_cells),
    " cells were removed from the curated reference without ",
    "a registered curation-exclusion reason."
  )
}

if (
  anyDuplicated(excluded_cells$cell_id) ||
    anyNA(excluded_cells$exclusion_reason)
) {
  stop(
    "Invalid excluded-cell audit information."
  )
}

rm(
  original_reference,
  original_reference_metadata
)

gc(
  verbose = FALSE
)

# ------------------------------------------------------------
# Load and validate all transferred predictions
# ------------------------------------------------------------

message(
  "Loading ",
  length(chunk_files),
  " transfer chunks..."
)

prediction_chunks <- lapply(
  chunk_files,
  readRDS
)

required_prediction_columns <- c(
  "cell_id",
  "predicted.id",
  "prediction.score.max",
  "sample",
  "donor",
  "group",
  "cohort",
  "V_10X",
  "exp_name"
)

invalid_chunk_columns <- vapply(
  prediction_chunks,
  function(x) {
    !is.data.frame(x) ||
      !all(
        required_prediction_columns %in%
          colnames(x)
      )
  },
  logical(1)
)

if (any(invalid_chunk_columns)) {
  stop(
    "One or more prediction chunks have an invalid structure: ",
    paste(
      basename(
        chunk_files[
          invalid_chunk_columns
        ]
      ),
      collapse = ", "
    )
  )
}

standardize_optional_column <- function(
  x,
  column_name,
  default_value
) {

  if (!column_name %in% colnames(x)) {
    x[[column_name]] <- default_value
  }

  x
}

prediction_chunks <- lapply(
  prediction_chunks,
  function(x) {

    x <- standardize_optional_column(
      x,
      "transfer_k_anchor",
      NA_integer_
    )

    x <- standardize_optional_column(
      x,
      "transfer_k_score",
      NA_integer_
    )

    x <- standardize_optional_column(
      x,
      "transfer_k_weight",
      NA_integer_
    )

    x <- standardize_optional_column(
      x,
      "available_query_anchor_cells",
      NA_integer_
    )

    x <- standardize_optional_column(
      x,
      "used_exact_tie_repair",
      FALSE
    )

    x <- standardize_optional_column(
      x,
      "cell_received_tie_repair",
      FALSE
    )

    x
  }
)

transferred_predictions <- bind_rows(
  prediction_chunks
) %>%
  transmute(
    cell_id =
      as.character(cell_id),
    transfer_candidate_state_id =
      as.character(predicted.id),
    transfer_prediction_score =
      as.numeric(prediction.score.max),
    sample =
      as.character(sample),
    donor =
      as.character(donor),
    group =
      as.character(group),
    cohort =
      as.character(cohort),
    V_10X =
      as.character(V_10X),
    exp_name =
      as.character(exp_name),
    transfer_k_anchor =
      as.integer(transfer_k_anchor),
    transfer_k_score =
      as.integer(transfer_k_score),
    transfer_k_weight =
      as.integer(transfer_k_weight),
    available_query_anchor_cells =
      as.integer(available_query_anchor_cells),
    used_exact_tie_repair =
      as.logical(used_exact_tie_repair),
    cell_received_tie_repair =
      as.logical(cell_received_tie_repair)
  )

rm(
  prediction_chunks
)

gc(
  verbose = FALSE
)

if (
  anyDuplicated(transferred_predictions$cell_id) ||
    anyNA(transferred_predictions$cell_id)
) {
  stop(
    "Duplicated or missing cell IDs were found among predictions."
  )
}

if (
  anyNA(transferred_predictions$transfer_candidate_state_id) ||
    anyNA(transferred_predictions$transfer_prediction_score)
) {
  stop(
    "Missing transferred labels or scores were found."
  )
}

if (
  !all(
    transferred_predictions$transfer_candidate_state_id %in%
      expected_state_ids
  )
) {
  stop(
    "Unexpected immune-state IDs were found among predictions."
  )
}

score_tolerance <- sqrt(
  .Machine$double.eps
)

if (
  any(
    !is.finite(
      transferred_predictions$transfer_prediction_score
    )
  ) ||
    any(
      transferred_predictions$transfer_prediction_score <
        -score_tolerance |
        transferred_predictions$transfer_prediction_score >
        1 + score_tolerance
    )
) {
  stop(
    "Invalid transfer scores were found."
  )
}

# Clip harmless floating-point excursions above one.
transferred_predictions <-
  transferred_predictions %>%
  mutate(
    transfer_prediction_score = pmin(
      1,
      pmax(
        0,
        transfer_prediction_score
      )
    )
  )

# ------------------------------------------------------------
# Load the complete T-cell object and define the full V3 set
# ------------------------------------------------------------

message("Loading the complete T-cell object...")

full_tcells <- readRDS(
  full_tcells_file
)

if (
  !inherits(full_tcells, "Seurat") ||
    !"V_10X" %in%
      colnames(full_tcells@meta.data)
) {
  stop(
    "The complete T-cell object is invalid or lacks V_10X."
  )
}

all_v3_cells <- colnames(
  full_tcells
)[
  !is.na(full_tcells$V_10X) &
    as.character(
      full_tcells$V_10X
    ) == REFERENCE_CHEMISTRY
]

expected_query_cells <- setdiff(
  all_v3_cells,
  reference_cell_ids
)

if (
  !setequal(
    transferred_predictions$cell_id,
    expected_query_cells
  )
) {
  stop(
    "The transferred prediction cell set does not exactly match ",
    "the expected non-reference V3 cell set."
  )
}

if (
  !all(
    excluded_cells$cell_id %in%
      transferred_predictions$cell_id
  )
) {
  stop(
    "One or more curation-excluded cells are absent from the ",
    "transfer predictions."
  )
}

# Verify that chunk metadata match the complete object.
full_query_metadata <- full_tcells@meta.data[
  transferred_predictions$cell_id,
  c(
    "sample",
    "donor",
    "group",
    "cohort",
    "V_10X",
    "exp_name"
  ),
  drop = FALSE
]

metadata_match <- all(
  mapply(
    function(chunk_value, object_value) {
      identical(
        as.character(chunk_value),
        as.character(object_value)
      )
    },
    transferred_predictions[
      c(
        "sample",
        "donor",
        "group",
        "cohort",
        "V_10X",
        "exp_name"
      )
    ],
    full_query_metadata[
      c(
        "sample",
        "donor",
        "group",
        "cohort",
        "V_10X",
        "exp_name"
      )
    ]
  )
)

if (!metadata_match) {
  stop(
    "Prediction-chunk metadata do not match the complete object."
  )
}

# ------------------------------------------------------------
# Construct accepted and excluded assignments
# ------------------------------------------------------------

excluded_lookup <- excluded_cells %>%
  select(
    cell_id,
    exclusion_reason
  )

query_assignments <- transferred_predictions %>%
  left_join(
    excluded_lookup,
    by = "cell_id"
  ) %>%
  mutate(
    is_curation_excluded =
      !is.na(exclusion_reason),
    curated_immune_state_id = if_else(
      is_curation_excluded,
      NA_character_,
      transfer_candidate_state_id
    ),
    state_assignment_origin = if_else(
      is_curation_excluded,
      "curation_excluded",
      "label_transfer"
    ),
    atlas_inclusion_status = if_else(
      is_curation_excluded,
      "excluded",
      "included"
    ),
    curation_exclusion_reason =
      exclusion_reason
  ) %>%
  select(
    cell_id,
    curated_immune_state_id,
    transfer_candidate_state_id,
    transfer_prediction_score,
    state_assignment_origin,
    atlas_inclusion_status,
    curation_exclusion_reason,
    transfer_k_anchor,
    transfer_k_score,
    transfer_k_weight,
    available_query_anchor_cells,
    used_exact_tie_repair,
    cell_received_tie_repair
  )

all_assignments <- bind_rows(
  reference_assignments,
  query_assignments
)

if (
  nrow(all_assignments) !=
    length(all_v3_cells) ||
    anyDuplicated(all_assignments$cell_id) ||
    !setequal(
      all_assignments$cell_id,
      all_v3_cells
    )
) {
  stop(
    "The assembled assignments do not exactly cover the full V3 set."
  )
}

all_assignments <- all_assignments %>%
  left_join(
    atlas_lookup %>%
      rename(
        curated_immune_state_id =
          immune_state_id,
        curated_immune_state_name =
          immune_state_name,
        curated_state_family =
          state_family,
        curated_state_definition_confidence =
          confidence,
        curated_state_interpretation =
          interpretation
      ),
    by = "curated_immune_state_id"
  ) %>%
  mutate(
    transfer_score_band = case_when(
      state_assignment_origin ==
        "curated_reference" ~
        "Not applicable: reference",
      atlas_inclusion_status ==
        "excluded" ~
        "Not applicable: excluded",
      transfer_prediction_score >= 0.90 ~
        "Score >= 0.90",
      transfer_prediction_score >= 0.75 ~
        "0.75 <= score < 0.90",
      transfer_prediction_score >= 0.50 ~
        "0.50 <= score < 0.75",
      TRUE ~
        "Score < 0.50"
    )
  )

included_assignments <- all_assignments %>%
  filter(
    atlas_inclusion_status ==
      "included"
  )

if (
  anyNA(
    included_assignments$curated_immune_state_id
  ) ||
    anyNA(
      included_assignments$curated_immune_state_name
    ) ||
    anyNA(
      included_assignments$curated_state_family
    )
) {
  stop(
    "One or more included cells lack a complete curated-state annotation."
  )
}

# ------------------------------------------------------------
# Add assignments to a full V3 Seurat object
# ------------------------------------------------------------

message("Constructing the full mapped V3 Seurat object...")

full_v3 <- subset(
  full_tcells,
  cells = all_v3_cells
)

rm(
  full_tcells,
  full_query_metadata
)

gc(
  verbose = FALSE
)

metadata_to_add <- all_assignments[
  match(
    colnames(full_v3),
    all_assignments$cell_id
  ),
  setdiff(
    colnames(all_assignments),
    "cell_id"
  ),
  drop = FALSE
]

rownames(metadata_to_add) <- colnames(
  full_v3
)

if (anyNA(rownames(metadata_to_add))) {
  stop(
    "Failed to align assignment metadata to the full V3 object."
  )
}

full_v3 <- AddMetaData(
  object = full_v3,
  metadata = metadata_to_add
)

# ------------------------------------------------------------
# QC summaries
# ------------------------------------------------------------

assignment_origin_counts <- all_assignments %>%
  count(
    state_assignment_origin,
    atlas_inclusion_status,
    name = "n_cells"
  ) %>%
  mutate(
    percent_full_V3 =
      100 * n_cells / nrow(all_assignments)
  )

score_band_counts <- all_assignments %>%
  count(
    transfer_score_band,
    name = "n_cells"
  ) %>%
  mutate(
    percent_full_V3 =
      100 * n_cells / nrow(all_assignments)
  )

assignments_with_sample <- all_assignments %>%
  left_join(
    full_v3@meta.data %>%
      tibble::rownames_to_column(
        "cell_id"
      ) %>%
      select(
        cell_id,
        sample,
        donor,
        group,
        cohort,
        exp_name,
        V_10X
      ),
    by = "cell_id"
  )

transfer_qc_by_sample <-
  assignments_with_sample %>%
  group_by(
    sample,
    donor,
    group,
    cohort
  ) %>%
  summarise(
    total_V3_cells = n(),
    curated_reference_cells = sum(
      state_assignment_origin ==
        "curated_reference"
    ),
    transferred_cells = sum(
      state_assignment_origin ==
        "label_transfer"
    ),
    curation_excluded_cells = sum(
      atlas_inclusion_status ==
        "excluded"
    ),
    median_transfer_score = if (
      any(
        state_assignment_origin ==
          "label_transfer"
      )
    ) {
      median(
        transfer_prediction_score[
          state_assignment_origin ==
            "label_transfer"
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    mean_transfer_score = if (
      any(
        state_assignment_origin ==
          "label_transfer"
      )
    ) {
      mean(
        transfer_prediction_score[
          state_assignment_origin ==
            "label_transfer"
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    pct_transfer_below_050 = if (
      any(
        state_assignment_origin ==
          "label_transfer"
      )
    ) {
      100 * mean(
        transfer_prediction_score[
          state_assignment_origin ==
            "label_transfer"
        ] < 0.50,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    pct_transfer_below_075 = if (
      any(
        state_assignment_origin ==
          "label_transfer"
      )
    ) {
      100 * mean(
        transfer_prediction_score[
          state_assignment_origin ==
            "label_transfer"
        ] < 0.75,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    pct_transfer_below_090 = if (
      any(
        state_assignment_origin ==
          "label_transfer"
      )
    ) {
      100 * mean(
        transfer_prediction_score[
          state_assignment_origin ==
            "label_transfer"
        ] < 0.90,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    tie_repaired_cells = sum(
      cell_received_tie_repair,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

reference_state_counts <- all_assignments %>%
  filter(
    state_assignment_origin ==
      "curated_reference"
  ) %>%
  count(
    curated_immune_state_id,
    name = "reference_cells"
  )

transferred_state_qc <- all_assignments %>%
  filter(
    state_assignment_origin ==
      "label_transfer"
  ) %>%
  group_by(
    curated_immune_state_id
  ) %>%
  summarise(
    transferred_cells = n(),
    mean_transfer_score = mean(
      transfer_prediction_score
    ),
    median_transfer_score = median(
      transfer_prediction_score
    ),
    q25_transfer_score = quantile(
      transfer_prediction_score,
      0.25,
      names = FALSE
    ),
    pct_below_050 = 100 * mean(
      transfer_prediction_score < 0.50
    ),
    pct_below_075 = 100 * mean(
      transfer_prediction_score < 0.75
    ),
    pct_below_090 = 100 * mean(
      transfer_prediction_score < 0.90
    ),
    tie_repaired_cells = sum(
      cell_received_tie_repair,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

transfer_qc_by_state <- atlas_lookup %>%
  select(
    immune_state_id,
    immune_state_name,
    state_family
  ) %>%
  rename(
    curated_immune_state_id =
      immune_state_id
  ) %>%
  left_join(
    reference_state_counts,
    by = "curated_immune_state_id"
  ) %>%
  left_join(
    transferred_state_qc,
    by = "curated_immune_state_id"
  ) %>%
  mutate(
    across(
      c(
        reference_cells,
        transferred_cells,
        tie_repaired_cells
      ),
      ~ replace_na(
        .x,
        0L
      )
    ),
    combined_cells =
      reference_cells +
      transferred_cells
  )

sample_totals <- assignments_with_sample %>%
  group_by(
    sample,
    donor,
    group,
    cohort
  ) %>%
  summarise(
    total_V3_cells = n(),
    included_atlas_cells = sum(
      atlas_inclusion_status ==
        "included"
    ),
    excluded_cells = sum(
      atlas_inclusion_status ==
        "excluded"
    ),
    .groups = "drop"
  )

sample_by_state_composition <-
  assignments_with_sample %>%
  filter(
    atlas_inclusion_status ==
      "included"
  ) %>%
  count(
    sample,
    donor,
    group,
    cohort,
    curated_immune_state_id,
    curated_immune_state_name,
    curated_state_family,
    name = "state_cells"
  ) %>%
  left_join(
    sample_totals,
    by = c(
      "sample",
      "donor",
      "group",
      "cohort"
    )
  ) %>%
  mutate(
    percent_included_atlas_cells =
      100 * state_cells /
      included_atlas_cells,
    percent_total_V3_cells =
      100 * state_cells /
      total_V3_cells
  )

excluded_cell_audit <- assignments_with_sample %>%
  filter(
    atlas_inclusion_status ==
      "excluded"
  ) %>%
  select(
    cell_id,
    sample,
    donor,
    group,
    cohort,
    curation_exclusion_reason,
    transfer_candidate_state_id,
    transfer_prediction_score
  )

final_validation <- data.frame(
  metric = c(
    "chunk_files",
    "expected_query_samples",
    "curated_reference_cells",
    "transferred_prediction_cells",
    "unique_transferred_cell_ids",
    "duplicated_transferred_cell_ids",
    "full_V3_cells",
    "curation_excluded_cells",
    "included_atlas_cells",
    "curated_states",
    "missing_included_state_ids",
    "missing_transfer_scores",
    "tie_repaired_cells",
    "exact_query_cell_set_match",
    "chunk_metadata_match"
  ),
  value = as.character(
    c(
      length(chunk_files),
      EXPECTED_QUERY_SAMPLES,
      nrow(reference_assignments),
      nrow(transferred_predictions),
      length(
        unique(
          transferred_predictions$cell_id
        )
      ),
      sum(
        duplicated(
          transferred_predictions$cell_id
        )
      ),
      length(all_v3_cells),
      nrow(excluded_cell_audit),
      nrow(included_assignments),
      length(
        unique(
          included_assignments$curated_immune_state_id
        )
      ),
      sum(
        is.na(
          included_assignments$curated_immune_state_id
        )
      ),
      sum(
        is.na(
          transferred_predictions$transfer_prediction_score
        )
      ),
      sum(
        all_assignments$cell_received_tie_repair,
        na.rm = TRUE
      ),
      setequal(
        transferred_predictions$cell_id,
        expected_query_cells
      ),
      metadata_match
    )
  )
)

# ------------------------------------------------------------
# Save tables and objects
# ------------------------------------------------------------

message("Saving final QC tables...")

write.csv(
  final_validation,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_final_validation.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  assignment_origin_counts,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_assignment_origin_counts.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  score_band_counts,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_score_band_counts.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  transfer_qc_by_sample,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_QC_by_sample.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  transfer_qc_by_state,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_QC_by_state.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  sample_by_state_composition,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_sample_by_state_composition.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  excluded_cell_audit,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_excluded_cell_audit.csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  all_assignments,
  output_assignments_file,
  compress = TRUE
)

message("Saving the full mapped V3 Seurat object...")

saveRDS(
  full_v3,
  output_object_file,
  compress = TRUE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

message("\nFinal curated V3 atlas validation:")
print(
  final_validation,
  row.names = FALSE
)

message("\nAssignment origins:")
print(
  assignment_origin_counts,
  row.names = FALSE
)

message("\nTransfer score summary:")
print(
  summary(
    transferred_predictions$transfer_prediction_score
  )
)

message("\nCurated-state counts:")
print(
  transfer_qc_by_state %>%
    select(
      curated_immune_state_id,
      reference_cells,
      transferred_cells,
      combined_cells,
      median_transfer_score
    ),
  n = Inf
)

message(
  "\nSaved full mapped V3 object: ",
  output_object_file
)

message(
  "Saved full assignment table: ",
  output_assignments_file
)

message(
  "Curated V3 atlas assembly and validation completed."
)
