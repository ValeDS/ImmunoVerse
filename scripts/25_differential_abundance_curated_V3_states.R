# ============================================================
# 25_differential_abundance_curated_V3_states.R
# Matched differential-abundance analysis of the 11 curated
# ImmunoVerse V3 immune states
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# Paths and reproducibility checks
# ------------------------------------------------------------

FULL_MAPPED_OBJECT <- paste0(
  "data/processed/",
  "GSE144744_Tcells_V3_curated_atlas_full_mapped.rds"
)

RAW_SAMPLE_METADATA <- paste0(
  "data/raw/GSE144744/",
  "GSE144744_metadata_per_sample.csv.gz"
)

ATLAS_TABLE <- paste0(
  "results/tables/",
  "GSE144744_Tcells_V3_Harmony_reference_",
  "curated_immune_state_atlas.csv"
)

OUTPUT_TABLE_DIR <- paste0(
  "results/tables/",
  "curated_V3_differential_abundance"
)

OUTPUT_FIGURE_DIR <- "results/figures"
OUTPUT_MODEL_DIR <- "results/models"

dir.create(
  OUTPUT_TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  FULL_MAPPED_OBJECT,
  RAW_SAMPLE_METADATA,
  ATLAS_TABLE
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

required_packages <- c(
  "Seurat",
  "speckle",
  "limma",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork"
)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package_name) {
      as.character(packageVersion(package_name))
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_package_versions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Load the full mapped atlas and retain compact metadata
# ------------------------------------------------------------

message("Loading the full mapped V3 atlas...")

full_mapped <- readRDS(FULL_MAPPED_OBJECT)

required_metadata <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "curated_immune_state_id",
  "state_assignment_origin",
  "atlas_inclusion_status",
  "transfer_prediction_score"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(full_mapped[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "The full mapped object is missing metadata column(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

cell_metadata <- full_mapped[[]] %>%
  tibble::rownames_to_column("cell_id") %>%
  select(all_of(c("cell_id", required_metadata))) %>%
  mutate(
    across(
      c(
        sample,
        donor,
        group,
        cohort,
        curated_immune_state_id,
        state_assignment_origin,
        atlas_inclusion_status
      ),
      as.character
    ),
    transfer_prediction_score = as.numeric(
      transfer_prediction_score
    )
  )

rm(full_mapped)
invisible(gc())

saveRDS(
  cell_metadata,
  paste0(
    "data/processed/",
    "GSE144744_Tcells_V3_curated_DA_metadata.rds"
  ),
  compress = "xz"
)

# ------------------------------------------------------------
# Load the curated atlas definitions
# ------------------------------------------------------------

atlas_table <- read.csv(
  ATLAS_TABLE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_atlas_columns <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family",
  "confidence"
)

missing_atlas_columns <- setdiff(
  required_atlas_columns,
  colnames(atlas_table)
)

if (length(missing_atlas_columns) > 0L) {
  stop(
    "The curated atlas table is missing column(s): ",
    paste(missing_atlas_columns, collapse = ", ")
  )
}

atlas_definitions <- atlas_table %>%
  transmute(
    immune_state_id = as.character(immune_state_id),
    immune_state_name = as.character(immune_state_name),
    state_family = as.character(state_family),
    annotation_confidence = as.character(confidence)
  ) %>%
  distinct()

state_levels <- atlas_definitions$immune_state_id

if (
  length(state_levels) != 11L ||
  anyDuplicated(state_levels)
) {
  stop(
    "Expected 11 unique curated immune-state identifiers; found ",
    length(unique(state_levels)),
    "."
  )
}

unknown_state_ids <- setdiff(
  unique(
    cell_metadata$curated_immune_state_id[
      cell_metadata$atlas_inclusion_status == "included"
    ]
  ),
  state_levels
)

unknown_state_ids <- unknown_state_ids[
  !is.na(unknown_state_ids)
]

if (length(unknown_state_ids) > 0L) {
  stop(
    "Included cells contain unknown curated state identifiers: ",
    paste(unknown_state_ids, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Reconstruct sample-level matching information
# ------------------------------------------------------------

original_sample_metadata <- read.csv(
  RAW_SAMPLE_METADATA,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_sample_columns <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "match",
  "batch_pair",
  "age_sampling",
  "sex"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(original_sample_metadata)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "The original sample metadata are missing column(s): ",
    paste(missing_sample_columns, collapse = ", ")
  )
}

matching_metadata <- original_sample_metadata %>%
  select(all_of(required_sample_columns)) %>%
  distinct()

if (anyDuplicated(matching_metadata$sample)) {
  duplicated_samples <- matching_metadata$sample[
    duplicated(matching_metadata$sample)
  ]

  stop(
    "Sample-level metadata are not unique for: ",
    paste(unique(duplicated_samples), collapse = ", ")
  )
}

sample_info <- cell_metadata %>%
  distinct(sample, donor, group, cohort) %>%
  left_join(
    matching_metadata,
    by = c(
      "sample",
      "donor",
      "group",
      "cohort"
    )
  )

if (anyDuplicated(sample_info$sample)) {
  stop("Sample identifiers are not unique after metadata reconstruction.")
}

# ------------------------------------------------------------
# Construct primary and sensitivity-analysis cell sets
# ------------------------------------------------------------

included_cells <- cell_metadata %>%
  filter(
    atlas_inclusion_status == "included",
    !is.na(curated_immune_state_id)
  )

analysis_cell_sets <- list(
  full_logit = included_cells,

  full_asin = included_cells,

  confidence_050_logit = included_cells %>%
    filter(
      state_assignment_origin == "curated_reference" |
        (
          state_assignment_origin == "label_transfer" &
            !is.na(transfer_prediction_score) &
            transfer_prediction_score >= 0.50
        )
    ),

  reference_only_logit = included_cells %>%
    filter(
      state_assignment_origin == "curated_reference"
    )
)

analysis_configuration <- tibble::tribble(
  ~analysis_set, ~transformation, ~is_primary,
  "full_logit", "logit", TRUE,
  "full_asin", "asin", FALSE,
  "confidence_050_logit", "logit", FALSE,
  "reference_only_logit", "logit", FALSE
)

analysis_set_summary <- bind_rows(
  lapply(
    names(analysis_cell_sets),
    function(analysis_name) {
      analysis_cell_sets[[analysis_name]] %>%
        summarise(
          analysis_set = analysis_name,
          n_cells = n(),
          n_samples = n_distinct(sample),
          n_donors = n_distinct(donor),
          n_states = n_distinct(curated_immune_state_id)
        )
    }
  )
)

write.csv(
  analysis_set_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_analysis_set_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Matched Propeller analysis
# ------------------------------------------------------------

run_matched_propeller <- function(
    cell_metadata,
    sample_metadata,
    cohort_name,
    reference_group,
    test_group,
    block_column,
    contrast_label,
    analysis_set,
    transform_method = "logit") {

  if (!block_column %in% colnames(sample_metadata)) {
    stop(
      "Block column '",
      block_column,
      "' was not found for ",
      contrast_label,
      "."
    )
  }

  analysis_samples <- sample_metadata %>%
    filter(
      cohort == cohort_name,
      group %in% c(reference_group, test_group)
    ) %>%
    mutate(
      .block = as.character(.data[[block_column]])
    ) %>%
    filter(
      !is.na(.block),
      nzchar(.block)
    ) %>%
    group_by(.block) %>%
    filter(
      n() == 2L,
      n_distinct(sample) == 2L,
      n_distinct(group) == 2L,
      all(c(reference_group, test_group) %in% group)
    ) %>%
    ungroup() %>%
    arrange(.block, factor(
      group,
      levels = c(reference_group, test_group)
    ))

  if (nrow(analysis_samples) == 0L) {
    stop(
      "No complete matched pairs were found for ",
      contrast_label,
      " in ",
      analysis_set,
      "."
    )
  }

  group_counts <- table(analysis_samples$group)

  if (
    length(group_counts) != 2L ||
    unname(group_counts[[reference_group]]) !=
      unname(group_counts[[test_group]])
  ) {
    stop(
      "Unequal group sizes after matched-pair selection for ",
      contrast_label,
      " in ",
      analysis_set,
      "."
    )
  }

  cell_subset <- cell_metadata %>%
    semi_join(
      analysis_samples %>% select(sample),
      by = "sample"
    )

  samples_without_cells <- setdiff(
    analysis_samples$sample,
    unique(cell_subset$sample)
  )

  if (length(samples_without_cells) > 0L) {
    stop(
      "No eligible cells were found for sample(s) ",
      paste(samples_without_cells, collapse = ", "),
      " in ",
      analysis_set,
      "."
    )
  }

  cell_subset <- cell_subset %>%
    mutate(
      curated_immune_state_id = factor(
        curated_immune_state_id,
        levels = state_levels
      )
    )

  if (anyNA(cell_subset$curated_immune_state_id)) {
    stop(
      "Unrecognized state identifiers were encountered in ",
      contrast_label,
      " (",
      analysis_set,
      ")."
    )
  }

  propeller_input <- speckle::getTransformedProps(
    clusters = cell_subset$curated_immune_state_id,
    sample = cell_subset$sample,
    transform = transform_method
  )

  sample_order <- colnames(
    propeller_input$TransformedProps
  )

  design_metadata <- analysis_samples[
    match(sample_order, analysis_samples$sample),
    ,
    drop = FALSE
  ]

  if (
    anyNA(design_metadata$sample) ||
    !identical(sample_order, design_metadata$sample)
  ) {
    stop(
      "The Propeller sample order could not be matched to the ",
      "design metadata for ",
      contrast_label,
      " (",
      analysis_set,
      ")."
    )
  }

  design_metadata$group_factor <- factor(
    design_metadata$group,
    levels = c(reference_group, test_group)
  )

  design_metadata$block_factor <- factor(
    design_metadata$.block
  )

  design <- model.matrix(
    ~ 0 + group_factor + block_factor,
    data = design_metadata
  )

  colnames(design) <- sub(
    "^group_factor",
    "group",
    colnames(design)
  )

  colnames(design) <- sub(
    "^block_factor",
    "block",
    colnames(design)
  )

  colnames(design) <- make.names(
    colnames(design),
    unique = TRUE
  )

  if (qr(design)$rank != ncol(design)) {
    stop(
      "The design matrix is not full rank for ",
      contrast_label,
      " (",
      analysis_set,
      ")."
    )
  }

  contrast_expression <- paste0(
    "group",
    make.names(test_group),
    "-group",
    make.names(reference_group)
  )

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_expression,
    levels = design
  )

  propeller_result <- speckle::propeller.ttest(
    prop.list = propeller_input,
    design = design,
    contrasts = contrast_matrix,
    robust = TRUE,
    trend = FALSE,
    sort = FALSE
  )

  propeller_result <- as.data.frame(
    propeller_result,
    stringsAsFactors = FALSE
  )

  propeller_result$immune_state_id <- rownames(
    propeller_result
  )

  required_result_columns <- c("P.Value", "FDR")

  if (
    !all(
      required_result_columns %in%
        colnames(propeller_result)
    )
  ) {
    stop(
      "Propeller did not return the expected P.Value and FDR ",
      "columns for ",
      contrast_label,
      " (",
      analysis_set,
      ")."
    )
  }

  proportion_matrix <- propeller_input$Proportions

  reference_columns <- (
    design_metadata$group == reference_group
  )

  test_columns <- (
    design_metadata$group == test_group
  )

  effect_sizes <- data.frame(
    immune_state_id = rownames(proportion_matrix),
    mean_reference_percent = 100 * rowMeans(
      proportion_matrix[
        ,
        reference_columns,
        drop = FALSE
      ]
    ),
    mean_test_percent = 100 * rowMeans(
      proportion_matrix[
        ,
        test_columns,
        drop = FALSE
      ]
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      difference_percentage_points =
        mean_test_percent - mean_reference_percent,
      proportion_ratio = ifelse(
        mean_reference_percent > 0,
        mean_test_percent / mean_reference_percent,
        NA_real_
      )
    )

  result <- effect_sizes %>%
    left_join(
      propeller_result,
      by = "immune_state_id"
    ) %>%
    left_join(
      atlas_definitions,
      by = "immune_state_id"
    ) %>%
    mutate(
      contrast = contrast_label,
      cohort = cohort_name,
      reference_group = reference_group,
      test_group = test_group,
      block_variable = block_column,
      n_pairs = nrow(analysis_samples) / 2,
      transformation = transform_method,
      analysis_set = analysis_set,
      significance_class = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        P.Value < 0.05 ~ "Nominal P < 0.05",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(
      factor(immune_state_id, levels = state_levels)
    )

  if (nrow(result) != length(state_levels)) {
    stop(
      "Expected ",
      length(state_levels),
      " tested states but obtained ",
      nrow(result),
      " for ",
      contrast_label,
      " (",
      analysis_set,
      ")."
    )
  }

  # `Proportions` is a table rather than a plain matrix. Converting
  # its transpose directly with as.data.frame() invokes the table
  # method and already creates a `sample` column, which would clash
  # with rownames_to_column(). Convert the table directly to its
  # canonical long representation instead.
  sample_proportions <- as.data.frame(
    proportion_matrix,
    stringsAsFactors = FALSE,
    responseName = "proportion"
  ) %>%
    rename(
      immune_state_id = clusters
    ) %>%
    mutate(
      sample = as.character(sample),
      immune_state_id = as.character(
        immune_state_id
      )
    ) %>%
    mutate(
      percent = 100 * proportion
    ) %>%
    left_join(
      design_metadata %>%
        select(
          sample,
          donor,
          group,
          cohort,
          .block
        ),
      by = "sample"
    )

  list(
    results = result,
    sample_proportions = sample_proportions,
    design = design,
    design_metadata = design_metadata,
    contrast_matrix = contrast_matrix,
    transformed_proportions =
      propeller_input$TransformedProps
  )
}

# ------------------------------------------------------------
# Define the four biologically valid within-cohort contrasts
# ------------------------------------------------------------

contrast_definitions <- tibble::tribble(
  ~contrast, ~cohort, ~reference_group, ~test_group, ~block_column,
  "MS1_vs_HI1", "NAT_HI", "HI1", "MS1", "batch_pair",
  "MS1nat_vs_MS1", "NAT_HI", "MS1", "MS1_nat", "donor",
  "MS2_vs_HI2", "RRMS_HI", "HI2", "MS2", "batch_pair",
  "PPMS_vs_HI3", "PPMS_HI", "HI3", "PPMS", "batch_pair"
)

write.csv(
  contrast_definitions,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_contrast_definitions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Run primary and sensitivity analyses
# ------------------------------------------------------------

all_models <- list()
all_results <- list()
all_primary_sample_proportions <- list()

for (
  configuration_index in
  seq_len(nrow(analysis_configuration))
) {
  configuration <- analysis_configuration[
    configuration_index,
    ,
    drop = FALSE
  ]

  analysis_name <- configuration$analysis_set[[1]]
  transformation_name <- configuration$transformation[[1]]

  message(
    "\nAnalysis set: ",
    analysis_name,
    " (",
    transformation_name,
    ")"
  )

  analysis_cells <- analysis_cell_sets[[analysis_name]]

  for (
    contrast_index in
    seq_len(nrow(contrast_definitions))
  ) {
    contrast_definition <- contrast_definitions[
      contrast_index,
      ,
      drop = FALSE
    ]

    contrast_name <- contrast_definition$contrast[[1]]

    message("  Running ", contrast_name, "...")

    model_key <- paste(
      analysis_name,
      contrast_name,
      sep = "__"
    )

    fitted_analysis <- run_matched_propeller(
      cell_metadata = analysis_cells,
      sample_metadata = sample_info,
      cohort_name =
        contrast_definition$cohort[[1]],
      reference_group =
        contrast_definition$reference_group[[1]],
      test_group =
        contrast_definition$test_group[[1]],
      block_column =
        contrast_definition$block_column[[1]],
      contrast_label = contrast_name,
      analysis_set = analysis_name,
      transform_method = transformation_name
    )

    all_models[[model_key]] <- fitted_analysis[
      c(
        "design",
        "design_metadata",
        "contrast_matrix",
        "transformed_proportions"
      )
    ]

    all_results[[model_key]] <- fitted_analysis$results

    if (isTRUE(configuration$is_primary[[1]])) {
      all_primary_sample_proportions[[contrast_name]] <-
        fitted_analysis$sample_proportions %>%
        mutate(
          contrast = contrast_name
        )
    }
  }
}

combined_results <- bind_rows(all_results)

primary_results <- combined_results %>%
  filter(analysis_set == "full_logit") %>%
  mutate(
    immune_state_id = factor(
      immune_state_id,
      levels = state_levels
    ),
    contrast = factor(
      contrast,
      levels = contrast_definitions$contrast
    )
  ) %>%
  arrange(contrast, immune_state_id) %>%
  mutate(
    immune_state_id = as.character(
      immune_state_id
    ),
    contrast = as.character(contrast)
  )

primary_sample_proportions <- bind_rows(
  all_primary_sample_proportions
)

# ------------------------------------------------------------
# Summaries and sensitivity diagnostics
# ------------------------------------------------------------

significance_summary <- primary_results %>%
  group_by(
    contrast,
    cohort,
    reference_group,
    test_group,
    block_variable,
    n_pairs
  ) %>%
  summarise(
    tested_states = n(),
    nominal_p_below_005 = sum(
      P.Value < 0.05,
      na.rm = TRUE
    ),
    FDR_below_005 = sum(
      FDR < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

top_primary_results <- primary_results %>%
  group_by(contrast) %>%
  arrange(P.Value, FDR, .by_group = TRUE) %>%
  slice_head(n = 5L) %>%
  ungroup() %>%
  select(
    contrast,
    immune_state_id,
    immune_state_name,
    state_family,
    mean_reference_percent,
    mean_test_percent,
    difference_percentage_points,
    proportion_ratio,
    P.Value,
    FDR
  )

primary_effects <- primary_results %>%
  select(
    contrast,
    immune_state_id,
    primary_difference_percentage_points =
      difference_percentage_points
  )

sensitivity_comparison_long <- combined_results %>%
  filter(analysis_set != "full_logit") %>%
  select(
    analysis_set,
    contrast,
    immune_state_id,
    sensitivity_difference_percentage_points =
      difference_percentage_points,
    sensitivity_P_Value = P.Value,
    sensitivity_FDR = FDR
  ) %>%
  left_join(
    primary_effects,
    by = c("contrast", "immune_state_id")
  ) %>%
  mutate(
    direction_agreement = case_when(
      primary_difference_percentage_points == 0 &
        sensitivity_difference_percentage_points == 0 ~ TRUE,
      primary_difference_percentage_points == 0 |
        sensitivity_difference_percentage_points == 0 ~ FALSE,
      TRUE ~ sign(
        primary_difference_percentage_points
      ) == sign(
        sensitivity_difference_percentage_points
      )
    )
  )

safe_correlation <- function(x, y, method) {
  complete <- complete.cases(x, y)

  if (
    sum(complete) < 3L ||
    stats::sd(x[complete]) == 0 ||
    stats::sd(y[complete]) == 0
  ) {
    return(NA_real_)
  }

  stats::cor(
    x[complete],
    y[complete],
    method = method
  )
}

sensitivity_summary <- sensitivity_comparison_long %>%
  group_by(analysis_set, contrast) %>%
  summarise(
    tested_states = n(),
    pearson_effect_correlation = safe_correlation(
      primary_difference_percentage_points,
      sensitivity_difference_percentage_points,
      method = "pearson"
    ),
    spearman_effect_correlation = safe_correlation(
      primary_difference_percentage_points,
      sensitivity_difference_percentage_points,
      method = "spearman"
    ),
    direction_agreement_percent = 100 * mean(
      direction_agreement,
      na.rm = TRUE
    ),
    nominal_p_below_005 = sum(
      sensitivity_P_Value < 0.05,
      na.rm = TRUE
    ),
    FDR_below_005 = sum(
      sensitivity_FDR < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Save complete results and fitted objects
# ------------------------------------------------------------

write.csv(
  primary_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_primary_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_all_analyses.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_sample_proportions,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_primary_sample_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  significance_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_significance_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_primary_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_top_primary_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_comparison_long,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_sensitivity_state_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    models = all_models,
    contrast_definitions = contrast_definitions,
    analysis_configuration = analysis_configuration,
    package_versions = package_versions,
    session_info = sessionInfo()
  ),
  file.path(
    OUTPUT_MODEL_DIR,
    "GSE144744_curated_V3_DA_propeller_models.rds"
  ),
  compress = "xz"
)

# ------------------------------------------------------------
# Publication-oriented summary figure
# ------------------------------------------------------------

contrast_labels <- c(
  MS1_vs_HI1 = "Untreated RRMS\nvs matched HI",
  MS1nat_vs_MS1 = "Natalizumab\nvs pretreatment",
  MS2_vs_HI2 = "RRMS\nvs matched HI",
  PPMS_vs_HI3 = "PPMS\nvs matched HI"
)

analysis_labels <- c(
  full_asin = "Arcsin transformation",
  confidence_050_logit = "Transferred score >= 0.50",
  reference_only_logit = "Reference cells only"
)

plot_primary <- primary_results %>%
  mutate(
    immune_state_id = factor(
      immune_state_id,
      levels = rev(state_levels)
    ),
    contrast = factor(
      contrast,
      levels = contrast_definitions$contrast,
      labels = contrast_labels[
        contrast_definitions$contrast
      ]
    ),
    effect_label = paste0(
      sprintf(
        "%+.2f",
        difference_percentage_points
      ),
      case_when(
        FDR < 0.05 ~ "**",
        P.Value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  )

effect_limit <- max(
  abs(
    plot_primary$difference_percentage_points
  ),
  na.rm = TRUE
)

if (!is.finite(effect_limit) || effect_limit == 0) {
  effect_limit <- 1
}

panel_A <- ggplot(
  plot_primary,
  aes(
    x = contrast,
    y = immune_state_id,
    fill = difference_percentage_points
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = effect_label),
    size = 3.6,
    color = "#17202A"
  ) +
  scale_fill_gradient2(
    low = "#3C78B5",
    mid = "#F7F7F7",
    high = "#D64C4C",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    name = "Difference\n(percentage points)"
  ) +
  labs(
    title = "Curated immune-state differential abundance",
    subtitle = paste0(
      "Positive values indicate enrichment in the disease or ",
      "treatment group; * nominal P < 0.05, ** FDR < 0.05"
    ),
    x = NULL,
    y = "Curated immune state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      color = "#5B5B5B",
      size = 10.5
    )
  )

plot_sensitivity <- sensitivity_summary %>%
  mutate(
    analysis_set = factor(
      analysis_set,
      levels = names(analysis_labels),
      labels = unname(analysis_labels)
    ),
    contrast = factor(
      contrast,
      levels = contrast_definitions$contrast,
      labels = contrast_labels[
        contrast_definitions$contrast
      ]
    ),
    sensitivity_label = paste0(
      "r = ",
      ifelse(
        is.na(pearson_effect_correlation),
        "NA",
        sprintf(
          "%.2f",
          pearson_effect_correlation
        )
      ),
      "\n",
      sprintf(
        "%.0f%% direction",
        direction_agreement_percent
      )
    )
  )

panel_B <- ggplot(
  plot_sensitivity,
  aes(
    x = contrast,
    y = analysis_set,
    fill = pearson_effect_correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sensitivity_label),
    size = 3.4,
    lineheight = 0.92,
    color = "#17202A"
  ) +
  scale_fill_gradient(
    low = "#F2F6FA",
    high = "#2E6FA8",
    limits = c(0, 1),
    na.value = "#E5E5E5",
    name = "Pearson r"
  ) +
  labs(
    title = "Sensitivity of abundance-effect estimates",
    subtitle = paste0(
      "Agreement with the primary full-atlas logit analysis"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      color = "#5B5B5B",
      size = 10.5
    )
  )

combined_figure <- (
  panel_A /
    panel_B
) +
  plot_layout(
    heights = c(1.8, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 17
      )
    )
  )

figure_pdf <- file.path(
  OUTPUT_FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_V3_curated_state_",
    "differential_abundance.pdf"
  )
)

figure_png <- file.path(
  OUTPUT_FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_V3_curated_state_",
    "differential_abundance.png"
  )
)

ggsave(
  filename = figure_pdf,
  plot = combined_figure,
  width = 13.5,
  height = 11,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = figure_png,
  plot = combined_figure,
  width = 13.5,
  height = 11,
  units = "in",
  dpi = 320,
  bg = "white"
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

message("\nPrimary differential-abundance summary:")
print(as.data.frame(significance_summary))

message("\nTop primary results:")
print(as.data.frame(top_primary_results))

message("\nSensitivity summary:")
print(as.data.frame(sensitivity_summary))

message("\nSaved primary result table: ")
message(
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_V3_DA_primary_results.csv"
  )
)

message("Saved summary PDF: ", figure_pdf)
message("Saved summary PNG: ", figure_png)
message(
  "Curated V3 immune-state differential-abundance ",
  "analysis completed."
)
