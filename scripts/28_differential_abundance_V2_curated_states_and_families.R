# ============================================================
# 28_differential_abundance_V2_curated_states_and_families.R
# Matched V2 differential-abundance analysis:
#   primary: six portable higher-order immune-state families
#   exploratory: eleven transferred curated immune states
# ============================================================

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

V2_PREDICTION_FILE <- paste0(
  "data/processed/",
  "GSE144744_Tcells_V3_to_V2_",
  "curated_state_predictions.rds"
)

RAW_SAMPLE_METADATA <- paste0(
  "data/raw/GSE144744/",
  "GSE144744_metadata_per_sample.csv.gz"
)

ATLAS_TABLE_FILE <- paste0(
  "results/tables/",
  "GSE144744_Tcells_V3_Harmony_reference_",
  "curated_immune_state_atlas.csv"
)

OUTPUT_TABLE_DIR <- paste0(
  "results/tables/",
  "V2_curated_differential_abundance"
)

OUTPUT_MODEL_DIR <- "results/models"
OUTPUT_FIGURE_DIR <- "results/figures"

dir.create(
  OUTPUT_TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  V2_PREDICTION_FILE,
  RAW_SAMPLE_METADATA,
  ATLAS_TABLE_FILE
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

# ------------------------------------------------------------
# Load V2 predictions and atlas definitions
# ------------------------------------------------------------

message("Loading cross-chemistry V2 predictions...")
v2_predictions <- readRDS(V2_PREDICTION_FILE)

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

missing_prediction_columns <- setdiff(
  required_prediction_columns,
  colnames(v2_predictions)
)

if (length(missing_prediction_columns) > 0L) {
  stop(
    "The V2 prediction table is missing column(s): ",
    paste(missing_prediction_columns, collapse = ", ")
  )
}

atlas_table <- read.csv(
  ATLAS_TABLE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_atlas_columns <- c(
  "immune_state_id",
  "immune_state_name",
  "state_family"
)

missing_atlas_columns <- setdiff(
  required_atlas_columns,
  colnames(atlas_table)
)

if (length(missing_atlas_columns) > 0L) {
  stop(
    "The atlas table is missing column(s): ",
    paste(missing_atlas_columns, collapse = ", ")
  )
}

atlas_definitions <- atlas_table %>%
  transmute(
    immune_state_id = as.character(immune_state_id),
    immune_state_name = as.character(immune_state_name),
    original_state_family = as.character(state_family)
  ) %>%
  distinct()

state_levels <- atlas_definitions$immune_state_id

if (
  length(state_levels) != 11L ||
    anyDuplicated(state_levels)
) {
  stop("Expected 11 unique curated immune states.")
}

higher_family_definitions <- data.frame(
  immune_state_id = state_levels,
  higher_order_family = c(
    rep("CD4 T-cell states", 4),
    "Cytotoxic/effector CD8",
    "Naive/memory CD8",
    "Cytotoxic/effector CD8",
    "Innate-like cytotoxic T",
    "Regulatory T",
    "Innate-like cytotoxic T",
    "Interferon-responsive T"
  ),
  stringsAsFactors = FALSE
)

family_levels <- c(
  "CD4 T-cell states",
  "Cytotoxic/effector CD8",
  "Naive/memory CD8",
  "Innate-like cytotoxic T",
  "Regulatory T",
  "Interferon-responsive T"
)

v2_metadata <- v2_predictions %>%
  transmute(
    cell_id = as.character(cell_id),
    sample = as.character(sample),
    donor = as.character(donor),
    group = as.character(group),
    cohort = as.character(cohort),
    V_10X = as.character(V_10X),
    immune_state_id = as.character(predicted.id),
    prediction_score = as.numeric(
      prediction.score.max
    )
  ) %>%
  left_join(
    higher_family_definitions,
    by = "immune_state_id"
  )

if (
  anyNA(v2_metadata$immune_state_id) ||
    anyNA(v2_metadata$higher_order_family) ||
    anyNA(v2_metadata$prediction_score)
) {
  stop(
    "Missing state, family, or prediction-score values were ",
    "found in the V2 metadata."
  )
}

if (
  any(v2_metadata$V_10X != "V2") ||
    !all(v2_metadata$cohort == "RRMS_HI") ||
    !all(v2_metadata$group %in% c("HI2", "MS2"))
) {
  stop(
    "The V2 prediction table does not contain the expected ",
    "RRMS_HI HI2/MS2 design."
  )
}

# ------------------------------------------------------------
# Reconstruct the matched-pair design
# ------------------------------------------------------------

raw_sample_metadata <- read.csv(
  RAW_SAMPLE_METADATA,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_sample_columns <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "batch_pair",
  "V_10X",
  "age_sampling",
  "sex"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(raw_sample_metadata)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "The raw sample metadata are missing column(s): ",
    paste(missing_sample_columns, collapse = ", ")
  )
}

v2_sample_metadata <- raw_sample_metadata %>%
  filter(V_10X == "V2") %>%
  select(all_of(required_sample_columns)) %>%
  distinct()

if (anyDuplicated(v2_sample_metadata$sample)) {
  stop(
    "The V2 sample metadata are not unique by sample."
  )
}

observed_samples <- v2_metadata %>%
  distinct(
    sample,
    donor,
    group,
    cohort,
    V_10X
  )

sample_info <- observed_samples %>%
  left_join(
    v2_sample_metadata,
    by = c(
      "sample",
      "donor",
      "group",
      "cohort",
      "V_10X"
    )
  )

if (
  anyNA(sample_info$batch_pair) ||
    any(sample_info$batch_pair == "")
) {
  stop(
    "One or more V2 samples lack a matched batch_pair."
  )
}

paired_samples <- sample_info %>%
  mutate(
    batch_pair = as.character(batch_pair)
  ) %>%
  group_by(batch_pair) %>%
  filter(
    n() == 2L,
    n_distinct(sample) == 2L,
    n_distinct(group) == 2L,
    all(c("HI2", "MS2") %in% group)
  ) %>%
  ungroup() %>%
  arrange(
    batch_pair,
    factor(group, levels = c("HI2", "MS2"))
  )

if (
  nrow(paired_samples) != 22L ||
    n_distinct(paired_samples$batch_pair) != 11L ||
    sum(paired_samples$group == "HI2") != 11L ||
    sum(paired_samples$group == "MS2") != 11L
) {
  stop(
    "Expected 11 complete V2 HI2/MS2 pairs; found ",
    n_distinct(paired_samples$batch_pair),
    "."
  )
}

write.csv(
  paired_samples,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_matched_sample_design.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Analysis sets
# ------------------------------------------------------------

analysis_sets <- list(
  all_cells = v2_metadata,
  score_050 = v2_metadata %>%
    filter(prediction_score >= 0.50)
)

analysis_set_summary <- bind_rows(
  lapply(
    names(analysis_sets),
    function(analysis_name) {
      analysis_sets[[analysis_name]] %>%
        summarise(
          analysis_set = analysis_name,
          n_cells = n(),
          n_samples = n_distinct(sample),
          n_states = n_distinct(immune_state_id),
          n_families = n_distinct(
            higher_order_family
          )
        )
    }
  )
)

write.csv(
  analysis_set_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_analysis_set_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Matched Propeller function
# ------------------------------------------------------------

run_v2_propeller <- function(
    cell_metadata,
    label_column,
    label_levels,
    resolution,
    analysis_set,
    transformation = "logit") {

  cell_subset <- cell_metadata %>%
    semi_join(
      paired_samples %>% select(sample),
      by = "sample"
    )

  samples_without_cells <- setdiff(
    paired_samples$sample,
    unique(cell_subset$sample)
  )

  if (length(samples_without_cells) > 0L) {
    stop(
      "No eligible cells remain for sample(s): ",
      paste(samples_without_cells, collapse = ", "),
      " in ",
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  cluster_labels <- factor(
    cell_subset[[label_column]],
    levels = label_levels
  )

  if (anyNA(cluster_labels)) {
    stop(
      "Unrecognized labels were found in ",
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  propeller_input <- speckle::getTransformedProps(
    clusters = cluster_labels,
    sample = cell_subset$sample,
    transform = transformation
  )

  sample_order <- colnames(
    propeller_input$TransformedProps
  )

  design_metadata <- paired_samples[
    match(sample_order, paired_samples$sample),
    ,
    drop = FALSE
  ]

  if (
    anyNA(design_metadata$sample) ||
      !identical(
        sample_order,
        design_metadata$sample
      )
  ) {
    stop(
      "Unable to match the Propeller sample order for ",
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  design_metadata$group_factor <- factor(
    design_metadata$group,
    levels = c("HI2", "MS2")
  )

  design_metadata$block_factor <- factor(
    design_metadata$batch_pair
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
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  contrast_matrix <- limma::makeContrasts(
    contrasts = "groupMS2-groupHI2",
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

  propeller_result$label <- rownames(
    propeller_result
  )

  if (
    !all(c("P.Value", "FDR") %in%
      colnames(propeller_result))
  ) {
    stop(
      "Propeller did not return P.Value and FDR for ",
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  proportion_matrix <- propeller_input$Proportions

  hi_columns <- design_metadata$group == "HI2"
  ms_columns <- design_metadata$group == "MS2"

  effects <- data.frame(
    label = rownames(proportion_matrix),
    mean_HI2_percent = 100 * rowMeans(
      proportion_matrix[
        ,
        hi_columns,
        drop = FALSE
      ]
    ),
    mean_MS2_percent = 100 * rowMeans(
      proportion_matrix[
        ,
        ms_columns,
        drop = FALSE
      ]
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      difference_percentage_points =
        mean_MS2_percent - mean_HI2_percent,
      proportion_ratio = ifelse(
        mean_HI2_percent > 0,
        mean_MS2_percent / mean_HI2_percent,
        NA_real_
      )
    )

  results <- effects %>%
    left_join(
      propeller_result,
      by = "label"
    ) %>%
    mutate(
      resolution = resolution,
      analysis_set = analysis_set,
      transformation = transformation,
      contrast = "MS2_vs_HI2",
      cohort = "RRMS_HI",
      n_pairs = 11L,
      significance_class = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        P.Value < 0.05 ~ "Nominal P < 0.05",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(
      factor(label, levels = label_levels)
    )

  if (nrow(results) != length(label_levels)) {
    stop(
      "Expected ",
      length(label_levels),
      " labels but obtained ",
      nrow(results),
      " in ",
      resolution,
      "/",
      analysis_set,
      "."
    )
  }

  sample_proportions <- as.data.frame(
    proportion_matrix,
    stringsAsFactors = FALSE,
    responseName = "proportion"
  ) %>%
    rename(label = clusters) %>%
    mutate(
      sample = as.character(sample),
      label = as.character(label),
      percent = 100 * proportion
    ) %>%
    left_join(
      design_metadata %>%
        select(
          sample,
          donor,
          group,
          cohort,
          batch_pair
        ),
      by = "sample"
    ) %>%
    mutate(
      resolution = resolution,
      analysis_set = analysis_set,
      transformation = transformation
    )

  list(
    results = results,
    sample_proportions = sample_proportions,
    design = design,
    design_metadata = design_metadata,
    contrast_matrix = contrast_matrix,
    transformed_proportions =
      propeller_input$TransformedProps
  )
}

# ------------------------------------------------------------
# Run primary and sensitivity analyses
# ------------------------------------------------------------

analysis_configuration <- tibble::tribble(
  ~resolution, ~analysis_set, ~transformation, ~is_primary,
  "family", "all_cells", "logit", TRUE,
  "family", "all_cells", "asin", FALSE,
  "family", "score_050", "logit", FALSE,
  "state", "all_cells", "logit", TRUE,
  "state", "all_cells", "asin", FALSE,
  "state", "score_050", "logit", FALSE
)

all_results <- list()
all_models <- list()
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

  resolution_name <- configuration$resolution[[1]]
  analysis_name <- configuration$analysis_set[[1]]
  transformation_name <-
    configuration$transformation[[1]]

  if (resolution_name == "family") {
    label_column <- "higher_order_family"
    label_levels <- family_levels
  } else {
    label_column <- "immune_state_id"
    label_levels <- state_levels
  }

  model_key <- paste(
    resolution_name,
    analysis_name,
    transformation_name,
    sep = "__"
  )

  message(
    "Running ",
    model_key,
    "..."
  )

  fitted_analysis <- run_v2_propeller(
    cell_metadata =
      analysis_sets[[analysis_name]],
    label_column = label_column,
    label_levels = label_levels,
    resolution = resolution_name,
    analysis_set = analysis_name,
    transformation = transformation_name
  )

  all_results[[model_key]] <-
    fitted_analysis$results

  all_models[[model_key]] <- fitted_analysis[
    c(
      "design",
      "design_metadata",
      "contrast_matrix",
      "transformed_proportions"
    )
  ]

  if (isTRUE(configuration$is_primary[[1]])) {
    all_primary_sample_proportions[[
      resolution_name
    ]] <- fitted_analysis$sample_proportions
  }
}

combined_results <- bind_rows(all_results)

primary_results <- combined_results %>%
  filter(
    analysis_set == "all_cells",
    transformation == "logit"
  )

primary_family_results <- primary_results %>%
  filter(resolution == "family") %>%
  mutate(
    label = factor(
      label,
      levels = family_levels
    )
  ) %>%
  arrange(label) %>%
  mutate(label = as.character(label))

primary_state_results <- primary_results %>%
  filter(resolution == "state") %>%
  left_join(
    atlas_definitions,
    by = c("label" = "immune_state_id")
  ) %>%
  mutate(
    label = factor(
      label,
      levels = state_levels
    )
  ) %>%
  arrange(label) %>%
  mutate(label = as.character(label))

primary_sample_proportions <- bind_rows(
  all_primary_sample_proportions
)

# ------------------------------------------------------------
# Sensitivity summaries
# ------------------------------------------------------------

primary_for_sensitivity <- primary_results %>%
  transmute(
    resolution,
    label,
    primary_effect =
      difference_percentage_points,
    primary_P_Value = P.Value,
    primary_FDR = FDR,
    primary_signed_evidence =
      sign(difference_percentage_points) *
      -log10(pmax(P.Value, 1e-300)),
    primary_nominal = P.Value < 0.05
  )

sensitivity_state_comparison <- combined_results %>%
  filter(
    !(analysis_set == "all_cells" &
      transformation == "logit")
  ) %>%
  transmute(
    resolution,
    analysis_set,
    transformation,
    label,
    sensitivity_effect =
      difference_percentage_points,
    sensitivity_P_Value = P.Value,
    sensitivity_FDR = FDR,
    sensitivity_signed_evidence =
      sign(difference_percentage_points) *
      -log10(pmax(P.Value, 1e-300)),
    sensitivity_nominal = P.Value < 0.05
  ) %>%
  left_join(
    primary_for_sensitivity,
    by = c("resolution", "label")
  ) %>%
  mutate(
    sensitivity_label = case_when(
      analysis_set == "all_cells" &
        transformation == "asin" ~
        "Arcsin transformation",
      analysis_set == "score_050" &
        transformation == "logit" ~
        "Prediction score >= 0.50",
      TRUE ~ paste(
        analysis_set,
        transformation
      )
    ),
    direction_agreement = case_when(
      primary_effect == 0 &
        sensitivity_effect == 0 ~ TRUE,
      primary_effect == 0 |
        sensitivity_effect == 0 ~ FALSE,
      TRUE ~ sign(primary_effect) ==
        sign(sensitivity_effect)
    ),
    nominal_agreement = (
      primary_nominal ==
        sensitivity_nominal
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

sensitivity_summary <- sensitivity_state_comparison %>%
  group_by(
    resolution,
    sensitivity_label
  ) %>%
  summarise(
    tested_labels = n(),
    effect_pearson = safe_correlation(
      primary_effect,
      sensitivity_effect,
      method = "pearson"
    ),
    signed_evidence_spearman =
      safe_correlation(
        primary_signed_evidence,
        sensitivity_signed_evidence,
        method = "spearman"
      ),
    direction_agreement_percent =
      100 * mean(
        direction_agreement,
        na.rm = TRUE
      ),
    nominal_agreement_percent =
      100 * mean(
        nominal_agreement,
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

significance_summary <- primary_results %>%
  group_by(resolution) %>%
  summarise(
    tested_labels = n(),
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

# ------------------------------------------------------------
# Save tables and fitted objects
# ------------------------------------------------------------

write.csv(
  primary_family_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_primary_family_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_state_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_exploratory_state_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_results,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_all_analyses.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_sample_proportions,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_primary_sample_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  significance_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_significance_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_state_comparison,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_sensitivity_label_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_DA_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    models = all_models,
    analysis_configuration =
      analysis_configuration,
    paired_samples = paired_samples,
    higher_family_definitions =
      higher_family_definitions,
    session_info = sessionInfo()
  ),
  file.path(
    OUTPUT_MODEL_DIR,
    "GSE144744_V2_curated_DA_propeller_models.rds"
  ),
  compress = "xz"
)

# ------------------------------------------------------------
# Publication-oriented summary figure
# ------------------------------------------------------------

family_palette <- c(
  "CD4 T-cell states" = "#4C78A8",
  "Cytotoxic/effector CD8" = "#E15759",
  "Naive/memory CD8" = "#76B7B2",
  "Innate-like cytotoxic T" = "#E7B416",
  "Regulatory T" = "#B279A2",
  "Interferon-responsive T" = "#2F4858"
)

family_sample_plot <-
  all_primary_sample_proportions$family %>%
  mutate(
    label = factor(
      label,
      levels = family_levels
    ),
    group = factor(
      group,
      levels = c("HI2", "MS2")
    )
  )

panel_A <- ggplot(
  family_sample_plot,
  aes(
    x = group,
    y = percent,
    group = batch_pair
  )
) +
  geom_line(
    color = "#A6A6A6",
    linewidth = 0.45,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = group),
    size = 1.8,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ label,
    scales = "free_y",
    ncol = 3
  ) +
  scale_color_manual(
    values = c(
      HI2 = "#4C78A8",
      MS2 = "#D64C4C"
    ),
    name = NULL
  ) +
  labs(
    title = "Paired V2 family abundances",
    subtitle = "Each line connects one matched HI2-MS2 pair",
    x = NULL,
    y = "T cells (%)"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(
      fill = "#F2F2F2",
      color = NA
    ),
    strip.text = element_text(
      face = "bold",
      size = 9.5
    ),
    legend.position = "top",
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

family_effect_plot <- primary_family_results %>%
  mutate(
    label = factor(
      label,
      levels = rev(family_levels)
    ),
    effect_label = paste0(
      sprintf(
        "%+.2f pp",
        difference_percentage_points
      ),
      case_when(
        FDR < 0.05 ~ " **",
        P.Value < 0.05 ~ " *",
        TRUE ~ ""
      )
    )
  )

panel_B <- ggplot(
  family_effect_plot,
  aes(
    x = difference_percentage_points,
    y = label,
    color = as.character(label)
  )
) +
  geom_vline(
    xintercept = 0,
    color = "#777777",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = difference_percentage_points,
      yend = label
    ),
    linewidth = 1
  ) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = effect_label),
    color = "#222222",
    hjust = ifelse(
      family_effect_plot$
        difference_percentage_points >= 0,
      -0.12,
      1.12
    ),
    size = 3.5
  ) +
  scale_color_manual(
    values = family_palette,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.25, 0.25))
  ) +
  labs(
    title = "Primary family-level differential abundance",
    subtitle = paste0(
      "Positive values indicate enrichment in MS2; ",
      "* nominal P < 0.05, ** FDR < 0.05"
    ),
    x = "MS2 - HI2 (percentage points)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

state_effect_plot <- primary_state_results %>%
  mutate(
    label = factor(
      label,
      levels = rev(state_levels)
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

panel_C <- ggplot(
  state_effect_plot,
  aes(
    x = difference_percentage_points,
    y = label,
    fill = difference_percentage_points
  )
) +
  geom_vline(
    xintercept = 0,
    color = "#777777",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_col(
    width = 0.68,
    color = "white"
  ) +
  geom_text(
    aes(label = effect_label),
    hjust = ifelse(
      state_effect_plot$
        difference_percentage_points >= 0,
      -0.12,
      1.12
    ),
    size = 3.2,
    color = "#222222"
  ) +
  scale_fill_gradient2(
    low = "#3C78B5",
    mid = "#F7F7F7",
    high = "#D64C4C",
    midpoint = 0,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.20, 0.20))
  ) +
  labs(
    title = "Exploratory state-level effects",
    subtitle = paste0(
      "Fine-grained states are less portable across chemistry"
    ),
    x = "MS2 - HI2 (percentage points)",
    y = "Curated immune state"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

sensitivity_plot <- sensitivity_summary %>%
  mutate(
    resolution = recode(
      resolution,
      family = "Higher-order families",
      state = "Curated states"
    ),
    sensitivity_label = factor(
      sensitivity_label,
      levels = c(
        "Prediction score >= 0.50",
        "Arcsin transformation"
      )
    ),
    tile_label = paste0(
      "rho evidence = ",
      sprintf(
        "%.2f",
        signed_evidence_spearman
      ),
      "\n",
      sprintf(
        "%.0f%% direction",
        direction_agreement_percent
      ),
      "; ",
      sprintf(
        "%.0f%% P-status",
        nominal_agreement_percent
      )
    )
  )

panel_D <- ggplot(
  sensitivity_plot,
  aes(
    x = sensitivity_label,
    y = resolution,
    fill = signed_evidence_spearman
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = tile_label),
    size = 3.3,
    lineheight = 0.92,
    color = "#17202A"
  ) +
  scale_fill_gradient(
    low = "#F2F6FA",
    high = "#2E6FA8",
    limits = c(0, 1),
    na.value = "#E5E5E5",
    name = paste0(
      "Spearman rho\n",
      "signed evidence"
    )
  ) +
  labs(
    title = "Sensitivity of inferential patterns",
    subtitle = "Agreement with the primary all-cell logit analysis",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

combined_figure <- (
  panel_A + panel_B
) / (
  panel_C + panel_D
) +
  plot_layout(
    widths = c(1.15, 1),
    heights = c(1.1, 1)
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
    "GSE144744_Tcells_V2_curated_",
    "differential_abundance.pdf"
  )
)

figure_png <- file.path(
  OUTPUT_FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_V2_curated_",
    "differential_abundance.png"
  )
)

ggsave(
  filename = figure_pdf,
  plot = combined_figure,
  width = 15,
  height = 11.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = figure_png,
  plot = combined_figure,
  width = 15,
  height = 11.5,
  units = "in",
  dpi = 320,
  bg = "white"
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

message("\nPrimary V2 significance summary:")
print(as.data.frame(significance_summary), row.names = FALSE)

message("\nPrimary family-level results:")
print(
  as.data.frame(
    primary_family_results %>%
      select(
        label,
        mean_HI2_percent,
        mean_MS2_percent,
        difference_percentage_points,
        proportion_ratio,
        P.Value,
        FDR
      )
  ),
  row.names = FALSE
)

message("\nExploratory state-level results:")
print(
  as.data.frame(
    primary_state_results %>%
      select(
        label,
        immune_state_name,
        mean_HI2_percent,
        mean_MS2_percent,
        difference_percentage_points,
        proportion_ratio,
        P.Value,
        FDR
      )
  ),
  row.names = FALSE
)

message("\nSensitivity summary:")
print(as.data.frame(sensitivity_summary), row.names = FALSE)

message("\nSaved V2 DA PDF: ", figure_pdf)
message("Saved V2 DA PNG: ", figure_png)
message(
  "V2 curated family/state differential-abundance ",
  "analysis completed."
)
