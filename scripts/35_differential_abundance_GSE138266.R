# ============================================================
# 35_differential_abundance_GSE138266.R
# Donor-aware external differential-abundance analysis
# ============================================================

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

PREDICTION_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_curated_state_predictions.rds"
)

OUTPUT_TABLE_DIR <- paste0(
  "results/tables/external_validation/GSE138266/",
  "differential_abundance"
)

OUTPUT_MODEL_DIR <- paste0(
  "results/models/external_validation/GSE138266"
)

OUTPUT_FIGURE_DIR <- paste0(
  "results/figures/external_validation/GSE138266"
)

dir.create(OUTPUT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PREDICTION_FILE)) {
  stop("Prediction file not found. Run scripts 33 and 34 first.")
}

message("Loading GSE138266 transferred-state predictions...")
predictions <- readRDS(PREDICTION_FILE)

required_columns <- c(
  "cell_id",
  "predicted.id",
  "prediction.score.max",
  "sample",
  "donor",
  "group",
  "compartment",
  "immune_state_name"
)

missing_columns <- setdiff(required_columns, colnames(predictions))

if (length(missing_columns) > 0L) {
  stop(
    "The prediction table is missing column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

all_state_levels <- sprintf("IVS-%02d", 1:11)
observed_state_levels <- setdiff(all_state_levels, "IVS-02")

higher_order_family_map <- c(
  "IVS-01" = "CD4 T-cell states",
  "IVS-02" = "CD4 T-cell states",
  "IVS-03" = "CD4 T-cell states",
  "IVS-04" = "CD4 T-cell states",
  "IVS-05" = "Cytotoxic/effector CD8",
  "IVS-06" = "Naive/memory CD8",
  "IVS-07" = "Cytotoxic/effector CD8",
  "IVS-08" = "Innate-like cytotoxic T",
  "IVS-09" = "Regulatory T",
  "IVS-10" = "Innate-like cytotoxic T",
  "IVS-11" = "Interferon-responsive T"
)

family_levels <- c(
  "CD4 T-cell states",
  "Cytotoxic/effector CD8",
  "Naive/memory CD8",
  "Innate-like cytotoxic T",
  "Regulatory T",
  "Interferon-responsive T"
)

state_lookup <- predictions %>%
  transmute(
    immune_state_id = as.character(predicted.id),
    immune_state_name = as.character(immune_state_name)
  ) %>%
  distinct() %>%
  arrange(immune_state_id)

cell_metadata <- predictions %>%
  transmute(
    cell_id = as.character(cell_id),
    sample = as.character(sample),
    donor = as.character(donor),
    group = as.character(group),
    compartment = as.character(compartment),
    immune_state_id = as.character(predicted.id),
    higher_order_family = unname(
      higher_order_family_map[as.character(predicted.id)]
    ),
    prediction_score = as.numeric(prediction.score.max)
  )

if (
  anyNA(cell_metadata$immune_state_id) ||
    anyNA(cell_metadata$higher_order_family) ||
    anyNA(cell_metadata$prediction_score)
) {
  stop("Missing state, family, or score values were found.")
}

sample_info <- cell_metadata %>%
  distinct(sample, donor, group, compartment)

if (
  nrow(sample_info) != 22L ||
    n_distinct(sample_info$donor) != 12L
) {
  stop("The external sample design failed structural validation.")
}

paired_donors <- sample_info %>%
  count(donor, name = "n_samples") %>%
  filter(n_samples == 2L) %>%
  pull(donor)

paired_sample_info <- sample_info %>%
  filter(donor %in% paired_donors) %>%
  group_by(donor) %>%
  filter(
    n() == 2L,
    n_distinct(compartment) == 2L,
    all(c("PBMC", "CSF") %in% compartment)
  ) %>%
  ungroup()

paired_donors <- sort(unique(paired_sample_info$donor))

if (
  length(paired_donors) != 10L ||
    nrow(paired_sample_info) != 20L
) {
  stop(
    "Expected 10 donors with paired PBMC and CSF samples; found ",
    length(paired_donors),
    "."
  )
}

write.csv(
  sample_info,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_sample_design.csv"),
  row.names = FALSE
)

write.csv(
  paired_sample_info,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_paired_compartment_design.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Analysis sets and comparison definitions
# ------------------------------------------------------------

analysis_sets <- list(
  all_cells = cell_metadata,
  score_050 = cell_metadata %>%
    filter(prediction_score >= 0.50)
)

analysis_set_summary <- bind_rows(lapply(
  names(analysis_sets),
  function(analysis_name) {
    analysis_sets[[analysis_name]] %>%
      summarise(
        analysis_set = analysis_name,
        n_cells = n(),
        n_samples = n_distinct(sample),
        n_donors = n_distinct(donor),
        n_states = n_distinct(immune_state_id),
        n_families = n_distinct(higher_order_family)
      )
  }
))

write.csv(
  analysis_set_summary,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_analysis_set_summary.csv"),
  row.names = FALSE
)

comparison_definitions <- list(
  PBMC_MS_vs_control = list(
    sample_filter = function(x) {
      x %>% filter(compartment == "PBMC")
    },
    reference_level = "Control",
    test_level = "MS",
    design_variable = "group",
    paired = FALSE,
    expected_samples = 10L,
    expected_donors = 10L,
    display_name = "PBMC: MS vs IIH control"
  ),
  CSF_MS_vs_control = list(
    sample_filter = function(x) {
      x %>% filter(compartment == "CSF")
    },
    reference_level = "Control",
    test_level = "MS",
    design_variable = "group",
    paired = FALSE,
    expected_samples = 12L,
    expected_donors = 12L,
    display_name = "CSF: MS vs IIH control"
  ),
  CSF_vs_PBMC_paired = list(
    sample_filter = function(x) {
      x %>% filter(donor %in% paired_donors)
    },
    reference_level = "PBMC",
    test_level = "CSF",
    design_variable = "compartment",
    paired = TRUE,
    expected_samples = 20L,
    expected_donors = 10L,
    display_name = "CSF vs paired PBMC"
  )
)

# ------------------------------------------------------------
# Propeller/limma fitting function
# ------------------------------------------------------------

run_external_propeller <- function(
    metadata,
    label_column,
    label_levels,
    resolution,
    analysis_set,
    transformation,
    comparison_name,
    comparison_definition) {

  cell_subset <- comparison_definition$sample_filter(metadata)

  design_metadata <- cell_subset %>%
    distinct(sample, donor, group, compartment)

  if (
    nrow(design_metadata) != comparison_definition$expected_samples ||
      n_distinct(design_metadata$donor) !=
        comparison_definition$expected_donors
  ) {
    stop(
      "Unexpected sample design for ", comparison_name,
      "/", resolution, "/", analysis_set, "."
    )
  }

  labels <- factor(cell_subset[[label_column]], levels = label_levels)

  if (anyNA(labels)) {
    stop("Unrecognized labels in ", comparison_name, "/", resolution, ".")
  }

  propeller_input <- speckle::getTransformedProps(
    clusters = labels,
    sample = cell_subset$sample,
    transform = transformation
  )

  sample_order <- colnames(propeller_input$TransformedProps)
  design_metadata <- design_metadata[
    match(sample_order, design_metadata$sample),
    ,
    drop = FALSE
  ]

  if (
    anyNA(design_metadata$sample) ||
      !identical(sample_order, design_metadata$sample)
  ) {
    stop("Unable to align Propeller samples for ", comparison_name, ".")
  }

  design_factor <- factor(
    design_metadata[[comparison_definition$design_variable]],
    levels = c(
      comparison_definition$reference_level,
      comparison_definition$test_level
    )
  )

  if (comparison_definition$paired) {
    donor_factor <- factor(design_metadata$donor)
    design <- model.matrix(~ 0 + design_factor + donor_factor)
  } else {
    design <- model.matrix(~ 0 + design_factor)
  }

  colnames(design) <- sub("^design_factor", "condition", colnames(design))
  colnames(design) <- sub("^donor_factor", "donor", colnames(design))
  colnames(design) <- make.names(colnames(design), unique = TRUE)

  if (qr(design)$rank != ncol(design)) {
    stop("The design matrix is not full rank for ", comparison_name, ".")
  }

  reference_column <- make.names(
    paste0("condition", comparison_definition$reference_level)
  )
  test_column <- make.names(
    paste0("condition", comparison_definition$test_level)
  )

  if (!all(c(reference_column, test_column) %in% colnames(design))) {
    stop("Unable to locate design columns for ", comparison_name, ".")
  }

  contrast_expression <- paste0(test_column, "-", reference_column)

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
  ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  propeller_result$label <- rownames(propeller_result)

  if (!all(c("P.Value", "FDR") %in% colnames(propeller_result))) {
    stop("Propeller did not return P.Value and FDR.")
  }

  proportion_matrix <- propeller_input$Proportions
  reference_columns <- design_factor == comparison_definition$reference_level
  test_columns <- design_factor == comparison_definition$test_level

  effects <- data.frame(
    label = rownames(proportion_matrix),
    mean_reference_percent = 100 * rowMeans(
      proportion_matrix[, reference_columns, drop = FALSE]
    ),
    mean_test_percent = 100 * rowMeans(
      proportion_matrix[, test_columns, drop = FALSE]
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

  results <- effects %>%
    left_join(propeller_result, by = "label") %>%
    mutate(
      resolution = resolution,
      analysis_set = analysis_set,
      transformation = transformation,
      contrast = comparison_name,
      contrast_display = comparison_definition$display_name,
      reference_level = comparison_definition$reference_level,
      test_level = comparison_definition$test_level,
      n_samples = nrow(design_metadata),
      n_donors = n_distinct(design_metadata$donor),
      paired_design = comparison_definition$paired,
      significance_class = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        P.Value < 0.05 ~ "Nominal P < 0.05",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(factor(label, levels = label_levels))

  if (nrow(results) != length(label_levels)) {
    stop(
      "Expected ", length(label_levels), " labels but obtained ",
      nrow(results), " for ", comparison_name, "/", resolution, "."
    )
  }

  list(
    results = results,
    design = design,
    design_metadata = design_metadata,
    contrast_matrix = contrast_matrix,
    transformed_proportions = propeller_input$TransformedProps
  )
}

# ------------------------------------------------------------
# Primary and sensitivity analyses
# ------------------------------------------------------------

analysis_configuration <- expand.grid(
  resolution = c("family", "state"),
  analysis_set = c("all_cells", "score_050"),
  transformation = c("logit", "asin"),
  contrast = names(comparison_definitions),
  stringsAsFactors = FALSE
) %>%
  arrange(contrast, resolution, analysis_set, transformation)

all_results <- list()
all_models <- list()

for (i in seq_len(nrow(analysis_configuration))) {
  configuration <- analysis_configuration[i, , drop = FALSE]

  resolution_name <- configuration$resolution[[1]]
  analysis_name <- configuration$analysis_set[[1]]
  transformation_name <- configuration$transformation[[1]]
  comparison_name <- configuration$contrast[[1]]

  if (resolution_name == "family") {
    label_column <- "higher_order_family"
    label_levels <- family_levels
  } else {
    label_column <- "immune_state_id"
    label_levels <- observed_state_levels
  }

  model_key <- paste(
    comparison_name,
    resolution_name,
    analysis_name,
    transformation_name,
    sep = "__"
  )

  message("Running ", model_key, "...")

  fitted <- run_external_propeller(
    metadata = analysis_sets[[analysis_name]],
    label_column = label_column,
    label_levels = label_levels,
    resolution = resolution_name,
    analysis_set = analysis_name,
    transformation = transformation_name,
    comparison_name = comparison_name,
    comparison_definition = comparison_definitions[[comparison_name]]
  )

  all_results[[model_key]] <- fitted$results
  all_models[[model_key]] <- fitted[
    c(
      "design",
      "design_metadata",
      "contrast_matrix",
      "transformed_proportions"
    )
  ]
}

combined_results <- bind_rows(all_results)

primary_results <- combined_results %>%
  filter(
    analysis_set == "all_cells",
    transformation == "logit"
  )

primary_state_results <- primary_results %>%
  filter(resolution == "state") %>%
  left_join(
    state_lookup,
    by = c("label" = "immune_state_id")
  )

primary_family_results <- primary_results %>%
  filter(resolution == "family")

# ------------------------------------------------------------
# Sensitivity summaries
# ------------------------------------------------------------

primary_effects <- primary_results %>%
  select(
    contrast,
    resolution,
    label,
    primary_effect = difference_percentage_points
  )

sensitivity_summary <- combined_results %>%
  filter(
    !(analysis_set == "all_cells" & transformation == "logit")
  ) %>%
  left_join(
    primary_effects,
    by = c("contrast", "resolution", "label")
  ) %>%
  group_by(contrast, contrast_display, resolution, analysis_set, transformation) %>%
  summarise(
    pearson_effect_correlation = cor(
      primary_effect,
      difference_percentage_points,
      method = "pearson",
      use = "complete.obs"
    ),
    spearman_effect_correlation = cor(
      primary_effect,
      difference_percentage_points,
      method = "spearman",
      use = "complete.obs"
    ),
    direction_agreement_percent = 100 * mean(
      sign(primary_effect) == sign(difference_percentage_points),
      na.rm = TRUE
    ),
    nominal_p_below_005 = sum(P.Value < 0.05, na.rm = TRUE),
    FDR_below_005 = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

significance_summary <- primary_results %>%
  group_by(contrast, contrast_display, resolution) %>%
  summarise(
    tested_labels = n(),
    nominal_p_below_005 = sum(P.Value < 0.05, na.rm = TRUE),
    FDR_below_005 = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  combined_results,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_all_results.csv"),
  row.names = FALSE
)

write.csv(
  primary_state_results,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_primary_state_results.csv"),
  row.names = FALSE
)

write.csv(
  primary_family_results,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_primary_family_results.csv"),
  row.names = FALSE
)

write.csv(
  sensitivity_summary,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_sensitivity_summary.csv"),
  row.names = FALSE
)

write.csv(
  significance_summary,
  file.path(OUTPUT_TABLE_DIR, "GSE138266_DA_significance_summary.csv"),
  row.names = FALSE
)

saveRDS(
  all_models,
  file.path(OUTPUT_MODEL_DIR, "GSE138266_DA_models.rds"),
  compress = "xz"
)

# ------------------------------------------------------------
# Publication-oriented summary figure
# ------------------------------------------------------------

contrast_levels <- c(
  "PBMC_MS_vs_control",
  "CSF_MS_vs_control",
  "CSF_vs_PBMC_paired"
)

contrast_labels <- c(
  "PBMC_MS_vs_control" = "PBMC: MS\nvs IIH control",
  "CSF_MS_vs_control" = "CSF: MS\nvs IIH control",
  "CSF_vs_PBMC_paired" = "CSF vs paired\nPBMC"
)

prepare_heatmap_data <- function(results, label_levels) {
  results %>%
    mutate(
      contrast = factor(contrast, levels = contrast_levels),
      label = factor(label, levels = rev(label_levels)),
      annotation = paste0(
        ifelse(difference_percentage_points > 0, "+", ""),
        sprintf("%.2f", difference_percentage_points),
        case_when(
          FDR < 0.05 ~ "**",
          P.Value < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    )
}

family_plot_data <- prepare_heatmap_data(
  primary_family_results,
  family_levels
)

state_plot_data <- prepare_heatmap_data(
  primary_state_results,
  observed_state_levels
)

effect_limit <- max(
  abs(c(
    family_plot_data$difference_percentage_points,
    state_plot_data$difference_percentage_points
  )),
  na.rm = TRUE
)

effect_limit <- ceiling(effect_limit)

heatmap_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9, face = "bold"),
    legend.position = "right"
  )

panel_A <- ggplot(
  family_plot_data,
  aes(x = contrast, y = label, fill = difference_percentage_points)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = annotation), size = 3.2) +
  scale_x_discrete(labels = contrast_labels) +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#E45756",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    name = "Difference\n(percentage points)"
  ) +
  labs(
    title = "External differential abundance at family level",
    subtitle = paste0(
      "Positive values indicate enrichment in MS or CSF; ",
      "* nominal P < 0.05, ** FDR < 0.05"
    )
  ) +
  heatmap_theme

panel_B <- ggplot(
  state_plot_data,
  aes(x = contrast, y = label, fill = difference_percentage_points)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = annotation), size = 3.0) +
  scale_x_discrete(labels = contrast_labels) +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#E45756",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    name = "Difference\n(percentage points)"
  ) +
  labs(
    title = "External differential abundance at immune-state level",
    subtitle = "IVS-02 was not recovered and was therefore not tested"
  ) +
  heatmap_theme

sensitivity_plot_data <- sensitivity_summary %>%
  mutate(
    contrast = factor(contrast, levels = contrast_levels),
    sensitivity = case_when(
      analysis_set == "all_cells" & transformation == "asin" ~
        "Arcsin transformation",
      analysis_set == "score_050" & transformation == "logit" ~
        "Transferred score >= 0.50",
      analysis_set == "score_050" & transformation == "asin" ~
        "Score >= 0.50 + arcsin",
      TRUE ~ paste(analysis_set, transformation)
    ),
    sensitivity = factor(
      sensitivity,
      levels = c(
        "Arcsin transformation",
        "Transferred score >= 0.50",
        "Score >= 0.50 + arcsin"
      )
    ),
    annotation = paste0(
      "r = ", sprintf("%.2f", pearson_effect_correlation),
      "\n", sprintf("%.0f%%", direction_agreement_percent),
      " direction"
    )
  )

panel_C <- ggplot(
  sensitivity_plot_data,
  aes(x = contrast, y = sensitivity, fill = pearson_effect_correlation)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = annotation), size = 3.0) +
  facet_wrap(~ resolution, ncol = 1, scales = "free_y") +
  scale_x_discrete(labels = contrast_labels) +
  scale_fill_gradient(
    low = "#E8EEF5",
    high = "#2F6FA7",
    limits = c(0, 1),
    name = "Pearson r"
  ) +
  labs(
    title = "Sensitivity of external abundance-effect estimates",
    subtitle = "Agreement with the primary all-cell logit analysis",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9, face = "bold"),
    strip.text = element_text(face = "bold")
  )

combined_figure <- (
  panel_A / panel_B / panel_C
) +
  plot_layout(heights = c(0.85, 1.2, 1.1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(
    OUTPUT_FIGURE_DIR,
    "GSE138266_external_differential_abundance.png"
  ),
  combined_figure,
  width = 13,
  height = 17,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    OUTPUT_FIGURE_DIR,
    "GSE138266_external_differential_abundance.pdf"
  ),
  combined_figure,
  width = 13,
  height = 17,
  device = grDevices::pdf,
  bg = "white"
)

message("")
message("GSE138266 external differential-abundance analysis completed.")
message("")
message("Primary significance summary:")
print(significance_summary, n = Inf)
message("")
message("Top primary results by nominal P value:")
top_primary_results <- primary_results %>%
  arrange(P.Value) %>%
  select(
    contrast,
    resolution,
    label,
    mean_reference_percent,
    mean_test_percent,
    difference_percentage_points,
    P.Value,
    FDR
  ) %>%
  head(20) %>%
  as.data.frame()

print(top_primary_results, row.names = FALSE)
