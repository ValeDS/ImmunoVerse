# ============================================================
# 27_validate_V3_to_V2_cross_chemistry_transfer.R
# QC and within-sample composition concordance for the
# cross-chemistry transfer of curated V3 states to V2 cells
# ============================================================

suppressPackageStartupMessages({
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

V3_METADATA_FILE <- paste0(
  "data/processed/",
  "GSE144744_Tcells_V3_curated_DA_metadata.rds"
)

ATLAS_TABLE_FILE <- paste0(
  "results/tables/",
  "GSE144744_Tcells_V3_Harmony_reference_",
  "curated_immune_state_atlas.csv"
)

OUTPUT_TABLE_DIR <- paste0(
  "results/tables/",
  "V3_to_V2_cross_chemistry_validation"
)

OUTPUT_FIGURE_DIR <- "results/figures"

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

required_files <- c(
  V2_PREDICTION_FILE,
  V3_METADATA_FILE,
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
# Load compact inputs
# ------------------------------------------------------------

message("Loading V2 cross-chemistry predictions...")
v2_predictions <- readRDS(V2_PREDICTION_FILE)

message("Loading curated V3 cell metadata...")
v3_metadata <- readRDS(V3_METADATA_FILE)

atlas_table <- read.csv(
  ATLAS_TABLE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_v2_columns <- c(
  "cell_id",
  "predicted.id",
  "prediction.score.max",
  "sample",
  "donor",
  "group",
  "cohort",
  "V_10X"
)

missing_v2_columns <- setdiff(
  required_v2_columns,
  colnames(v2_predictions)
)

if (length(missing_v2_columns) > 0L) {
  stop(
    "The V2 prediction table is missing column(s): ",
    paste(missing_v2_columns, collapse = ", ")
  )
}

required_v3_columns <- c(
  "cell_id",
  "sample",
  "donor",
  "group",
  "cohort",
  "curated_immune_state_id",
  "atlas_inclusion_status"
)

missing_v3_columns <- setdiff(
  required_v3_columns,
  colnames(v3_metadata)
)

if (length(missing_v3_columns) > 0L) {
  stop(
    "The V3 metadata table is missing column(s): ",
    paste(missing_v3_columns, collapse = ", ")
  )
}

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
    state_family = as.character(state_family)
  ) %>%
  distinct()

state_levels <- atlas_definitions$immune_state_id

if (
  length(state_levels) != 11L ||
    anyDuplicated(state_levels)
) {
  stop("Expected 11 unique curated immune-state definitions.")
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

higher_family_levels <- c(
  "CD4 T-cell states",
  "Cytotoxic/effector CD8",
  "Naive/memory CD8",
  "Innate-like cytotoxic T",
  "Regulatory T",
  "Interferon-responsive T"
)

v2_predictions <- v2_predictions %>%
  mutate(
    predicted.id = as.character(predicted.id),
    sample = as.character(sample),
    donor = as.character(donor),
    group = as.character(group),
    cohort = as.character(cohort),
    V_10X = as.character(V_10X),
    prediction.score.max = as.numeric(
      prediction.score.max
    )
  )

v3_included <- v3_metadata %>%
  filter(
    atlas_inclusion_status == "included",
    !is.na(curated_immune_state_id)
  ) %>%
  mutate(
    curated_immune_state_id = as.character(
      curated_immune_state_id
    ),
    sample = as.character(sample),
    donor = as.character(donor),
    group = as.character(group),
    cohort = as.character(cohort)
  )

# ------------------------------------------------------------
# Global and state-specific V2 transfer QC
# ------------------------------------------------------------

global_score_summary <- data.frame(
  predicted_cells = nrow(v2_predictions),
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

state_score_summary <- v2_predictions %>%
  group_by(predicted.id) %>%
  summarise(
    predicted_cells = n(),
    sample_prevalence = n_distinct(sample),
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
      0.25,
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
  rename(
    immune_state_id = predicted.id
  ) %>%
  left_join(
    atlas_definitions,
    by = "immune_state_id"
  ) %>%
  mutate(
    immune_state_id = factor(
      immune_state_id,
      levels = state_levels
    )
  ) %>%
  arrange(immune_state_id) %>%
  mutate(
    immune_state_id = as.character(
      immune_state_id
    )
  )

sample_score_summary <- v2_predictions %>%
  group_by(
    sample,
    donor,
    group,
    cohort
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
    pct_below_050 = 100 * mean(
      prediction.score.max < 0.50,
      na.rm = TRUE
    ),
    pct_below_075 = 100 * mean(
      prediction.score.max < 0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(sample)

# ------------------------------------------------------------
# Identify samples measured with both V2 and V3 chemistry
# ------------------------------------------------------------

v2_sample_metadata <- v2_predictions %>%
  distinct(
    sample,
    donor,
    group,
    cohort
  )

v3_sample_metadata <- v3_included %>%
  distinct(
    sample,
    donor,
    group,
    cohort
  )

overlap_metadata <- inner_join(
  v2_sample_metadata,
  v3_sample_metadata,
  by = c(
    "sample",
    "donor",
    "group",
    "cohort"
  )
) %>%
  arrange(sample)

expected_overlap_samples <- c(
  "HH-OX-15",
  "HH-OX-16",
  "HH-OX-19",
  "HH-OX-20"
)

if (
  nrow(overlap_metadata) != 4L ||
    !setequal(
      overlap_metadata$sample,
      expected_overlap_samples
    )
) {
  stop(
    "Expected four V2/V3-overlap samples (HH-OX-15, ",
    "HH-OX-16, HH-OX-19, HH-OX-20), but found: ",
    paste(overlap_metadata$sample, collapse = ", ")
  )
}

overlap_samples <- overlap_metadata$sample

# ------------------------------------------------------------
# Calculate V2 and V3 state compositions in overlap samples
# ------------------------------------------------------------

v2_overlap_counts <- v2_predictions %>%
  filter(sample %in% overlap_samples) %>%
  count(
    sample,
    predicted.id,
    name = "n_cells"
  ) %>%
  transmute(
    sample,
    chemistry = "V2",
    immune_state_id = as.character(predicted.id),
    n_cells
  )

v3_overlap_counts <- v3_included %>%
  filter(sample %in% overlap_samples) %>%
  count(
    sample,
    curated_immune_state_id,
    name = "n_cells"
  ) %>%
  transmute(
    sample,
    chemistry = "V3",
    immune_state_id = as.character(
      curated_immune_state_id
    ),
    n_cells
  )

composition_grid <- expand_grid(
  sample = overlap_samples,
  chemistry = c("V2", "V3"),
  immune_state_id = state_levels
)

overlap_composition <- composition_grid %>%
  left_join(
    bind_rows(
      v2_overlap_counts,
      v3_overlap_counts
    ),
    by = c(
      "sample",
      "chemistry",
      "immune_state_id"
    )
  ) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L)
  ) %>%
  group_by(sample, chemistry) %>%
  mutate(
    total_cells = sum(n_cells),
    percent = 100 * n_cells / total_cells
  ) %>%
  ungroup() %>%
  left_join(
    overlap_metadata,
    by = "sample"
  ) %>%
  left_join(
    atlas_definitions,
    by = "immune_state_id"
  )

composition_wide <- overlap_composition %>%
  select(
    sample,
    donor,
    group,
    cohort,
    immune_state_id,
    chemistry,
    percent
  ) %>%
  pivot_wider(
    names_from = chemistry,
    values_from = percent
  ) %>%
  mutate(
    difference_V2_minus_V3 =
      V2 - V3,
    absolute_difference =
      abs(difference_V2_minus_V3)
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

sample_concordance <- composition_wide %>%
  group_by(
    sample,
    donor,
    group,
    cohort
  ) %>%
  summarise(
    pearson = safe_correlation(
      V2,
      V3,
      method = "pearson"
    ),
    spearman = safe_correlation(
      V2,
      V3,
      method = "spearman"
    ),
    mean_absolute_difference =
      mean(
        absolute_difference,
        na.rm = TRUE
      ),
    maximum_absolute_difference =
      max(
        absolute_difference,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(sample)

overall_composition_concordance <- data.frame(
  overlap_samples = n_distinct(
    composition_wide$sample
  ),
  state_sample_pairs = nrow(
    composition_wide
  ),
  pooled_pearson = safe_correlation(
    composition_wide$V2,
    composition_wide$V3,
    method = "pearson"
  ),
  pooled_spearman = safe_correlation(
    composition_wide$V2,
    composition_wide$V3,
    method = "spearman"
  ),
  median_sample_pearson = median(
    sample_concordance$pearson,
    na.rm = TRUE
  ),
  median_sample_spearman = median(
    sample_concordance$spearman,
    na.rm = TRUE
  ),
  median_sample_MAD_percentage_points =
    median(
      sample_concordance$mean_absolute_difference,
      na.rm = TRUE
    ),
  maximum_sample_MAD_percentage_points =
    max(
      sample_concordance$mean_absolute_difference,
      na.rm = TRUE
    )
)

# ------------------------------------------------------------
# Higher-order family concordance
# ------------------------------------------------------------

family_composition_wide <- composition_wide %>%
  left_join(
    higher_family_definitions,
    by = "immune_state_id"
  )

if (anyNA(
  family_composition_wide$higher_order_family
)) {
  stop(
    "One or more immune states could not be assigned to a ",
    "higher-order family."
  )
}

family_composition_wide <-
  family_composition_wide %>%
  group_by(
    sample,
    donor,
    group,
    cohort,
    higher_order_family
  ) %>%
  summarise(
    V2 = sum(V2),
    V3 = sum(V3),
    .groups = "drop"
  ) %>%
  mutate(
    difference_V2_minus_V3 = V2 - V3,
    absolute_difference = abs(
      difference_V2_minus_V3
    )
  )

family_sample_concordance <-
  family_composition_wide %>%
  group_by(
    sample,
    donor,
    group,
    cohort
  ) %>%
  summarise(
    pearson = safe_correlation(
      V2,
      V3,
      method = "pearson"
    ),
    spearman = safe_correlation(
      V2,
      V3,
      method = "spearman"
    ),
    mean_absolute_difference =
      mean(
        absolute_difference,
        na.rm = TRUE
      ),
    maximum_absolute_difference =
      max(
        absolute_difference,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(sample)

family_overall_concordance <- data.frame(
  overlap_samples = n_distinct(
    family_composition_wide$sample
  ),
  family_sample_pairs = nrow(
    family_composition_wide
  ),
  pooled_pearson = safe_correlation(
    family_composition_wide$V2,
    family_composition_wide$V3,
    method = "pearson"
  ),
  pooled_spearman = safe_correlation(
    family_composition_wide$V2,
    family_composition_wide$V3,
    method = "spearman"
  ),
  median_sample_pearson = median(
    family_sample_concordance$pearson,
    na.rm = TRUE
  ),
  median_sample_spearman = median(
    family_sample_concordance$spearman,
    na.rm = TRUE
  ),
  median_sample_MAD_percentage_points =
    median(
      family_sample_concordance$
        mean_absolute_difference,
      na.rm = TRUE
    ),
  maximum_sample_MAD_percentage_points =
    max(
      family_sample_concordance$
        mean_absolute_difference,
      na.rm = TRUE
    )
)

# ------------------------------------------------------------
# Save validation tables
# ------------------------------------------------------------

write.csv(
  global_score_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V3_to_V2_global_score_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  state_score_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V3_to_V2_state_score_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_score_summary,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V3_to_V2_sample_score_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  overlap_composition,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_state_composition_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  composition_wide,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_state_composition_wide.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_concordance,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_sample_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_composition_concordance,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_overall_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  higher_family_definitions,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_curated_state_higher_order_family_map.csv"
  ),
  row.names = FALSE
)

write.csv(
  family_composition_wide,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_family_composition.csv"
  ),
  row.names = FALSE
)

write.csv(
  family_sample_concordance,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_family_sample_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  family_overall_concordance,
  file.path(
    OUTPUT_TABLE_DIR,
    "GSE144744_V2_V3_overlap_family_overall_concordance.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Publication-oriented validation figure
# ------------------------------------------------------------

state_palette <- c(
  "IVS-01" = "#3F77B5",
  "IVS-02" = "#78A6D0",
  "IVS-03" = "#4E9E60",
  "IVS-04" = "#F28E2B",
  "IVS-05" = "#E15759",
  "IVS-06" = "#76B7B2",
  "IVS-07" = "#9C6BAD",
  "IVS-08" = "#E7B416",
  "IVS-09" = "#C875A0",
  "IVS-10" = "#8C564B",
  "IVS-11" = "#2F4858"
)

plot_v2_predictions <- v2_predictions %>%
  mutate(
    plot_prediction_score = pmin(
      1,
      pmax(
        0,
        prediction.score.max
      )
    )
  )

panel_A <- ggplot(
  plot_v2_predictions,
  aes(x = plot_prediction_score)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 45,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.25
  ) +
  geom_density(
    color = "#173B57",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(0.50, 0.75, 0.90),
    linetype = "dashed",
    color = c("#C73E3A", "#E69500", "#2D8B45"),
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = c(0.50, 0.75, 0.90),
    y = Inf,
    label = c(
      sprintf(
        "%.1f%% < 0.50",
        global_score_summary$pct_below_050
      ),
      sprintf(
        "%.1f%% < 0.75",
        global_score_summary$pct_below_075
      ),
      sprintf(
        "%.1f%% < 0.90",
        global_score_summary$pct_below_090
      )
    ),
    hjust = c(1.05, 1.05, 1.05),
    vjust = 1.4,
    angle = 90,
    size = 3.2,
    color = c("#A8322E", "#B87400", "#247238")
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross-chemistry prediction confidence",
    subtitle = paste0(
      format(
        nrow(v2_predictions),
        big.mark = ","
      ),
      " V2 cells mapped from the curated V3 reference"
    ),
    x = "Prediction score",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

plot_v2_scores <- plot_v2_predictions %>%
  mutate(
    predicted.id = factor(
      predicted.id,
      levels = state_levels
    )
  )

panel_B <- ggplot(
  plot_v2_scores,
  aes(
    x = predicted.id,
    y = plot_prediction_score,
    fill = predicted.id
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "#555555",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "#444444",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 0.50,
    linetype = "dashed",
    color = "#C73E3A",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = state_palette,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Confidence varies across transferred states",
    subtitle = "Dashed line indicates a score of 0.50",
    x = NULL,
    y = "Prediction score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

family_palette <- c(
  "CD4 T-cell states" = "#4C78A8",
  "Cytotoxic/effector CD8" = "#E15759",
  "Naive/memory CD8" = "#76B7B2",
  "Innate-like cytotoxic T" = "#E7B416",
  "Regulatory T" = "#B279A2",
  "Interferon-responsive T" = "#2F4858"
)

plot_composition <- family_composition_wide %>%
  mutate(
    sample_label = paste0(
      sample,
      " (",
      group,
      ")"
    ),
    higher_order_family = factor(
      higher_order_family,
      levels = higher_family_levels
    )
  )

composition_limit <- max(
  c(
    plot_composition$V2,
    plot_composition$V3
  ),
  na.rm = TRUE
)

composition_limit <- ceiling(
  composition_limit / 5
) * 5

panel_C <- ggplot(
  plot_composition,
  aes(
    x = V3,
    y = V2,
    color = higher_order_family
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "#777777",
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.7,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ sample_label,
    ncol = 2
  ) +
  scale_color_manual(
    values = family_palette,
    drop = FALSE,
    name = "State family"
  ) +
  scale_x_continuous(
    limits = c(0, composition_limit),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0, composition_limit),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_fixed() +
  labs(
    title = "Within-sample V2-V3 family concordance",
    subtitle = paste0(
      "Each point represents one higher-order family; ",
      "the dashed line indicates identity"
    ),
    x = "V3 family abundance (%)",
    y = "V2 family abundance (%)"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(
      fill = "#F2F2F2",
      color = NA
    ),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      color = "#606060"
    )
  )

concordance_long <- family_sample_concordance %>%
  select(
    sample,
    group,
    pearson,
    spearman
  ) %>%
  pivot_longer(
    cols = c(pearson, spearman),
    names_to = "metric",
    values_to = "correlation"
  ) %>%
  mutate(
    sample_label = paste0(
      sample,
      "\n",
      group
    ),
    metric = recode(
      metric,
      pearson = "Pearson",
      spearman = "Spearman"
    )
  )

panel_D <- ggplot(
  concordance_long,
  aes(
    x = sample_label,
    y = correlation,
    color = metric,
    group = metric
  )
) +
  geom_hline(
    yintercept = 0.75,
    linetype = "dashed",
    color = "#E69500",
    linewidth = 0.7
  ) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      Pearson = "#2E6FA8",
      Spearman = "#C65A3A"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Family agreement in dual-chemistry samples",
    subtitle = paste0(
      "Median Pearson = ",
      sprintf(
        "%.3f",
        family_overall_concordance$
          median_sample_pearson
      ),
      "; median Spearman = ",
      sprintf(
        "%.3f",
        family_overall_concordance$
          median_sample_spearman
      )
    ),
    x = NULL,
    y = "Correlation"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
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
    widths = c(1.12, 1),
    heights = c(0.9, 1.25),
    guides = "collect"
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
    "GSE144744_Tcells_curated_atlas_",
    "V3_to_V2_cross_chemistry_validation.pdf"
  )
)

figure_png <- file.path(
  OUTPUT_FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_curated_atlas_",
    "V3_to_V2_cross_chemistry_validation.png"
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

message("\nGlobal V3-to-V2 transfer score summary:")
print(as.data.frame(global_score_summary), row.names = FALSE)

message("\nDual-chemistry sample concordance:")
print(as.data.frame(sample_concordance), row.names = FALSE)

message("\nOverall composition concordance:")
print(
  as.data.frame(overall_composition_concordance),
  row.names = FALSE
)

message("\nHigher-order family concordance:")
print(
  as.data.frame(family_sample_concordance),
  row.names = FALSE
)

message("\nOverall higher-order family concordance:")
print(
  as.data.frame(family_overall_concordance),
  row.names = FALSE
)

message("\nSaved validation PDF: ", figure_pdf)
message("Saved validation PNG: ", figure_png)
message(
  "Cross-chemistry transfer validation completed."
)
