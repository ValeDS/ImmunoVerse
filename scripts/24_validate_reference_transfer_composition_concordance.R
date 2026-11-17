# ============================================================
# 24_validate_reference_transfer_composition_concordance.R
# Within-sample concordance between curated-reference and
# transferred immune-state compositions
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

MIN_REFERENCE_CELLS <- 100L
MIN_TRANSFERRED_CELLS <- 100L

mapped_object_file <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_curated_atlas_full_mapped.rds"
  )
)

qc_dir <- file.path(
  TAB_DIR,
  "curated_V3_mapping_QC"
)

figure_pdf <- file.path(
  FIG_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_reference_transfer_composition_concordance.pdf"
  )
)

figure_png <- file.path(
  FIG_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_reference_transfer_composition_concordance.png"
  )
)

if (!file.exists(mapped_object_file)) {
  stop(
    "Missing mapped V3 object: ",
    mapped_object_file
  )
}

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load metadata only
# ------------------------------------------------------------

message("Loading the full mapped V3 object...")

full_v3 <- readRDS(
  mapped_object_file
)

required_metadata <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "curated_immune_state_id",
  "state_assignment_origin",
  "atlas_inclusion_status"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(full_v3@meta.data)
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing mapped-object metadata: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

meta <- full_v3@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  transmute(
    cell_id,
    sample =
      as.character(sample),
    donor =
      as.character(donor),
    group =
      as.character(group),
    cohort =
      as.character(cohort),
    immune_state_id =
      as.character(curated_immune_state_id),
    assignment_origin =
      as.character(state_assignment_origin),
    inclusion_status =
      as.character(atlas_inclusion_status)
  )

rm(
  full_v3
)

gc(
  verbose = FALSE
)

meta <- meta %>%
  filter(
    inclusion_status ==
      "included",
    assignment_origin %in%
      c(
        "curated_reference",
        "label_transfer"
      )
  )

state_order <- sprintf(
  "IVS-%02d",
  seq_len(11)
)

if (
  !setequal(
    unique(meta$immune_state_id),
    state_order
  )
) {
  stop(
    "The included metadata do not contain the expected 11 states."
  )
}

# ------------------------------------------------------------
# Build within-sample state compositions by assignment origin
# ------------------------------------------------------------

sample_information <- meta %>%
  select(
    sample,
    donor,
    group,
    cohort
  ) %>%
  distinct()

if (anyDuplicated(sample_information$sample)) {
  stop(
    "A sample is associated with inconsistent biological metadata."
  )
}

composition_long <- meta %>%
  count(
    sample,
    donor,
    group,
    cohort,
    assignment_origin,
    immune_state_id,
    name = "state_cells"
  ) %>%
  complete(
    nesting(
      sample,
      donor,
      group,
      cohort
    ),
    assignment_origin = c(
      "curated_reference",
      "label_transfer"
    ),
    immune_state_id =
      state_order,
    fill = list(
      state_cells = 0L
    )
  ) %>%
  group_by(
    sample,
    assignment_origin
  ) %>%
  mutate(
    origin_total_cells =
      sum(state_cells),
    state_percent = if_else(
      origin_total_cells > 0,
      100 * state_cells /
        origin_total_cells,
      NA_real_
    )
  ) %>%
  ungroup()

composition_wide <- composition_long %>%
  select(
    sample,
    donor,
    group,
    cohort,
    assignment_origin,
    immune_state_id,
    state_cells,
    origin_total_cells,
    state_percent
  ) %>%
  pivot_wider(
    names_from =
      assignment_origin,
    values_from = c(
      state_cells,
      origin_total_cells,
      state_percent
    ),
    names_sep = "__"
  ) %>%
  rename(
    reference_state_cells =
      state_cells__curated_reference,
    transferred_state_cells =
      state_cells__label_transfer,
    reference_total_cells =
      origin_total_cells__curated_reference,
    transferred_total_cells =
      origin_total_cells__label_transfer,
    reference_percent =
      state_percent__curated_reference,
    transferred_percent =
      state_percent__label_transfer
  ) %>%
  mutate(
    difference_percentage_points =
      transferred_percent -
      reference_percent,
    absolute_difference_percentage_points =
      abs(
        difference_percentage_points
      ),
    eligible_for_concordance =
      reference_total_cells >=
        MIN_REFERENCE_CELLS &
        transferred_total_cells >=
        MIN_TRANSFERRED_CELLS
  )

write.csv(
  composition_wide,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_V3_reference_transfer_state_composition.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Sample-level concordance
# ------------------------------------------------------------

safe_cor <- function(
  x,
  y,
  method
) {

  valid <- is.finite(x) &
    is.finite(y)

  if (
    sum(valid) < 3L ||
      stats::sd(x[valid]) == 0 ||
      stats::sd(y[valid]) == 0
  ) {
    return(
      NA_real_
    )
  }

  stats::cor(
    x[valid],
    y[valid],
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
    reference_cells =
      first(reference_total_cells),
    transferred_cells =
      first(transferred_total_cells),
    eligible_for_concordance =
      first(eligible_for_concordance),
    pearson = safe_cor(
      reference_percent,
      transferred_percent,
      "pearson"
    ),
    spearman = safe_cor(
      reference_percent,
      transferred_percent,
      "spearman"
    ),
    mean_absolute_difference =
      mean(
        absolute_difference_percentage_points,
        na.rm = TRUE
      ),
    maximum_absolute_difference = if (
      all(
        is.na(
          absolute_difference_percentage_points
        )
      )
    ) {
      NA_real_
    } else {
      max(
        absolute_difference_percentage_points,
        na.rm = TRUE
      )
    },
    .groups = "drop"
  ) %>%
  mutate(
    pearson = if_else(
      eligible_for_concordance,
      pearson,
      NA_real_
    ),
    spearman = if_else(
      eligible_for_concordance,
      spearman,
      NA_real_
    ),
    mean_absolute_difference = if_else(
      eligible_for_concordance,
      mean_absolute_difference,
      NA_real_
    ),
    maximum_absolute_difference = if_else(
      eligible_for_concordance,
      maximum_absolute_difference,
      NA_real_
    )
  )

write.csv(
  sample_concordance,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_V3_reference_transfer_concordance_by_sample.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# State-level concordance across eligible samples
# ------------------------------------------------------------

eligible_composition <- composition_wide %>%
  filter(
    eligible_for_concordance
  )

state_concordance <- eligible_composition %>%
  group_by(
    immune_state_id
  ) %>%
  summarise(
    evaluated_samples = n(),
    pearson = safe_cor(
      reference_percent,
      transferred_percent,
      "pearson"
    ),
    spearman = safe_cor(
      reference_percent,
      transferred_percent,
      "spearman"
    ),
    mean_reference_percent =
      mean(
        reference_percent
      ),
    mean_transferred_percent =
      mean(
        transferred_percent
      ),
    mean_difference_percentage_points =
      mean(
        difference_percentage_points
      ),
    mean_absolute_difference =
      mean(
        absolute_difference_percentage_points
      ),
    .groups = "drop"
  )

write.csv(
  state_concordance,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_V3_reference_transfer_concordance_by_state.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Summary metrics
# ------------------------------------------------------------

eligible_sample_concordance <-
  sample_concordance %>%
  filter(
    eligible_for_concordance
  )

concordance_summary <- data.frame(
  metric = c(
    "total_samples",
    "eligible_samples",
    "minimum_reference_cells",
    "minimum_transferred_cells",
    "median_sample_pearson",
    "median_sample_spearman",
    "median_sample_MAD_percentage_points",
    "maximum_sample_MAD_percentage_points"
  ),
  value = c(
    nrow(sample_concordance),
    nrow(eligible_sample_concordance),
    MIN_REFERENCE_CELLS,
    MIN_TRANSFERRED_CELLS,
    median(
      eligible_sample_concordance$pearson,
      na.rm = TRUE
    ),
    median(
      eligible_sample_concordance$spearman,
      na.rm = TRUE
    ),
    median(
      eligible_sample_concordance$mean_absolute_difference,
      na.rm = TRUE
    ),
    max(
      eligible_sample_concordance$mean_absolute_difference,
      na.rm = TRUE
    )
  )
)

write.csv(
  concordance_summary,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_V3_reference_transfer_concordance_summary.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Plotting setup
# ------------------------------------------------------------

state_palette <- c(
  "IVS-01" = "#3B6FB6",
  "IVS-02" = "#7BA7D1",
  "IVS-03" = "#4E9F62",
  "IVS-04" = "#F28E2B",
  "IVS-05" = "#D84A4A",
  "IVS-06" = "#63B7AF",
  "IVS-07" = "#9C6BA3",
  "IVS-08" = "#E8B923",
  "IVS-09" = "#C778A6",
  "IVS-10" = "#8C564B",
  "IVS-11" = "#2F4858"
)

cohort_palette <- c(
  "NAT_HI" = "#F28E2B",
  "RRMS_HI" = "#4C78A8",
  "PPMS_HI" = "#59A14F"
)

theme_publication <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.subtitle = element_text(
      colour = "#5B5B5B",
      size = 10,
      margin = margin(
        b = 8
      )
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "#222222"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.tag = element_text(
      face = "bold",
      size = 15
    )
  )

# ------------------------------------------------------------
# Panel A: all eligible sample-state pairs
# ------------------------------------------------------------

p_scatter <- ggplot(
  eligible_composition,
  aes(
    x = reference_percent,
    y = transferred_percent,
    colour = immune_state_id
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "#555555",
    linewidth = 0.55
  ) +
  geom_point(
    alpha = 0.68,
    size = 1.8
  ) +
  scale_colour_manual(
    values = state_palette,
    breaks = state_order,
    name = "Immune state"
  ) +
  coord_equal() +
  labs(
    title = "Reference and transferred compositions",
    subtitle = paste0(
      nrow(eligible_sample_concordance),
      " samples with at least ",
      MIN_TRANSFERRED_CELLS,
      " transferred cells"
    ),
    x = "Reference composition (%)",
    y = "Transferred composition (%)"
  ) +
  theme_publication +
  theme(
    legend.position = "right"
  )

# ------------------------------------------------------------
# Panel B: sample-level concordance versus mapping depth
# ------------------------------------------------------------

p_depth <- ggplot(
  eligible_sample_concordance,
  aes(
    x = transferred_cells,
    y = pearson,
    colour = cohort
  )
) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    colour = "#777777",
    linewidth = 0.5
  ) +
  geom_point(
    size = 3,
    alpha = 0.82
  ) +
  scale_colour_manual(
    values = cohort_palette,
    name = "Cohort"
  ) +
  scale_x_continuous(
    labels = scales::comma
  ) +
  scale_y_continuous(
    limits = c(
      0.78,
      1.01
    ),
    breaks = c(
      0.80,
      0.85,
      0.90,
      0.95,
      1.00
    )
  ) +
  labs(
    title = "Concordance as a function of mapping depth",
    subtitle = "Dashed line marks Pearson r = 0.80",
    x = "Transferred cells per sample",
    y = "Pearson correlation"
  ) +
  theme_publication

# ------------------------------------------------------------
# Panel C: sample-level concordance distribution
# ------------------------------------------------------------

sample_correlation_long <-
  eligible_sample_concordance %>%
  select(
    sample,
    cohort,
    pearson,
    spearman
  ) %>%
  pivot_longer(
    cols = c(
      pearson,
      spearman
    ),
    names_to = "metric",
    values_to = "correlation"
  ) %>%
  mutate(
    metric = recode(
      metric,
      pearson = "Pearson",
      spearman = "Spearman"
    )
  )

p_correlation <- ggplot(
  sample_correlation_long,
  aes(
    x = metric,
    y = correlation,
    fill = metric
  )
) +
  geom_violin(
    width = 0.78,
    trim = TRUE,
    colour = "#444444",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "#222222",
    linewidth = 0.4
  ) +
  geom_jitter(
    aes(
      colour = cohort
    ),
    width = 0.08,
    height = 0,
    size = 1.7,
    alpha = 0.72
  ) +
  scale_fill_manual(
    values = c(
      "Pearson" = "#4C78A8",
      "Spearman" = "#9C6BA3"
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = cohort_palette,
    name = "Cohort"
  ) +
  scale_y_continuous(
    limits = c(
      0.78,
      1.01
    ),
    breaks = c(
      0.80,
      0.85,
      0.90,
      0.95,
      1.00
    )
  ) +
  labs(
    title = "Within-sample composition concordance",
    subtitle = "Each point represents one eligible biological sample",
    x = NULL,
    y = "Correlation"
  ) +
  theme_publication

# ------------------------------------------------------------
# Panel D: state-by-sample composition differences
# ------------------------------------------------------------

sample_order <- eligible_sample_concordance %>%
  arrange(
    cohort,
    group,
    sample
  ) %>%
  pull(
    sample
  )

heatmap_data <- eligible_composition %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_order
    ),
    immune_state_id = factor(
      immune_state_id,
      levels = rev(
        state_order
      )
    )
  )

difference_limit <- max(
  abs(
    heatmap_data$difference_percentage_points
  ),
  na.rm = TRUE
)

difference_limit <- max(
  5,
  ceiling(
    difference_limit
  )
)

p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = sample,
    y = immune_state_id,
    fill = difference_percentage_points
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.18
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#D84A4A",
    midpoint = 0,
    limits = c(
      -difference_limit,
      difference_limit
    ),
    oob = scales::squish,
    name = expression(
      Delta *
        " percentage points"
    )
  ) +
  labs(
    title = "State-specific composition differences",
    subtitle = "Transferred minus reference percentage within each sample",
    x = "Biological sample",
    y = NULL
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.subtitle = element_text(
      colour = "#5B5B5B",
      size = 10,
      margin = margin(
        b = 8
      )
    ),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      size = 6.5
    ),
    axis.text.y = element_text(
      colour = "#222222"
    ),
    panel.grid = element_blank(),
    legend.title = element_text(
      face = "bold"
    ),
    plot.tag = element_text(
      face = "bold",
      size = 15
    )
  )

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

concordance_figure <- (
  p_scatter |
    p_depth
) / (
  p_correlation |
    p_heatmap
) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    widths = c(
      1,
      1.15
    ),
    heights = c(
      1,
      1.05
    )
  )

ggsave(
  filename = figure_pdf,
  plot = concordance_figure,
  width = 16,
  height = 11.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = figure_png,
  plot = concordance_figure,
  width = 16,
  height = 11.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("\nReference-transfer concordance summary:")
print(
  concordance_summary,
  row.names = FALSE
)

message("\nLowest eligible sample-level Pearson correlations:")
print(
  as.data.frame(
    eligible_sample_concordance %>%
      arrange(
        pearson
      ) %>%
      select(
        sample,
        group,
        cohort,
        reference_cells,
        transferred_cells,
        pearson,
        spearman,
        mean_absolute_difference
      ) %>%
      head(10)
  ),
  row.names = FALSE
)

message(
  "\nSaved concordance PDF: ",
  figure_pdf
)

message(
  "Saved concordance PNG: ",
  figure_png
)

message(
  "Reference-transfer composition validation completed."
)
