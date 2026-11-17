# ============================================================
# 23_visualize_curated_V3_transfer_QC.R
# Publication-ready QC of curated V3 label transfer
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_curated_atlas"
)

assignment_file <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    "_full_assignments.rds"
  )
)

qc_dir <- file.path(
  TAB_DIR,
  "curated_V3_mapping_QC"
)

sample_qc_file <- file.path(
  qc_dir,
  paste0(
    GEO_ACCESSION,
    "_curated_V3_transfer_QC_by_sample.csv"
  )
)

state_qc_file <- file.path(
  qc_dir,
  paste0(
    GEO_ACCESSION,
    "_curated_V3_transfer_QC_by_state.csv"
  )
)

figure_pdf <- file.path(
  FIG_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_curated_atlas_transfer_QC.pdf"
  )
)

figure_png <- file.path(
  FIG_DIR,
  paste0(
    GEO_ACCESSION,
    "_Tcells_V3_curated_atlas_transfer_QC.png"
  )
)

required_files <- c(
  assignment_file,
  sample_qc_file,
  state_qc_file
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

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load validated outputs
# ------------------------------------------------------------

assignments <- readRDS(
  assignment_file
)

sample_qc <- read.csv(
  sample_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

state_qc <- read.csv(
  state_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_assignment_columns <- c(
  "curated_immune_state_id",
  "transfer_prediction_score",
  "state_assignment_origin",
  "atlas_inclusion_status",
  "cell_received_tie_repair"
)

if (
  !all(
    required_assignment_columns %in%
      colnames(assignments)
  )
) {
  stop(
    "The assignment table lacks required QC columns."
  )
}

transferred <- assignments %>%
  filter(
    state_assignment_origin ==
      "label_transfer",
    atlas_inclusion_status ==
      "included"
  ) %>%
  transmute(
    immune_state_id =
      as.character(curated_immune_state_id),
    prediction_score =
      as.numeric(transfer_prediction_score),
    cell_received_tie_repair =
      as.logical(cell_received_tie_repair)
  )

if (
  nrow(transferred) == 0L ||
    anyNA(transferred$prediction_score) ||
    any(
      !is.finite(
        transferred$prediction_score
      )
    )
) {
  stop(
    "No valid included transferred predictions were found."
  )
}

state_order <- sprintf(
  "IVS-%02d",
  seq_len(11)
)

if (
  !setequal(
    unique(transferred$immune_state_id),
    state_order
  )
) {
  stop(
    "The transferred assignments do not contain the expected 11 states."
  )
}

transferred$immune_state_id <- factor(
  transferred$immune_state_id,
  levels = state_order
)

state_qc$curated_immune_state_id <- factor(
  state_qc$curated_immune_state_id,
  levels = state_order
)

# ------------------------------------------------------------
# Palettes and theme
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

group_order <- c(
  "HI1",
  "MS1",
  "MS1_nat",
  "HI2",
  "MS2",
  "HI3",
  "PPMS"
)

group_palette <- c(
  "HI1" = "#7BA7D1",
  "MS1" = "#F28E2B",
  "MS1_nat" = "#9C6BA3",
  "HI2" = "#76B7B2",
  "MS2" = "#D84A4A",
  "HI3" = "#59A14F",
  "PPMS" = "#8C564B"
)

sample_qc$group <- factor(
  sample_qc$group,
  levels = group_order
)

sample_qc$cohort <- factor(
  sample_qc$cohort,
  levels = c(
    "NAT_HI",
    "RRMS_HI",
    "PPMS_HI"
  )
)

theme_publication <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      margin = margin(
        b = 3
      )
    ),
    plot.subtitle = element_text(
      colour = "#5B5B5B",
      size = 10,
      margin = margin(
        b = 9
      )
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "#222222"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
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
# Global QC summaries
# ------------------------------------------------------------

global_qc <- data.frame(
  metric = c(
    "included_transferred_cells",
    "median_score",
    "mean_score",
    "percent_below_050",
    "percent_below_075",
    "percent_below_090",
    "tie_repaired_cells"
  ),
  value = c(
    nrow(transferred),
    median(
      transferred$prediction_score
    ),
    mean(
      transferred$prediction_score
    ),
    100 * mean(
      transferred$prediction_score < 0.50
    ),
    100 * mean(
      transferred$prediction_score < 0.75
    ),
    100 * mean(
      transferred$prediction_score < 0.90
    ),
    sum(
      transferred$cell_received_tie_repair,
      na.rm = TRUE
    )
  )
)

write.csv(
  global_qc,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_transfer_global_QC.csv"
    )
  ),
  row.names = FALSE
)

samples_without_transfer <- sample_qc %>%
  filter(
    transferred_cells == 0 |
      is.na(median_transfer_score)
  )

write.csv(
  samples_without_transfer,
  file.path(
    qc_dir,
    paste0(
      GEO_ACCESSION,
      "_curated_V3_samples_without_included_transfer_cells.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Panel A: global score distribution
# ------------------------------------------------------------

score_annotation <- data.frame(
  threshold = c(
    0.50,
    0.75,
    0.90
  ),
  label = c(
    paste0(
      round(
        100 * mean(
          transferred$prediction_score < 0.50
        ),
        1
      ),
      "% < 0.50"
    ),
    paste0(
      round(
        100 * mean(
          transferred$prediction_score < 0.75
        ),
        1
      ),
      "% < 0.75"
    ),
    paste0(
      round(
        100 * mean(
          transferred$prediction_score < 0.90
        ),
        1
      ),
      "% < 0.90"
    )
  )
)

p_global <- ggplot(
  transferred,
  aes(
    x = prediction_score
  )
) +
  geom_histogram(
    aes(
      y = after_stat(
        density
      )
    ),
    bins = 60,
    boundary = 0,
    fill = "#4C78A8",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_density(
    colour = "#173F5F",
    linewidth = 0.8,
    adjust = 1.1
  ) +
  geom_vline(
    data = score_annotation,
    aes(
      xintercept = threshold
    ),
    linetype = "dashed",
    colour = c(
      "#B22222",
      "#D98C00",
      "#2E7D32"
    ),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = c(
      0.48,
      0.73,
      0.88
    ),
    y = Inf,
    label = score_annotation$label,
    hjust = 1,
    vjust = c(
      1.5,
      3.0,
      4.5
    ),
    size = 3.1,
    colour = c(
      "#8B1A1A",
      "#A66300",
      "#1B5E20"
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  labs(
    title = "Global transfer-confidence distribution",
    subtitle = paste0(
      format(
        nrow(transferred),
        big.mark = ","
      ),
      " included transferred cells"
    ),
    x = "Prediction score",
    y = "Density"
  ) +
  theme_publication

# ------------------------------------------------------------
# Panel B: score distribution by curated state
# ------------------------------------------------------------

p_state <- ggplot(
  transferred,
  aes(
    x = immune_state_id,
    y = prediction_score,
    fill = immune_state_id
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.25,
    colour = "#333333"
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "#222222",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 0.75,
    linetype = "dashed",
    colour = "#D98C00",
    linewidth = 0.55
  ) +
  scale_fill_manual(
    values = state_palette,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    )
  ) +
  labs(
    title = "Confidence varies across immune states",
    subtitle = "Dashed line indicates a descriptive score of 0.75",
    x = NULL,
    y = "Prediction score"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# ------------------------------------------------------------
# Panel C: sample-level confidence by clinical group
# ------------------------------------------------------------

sample_qc_plot <- sample_qc %>%
  filter(
    transferred_cells > 0,
    !is.na(median_transfer_score)
  )

p_sample <- ggplot(
  sample_qc_plot,
  aes(
    x = group,
    y = median_transfer_score,
    colour = group
  )
) +
  geom_boxplot(
    aes(
      group = group
    ),
    width = 0.58,
    outlier.shape = NA,
    colour = "#555555",
    fill = "white",
    linewidth = 0.45
  ) +
  geom_point(
    aes(
      size = transferred_cells
    ),
    alpha = 0.82,
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 123
    )
  ) +
  geom_hline(
    yintercept = 0.75,
    linetype = "dashed",
    colour = "#D98C00",
    linewidth = 0.55
  ) +
  facet_grid(
    cols = vars(
      cohort
    ),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_colour_manual(
    values = group_palette,
    drop = FALSE,
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(
      1.8,
      5.5
    ),
    name = "Transferred cells"
  ) +
  scale_y_continuous(
    limits = c(
      0.45,
      1
    ),
    breaks = seq(
      0.5,
      1,
      0.1
    )
  ) +
  labs(
    title = "Sample-level mapping confidence",
    subtitle = "Each point represents one biological sample",
    x = NULL,
    y = "Median prediction score"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# Panel D: state-level threshold sensitivity
# ------------------------------------------------------------

threshold_qc_long <- state_qc %>%
  select(
    curated_immune_state_id,
    pct_below_050,
    pct_below_075,
    pct_below_090
  ) %>%
  pivot_longer(
    cols = starts_with(
      "pct_below_"
    ),
    names_to = "threshold",
    values_to = "percent_cells"
  ) %>%
  mutate(
    threshold = recode(
      threshold,
      pct_below_050 = "Score < 0.50",
      pct_below_075 = "Score < 0.75",
      pct_below_090 = "Score < 0.90"
    ),
    threshold = factor(
      threshold,
      levels = c(
        "Score < 0.50",
        "Score < 0.75",
        "Score < 0.90"
      )
    ),
    curated_immune_state_id = factor(
      curated_immune_state_id,
      levels = rev(
        state_order
      )
    )
  )

p_threshold <- ggplot(
  threshold_qc_long,
  aes(
    x = threshold,
    y = curated_immune_state_id,
    fill = percent_cells
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percent_cells
      )
    ),
    size = 3.1,
    colour = "#111111"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#F7FBFF",
      "#9ECAE1",
      "#4292C6",
      "#084594"
    ),
    limits = c(
      0,
      100
    ),
    name = "% cells"
  ) +
  labs(
    title = "State-level confidence sensitivity",
    subtitle = "Percentage of transferred cells below each score threshold",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
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
        b = 9
      )
    ),
    axis.text = element_text(
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
# Assemble and save the publication figure
# ------------------------------------------------------------

qc_figure <- (
  p_global |
    p_state
) / (
  p_sample |
    p_threshold
) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    widths = c(
      1,
      1.12
    ),
    heights = c(
      1,
      1.05
    )
  )

ggsave(
  filename = figure_pdf,
  plot = qc_figure,
  width = 15.5,
  height = 11,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  filename = figure_png,
  plot = qc_figure,
  width = 15.5,
  height = 11,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("\nCurated V3 transfer QC summary:")
print(
  global_qc,
  row.names = FALSE
)

message(
  "\nSamples without included transferred cells: ",
  nrow(samples_without_transfer)
)

if (nrow(samples_without_transfer) > 0L) {
  print(
    as.data.frame(
      samples_without_transfer %>%
        select(
          sample,
          donor,
          group,
          cohort,
          curated_reference_cells,
          curation_excluded_cells
        )
    ),
    row.names = FALSE
  )
}

message(
  "\nSaved transfer-QC PDF: ",
  figure_pdf
)

message(
  "Saved transfer-QC PNG: ",
  figure_png
)

message(
  "Curated V3 transfer-QC visualization completed."
)
