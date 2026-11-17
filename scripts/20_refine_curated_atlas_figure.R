# ============================================================
# 20_refine_curated_atlas_figure.R
# Publication-oriented visualization of the curated V3 atlas
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load curated atlas
# ------------------------------------------------------------

curated_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas.rds"
  )
)

if (!file.exists(curated_rds)) {
  stop(
    "Curated atlas not found: ",
    curated_rds
  )
}

message("Loading curated V3 atlas...")

tcells_v3_curated <- readRDS(
  curated_rds
)

DefaultAssay(tcells_v3_curated) <- "RNA"

if (
  !"immune_state_id" %in%
    colnames(tcells_v3_curated@meta.data)
) {
  stop(
    "Metadata column immune_state_id is missing."
  )
}

if (
  !"umap_harmony" %in%
    Reductions(tcells_v3_curated)
) {
  stop(
    "Reduction umap_harmony is missing."
  )
}

Idents(tcells_v3_curated) <-
  "immune_state_id"

# ------------------------------------------------------------
# Define higher-order state families
# ------------------------------------------------------------

major_family_map <- c(
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

major_family_levels <- c(
  "CD4 T-cell states",
  "Cytotoxic/effector CD8",
  "Naive/memory CD8",
  "Innate-like cytotoxic T",
  "Regulatory T",
  "Interferon-responsive T"
)

state_ids <- as.character(
  tcells_v3_curated$immune_state_id
)

unexpected_state_ids <- setdiff(
  unique(state_ids),
  names(major_family_map)
)

if (length(unexpected_state_ids) > 0) {
  stop(
    "No major-family definition was found for: ",
    paste(
      unexpected_state_ids,
      collapse = ", "
    )
  )
}

tcells_v3_curated$major_state_family <- unname(
  major_family_map[
    state_ids
  ]
)

if (
  anyNA(
    tcells_v3_curated$major_state_family
  )
) {
  stop(
    "Missing major-state family assignments."
  )
}

tcells_v3_curated$major_state_family <- factor(
  tcells_v3_curated$major_state_family,
  levels = major_family_levels
)

# Save the state-to-family mapping used in the figure.

major_family_assignment <- data.frame(
  immune_state_id = names(major_family_map),
  major_state_family = unname(
    major_family_map
  ),
  stringsAsFactors = FALSE
)

write.csv(
  major_family_assignment,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_major_state_family_mapping.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Reproducible plotting limits
# ------------------------------------------------------------

umap_coordinates <- Embeddings(
  tcells_v3_curated,
  reduction = "umap_harmony"
)

if (ncol(umap_coordinates) < 2) {
  stop(
    "The Harmony UMAP must contain at least two dimensions."
  )
}

# A very small number of extreme UMAP coordinates would otherwise
# compress the main atlas. The underlying cells are not removed;
# only the displayed coordinate window is restricted.

plot_quantiles <- c(
  0.0005,
  0.9995
)

x_limits <- as.numeric(
  quantile(
    umap_coordinates[, 1],
    probs = plot_quantiles,
    na.rm = TRUE,
    names = FALSE
  )
)

y_limits <- as.numeric(
  quantile(
    umap_coordinates[, 2],
    probs = plot_quantiles,
    na.rm = TRUE,
    names = FALSE
  )
)

x_padding <- 0.03 * diff(x_limits)
y_padding <- 0.03 * diff(y_limits)

x_limits <- c(
  x_limits[1] - x_padding,
  x_limits[2] + x_padding
)

y_limits <- c(
  y_limits[1] - y_padding,
  y_limits[2] + y_padding
)

plot_limits <- data.frame(
  dimension = c(
    "UMAP_1",
    "UMAP_2"
  ),
  lower_limit = c(
    x_limits[1],
    y_limits[1]
  ),
  upper_limit = c(
    x_limits[2],
    y_limits[2]
  ),
  lower_quantile = plot_quantiles[1],
  upper_quantile = plot_quantiles[2]
)

write.csv(
  plot_limits,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_curated_atlas_plot_limits.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Colour palettes
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

major_family_palette <- c(
  "CD4 T-cell states" =
    "#4C78A8",
  "Cytotoxic/effector CD8" =
    "#E15759",
  "Naive/memory CD8" =
    "#76B7B2",
  "Innate-like cytotoxic T" =
    "#E8B923",
  "Regulatory T" =
    "#B279A2",
  "Interferon-responsive T" =
    "#2F4858"
)

missing_state_colours <- setdiff(
  unique(state_ids),
  names(state_palette)
)

if (length(missing_state_colours) > 0) {
  stop(
    "Missing state colours for: ",
    paste(
      missing_state_colours,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Common publication theme
# ------------------------------------------------------------

atlas_theme <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      colour = "grey35",
      margin = margin(
        b = 6
      )
    ),
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      size = 9,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.35
    ),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

# ------------------------------------------------------------
# Panel A: curated immune states
# ------------------------------------------------------------

p_states_clean <- DimPlot(
  tcells_v3_curated,
  reduction = "umap_harmony",
  group.by = "immune_state_id",
  cols = state_palette,
  pt.size = 0.12,
  shuffle = TRUE,
  seed = 123,
  label = TRUE,
  repel = TRUE,
  label.size = 3.5,
  label.box = TRUE,
  raster = TRUE
) +
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE,
    clip = "on"
  ) +
  labs(
    title =
      "Curated ImmunoVerse T-cell states",
    subtitle =
      "Eleven transcriptionally and functionally defined immune states",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  atlas_theme +
  theme(
    legend.position = "none"
  )

# ------------------------------------------------------------
# Panel B: higher-order state families
# ------------------------------------------------------------

p_families_clean <- DimPlot(
  tcells_v3_curated,
  reduction = "umap_harmony",
  group.by = "major_state_family",
  cols = major_family_palette,
  pt.size = 0.12,
  shuffle = TRUE,
  seed = 123,
  label = FALSE,
  raster = TRUE
) +
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE,
    clip = "on"
  ) +
  labs(
    title =
      "Higher-order immune-state families",
    subtitle =
      "Functional organization of the curated T-cell atlas",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "State family"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  atlas_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.spacing.x = grid::unit(
      0.15,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.42,
      "cm"
    ),
    legend.key.width = grid::unit(
      0.45,
      "cm"
    ),
    legend.margin = margin(
      t = 4
    )
  )

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

p_curated_atlas_publication <- (
  p_states_clean |
    p_families_clean
) +
  plot_layout(
    widths = c(
      1,
      1.08
    )
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 15
      )
    )
  )

# Display in the active graphics device.

print(
  p_curated_atlas_publication
)

# ------------------------------------------------------------
# Save publication-oriented figure
# ------------------------------------------------------------

publication_pdf <- file.path(
  FIG_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas_publication_UMAP.pdf"
  )
)

publication_png <- file.path(
  FIG_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas_publication_UMAP.png"
  )
)

ggsave(
  filename = publication_pdf,
  plot = p_curated_atlas_publication,
  width = 15,
  height = 7.8,
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = publication_png,
  plot = p_curated_atlas_publication,
  width = 15,
  height = 7.8,
  dpi = 400,
  bg = "white"
)

# ------------------------------------------------------------
# Final validation and output
# ------------------------------------------------------------

major_family_counts <- tcells_v3_curated@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  count(
    immune_state_id,
    major_state_family,
    name = "n_cells"
  ) %>%
  arrange(
    immune_state_id
  )

write.csv(
  major_family_counts,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_major_state_family_counts.csv"
    )
  ),
  row.names = FALSE
)

message(
  "Publication-oriented curated atlas figure completed."
)

message(
  "PDF saved to: ",
  publication_pdf
)

message(
  "PNG saved to: ",
  publication_png
)

