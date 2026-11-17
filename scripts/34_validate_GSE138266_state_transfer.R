# ============================================================
# 34_validate_GSE138266_state_transfer.R
# QC, recurrence, marker preservation, and visualization
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

QUERY_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_Tcells_external_validation.rds"
)

PREDICTION_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_curated_state_predictions.rds"
)

MAPPED_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_Tcells_curated_state_mapped.rds"
)

TABLE_DIR <- paste0(
  "results/tables/external_validation/GSE138266/",
  "state_transfer"
)

FIGURE_DIR <- paste0(
  "results/figures/external_validation/GSE138266"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(QUERY_FILE, PREDICTION_FILE)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Loading external T cells and transferred-state predictions...")
t_cells <- readRDS(QUERY_FILE)
predictions <- readRDS(PREDICTION_FILE)

if (!inherits(t_cells, "Seurat") || !is.data.frame(predictions)) {
  stop("The input objects have unexpected classes.")
}

if (
  anyDuplicated(predictions$cell_id) ||
    !setequal(predictions$cell_id, colnames(t_cells))
) {
  stop("Prediction cell identifiers do not match the T-cell object.")
}

predictions <- predictions[
  match(colnames(t_cells), predictions$cell_id),
  ,
  drop = FALSE
]

if (!identical(predictions$cell_id, colnames(t_cells))) {
  stop("Unable to align predictions to the T-cell object.")
}

state_ids <- sprintf("IVS-%02d", 1:11)

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

higher_order_family_levels <- c(
  "CD4 T-cell states",
  "Cytotoxic/effector CD8",
  "Naive/memory CD8",
  "Innate-like cytotoxic T",
  "Regulatory T",
  "Interferon-responsive T"
)

mapped_metadata <- data.frame(
  predicted_state_id = factor(
    as.character(predictions$predicted.id),
    levels = state_ids
  ),
  predicted_state_name = predictions$immune_state_name,
  reference_state_family = predictions$state_family,
  higher_order_family = factor(
    unname(
      higher_order_family_map[
        as.character(predictions$predicted.id)
      ]
    ),
    levels = higher_order_family_levels
  ),
  prediction_score = predictions$prediction.score.max,
  stringsAsFactors = FALSE,
  row.names = predictions$cell_id
)

if (anyNA(mapped_metadata$higher_order_family)) {
  stop("One or more transferred states lack a higher-order family.")
}

t_cells <- AddMetaData(t_cells, metadata = mapped_metadata)
DefaultAssay(t_cells) <- "RNA"

# ------------------------------------------------------------
# T-cell-specific representation for visualization only
# ------------------------------------------------------------

if (!all(c("pca_external", "harmony_external", "umap_external") %in%
  Reductions(t_cells))) {

  message("Building a T-cell-specific visualization representation...")

  t_cells <- FindVariableFeatures(
    t_cells,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  variable_features <- VariableFeatures(t_cells)
  variable_features <- variable_features[
    !grepl("^(MT-|RPL|RPS|HBA|HBB|HBD)", variable_features)
  ]
  VariableFeatures(t_cells) <- variable_features

  t_cells <- ScaleData(
    t_cells,
    features = variable_features,
    verbose = FALSE
  )

  t_cells <- RunPCA(
    t_cells,
    features = variable_features,
    npcs = 30,
    reduction.name = "pca_external",
    reduction.key = "PCAext_",
    verbose = FALSE
  )

  # Harmony is used only to produce a sample-aware visualization.
  # All abundance tests use the original sample-level counts.
  t_cells <- RunHarmony(
    object = t_cells,
    group.by.vars = "sample",
    reduction.use = "pca_external",
    dims.use = 1:30,
    reduction.save = "harmony_external",
    project.dim = FALSE,
    verbose = FALSE
  )

  t_cells <- RunUMAP(
    t_cells,
    reduction = "harmony_external",
    dims = 1:30,
    reduction.name = "umap_external",
    reduction.key = "UMAPext_",
    seed.use = 123,
    verbose = FALSE
  )
}

saveRDS(t_cells, MAPPED_FILE, compress = FALSE)

# ------------------------------------------------------------
# Recurrence, composition, and confidence tables
# ------------------------------------------------------------

cell_metadata <- t_cells[[]] %>%
  mutate(cell_id_join = rownames(t_cells[[]]))

state_composition <- cell_metadata %>%
  count(
    sample,
    donor,
    group,
    compartment,
    predicted_state_id,
    name = "n_cells"
  ) %>%
  group_by(sample) %>%
  complete(
    predicted_state_id = factor(state_ids, levels = state_ids),
    fill = list(n_cells = 0L)
  ) %>%
  mutate(
    percent = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

family_composition <- cell_metadata %>%
  count(
    sample,
    donor,
    group,
    compartment,
    higher_order_family,
    name = "n_cells"
  ) %>%
  group_by(sample) %>%
  complete(
    higher_order_family = factor(
      higher_order_family_levels,
      levels = higher_order_family_levels
    ),
    fill = list(n_cells = 0L)
  ) %>%
  mutate(
    percent = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

state_recurrence <- cell_metadata %>%
  group_by(
    predicted_state_id,
    predicted_state_name,
    higher_order_family
  ) %>%
  summarise(
    n_cells = n(),
    percent_all_T_cells = 100 * n() / nrow(cell_metadata),
    samples_with_cells = n_distinct(sample),
    donors_with_cells = n_distinct(donor),
    PBMC_samples_with_cells = n_distinct(sample[compartment == "PBMC"]),
    CSF_samples_with_cells = n_distinct(sample[compartment == "CSF"]),
    median_score = median(prediction_score),
    .groups = "drop"
  ) %>%
  complete(
    predicted_state_id = factor(state_ids, levels = state_ids),
    fill = list(
      n_cells = 0L,
      percent_all_T_cells = 0,
      samples_with_cells = 0L,
      donors_with_cells = 0L,
      PBMC_samples_with_cells = 0L,
      CSF_samples_with_cells = 0L
    )
  ) %>%
  mutate(
    predicted_state_name = ifelse(
      as.character(predicted_state_id) == "IVS-02" &
        is.na(predicted_state_name),
      "ANKRD55+ naive/central-memory CD4 state",
      predicted_state_name
    ),
    higher_order_family = factor(
      ifelse(
        as.character(predicted_state_id) == "IVS-02" &
          is.na(as.character(higher_order_family)),
        "CD4 T-cell states",
        as.character(higher_order_family)
      ),
      levels = higher_order_family_levels
    )
  ) %>%
  arrange(predicted_state_id)

compartment_confidence <- cell_metadata %>%
  group_by(compartment) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(prediction_score),
    median_score = median(prediction_score),
    q25_score = quantile(prediction_score, 0.25, names = FALSE),
    q75_score = quantile(prediction_score, 0.75, names = FALSE),
    pct_below_050 = 100 * mean(prediction_score < 0.50),
    pct_below_075 = 100 * mean(prediction_score < 0.75),
    pct_below_090 = 100 * mean(prediction_score < 0.90),
    .groups = "drop"
  )

sample_confidence <- cell_metadata %>%
  group_by(sample, donor, group, compartment) %>%
  summarise(
    n_cells = n(),
    median_score = median(prediction_score),
    mean_score = mean(prediction_score),
    .groups = "drop"
  )

write.csv(
  state_composition,
  file.path(TABLE_DIR, "GSE138266_state_composition_by_sample.csv"),
  row.names = FALSE
)

write.csv(
  family_composition,
  file.path(TABLE_DIR, "GSE138266_family_composition_by_sample.csv"),
  row.names = FALSE
)

write.csv(
  state_recurrence,
  file.path(TABLE_DIR, "GSE138266_state_recurrence.csv"),
  row.names = FALSE
)

write.csv(
  compartment_confidence,
  file.path(TABLE_DIR, "GSE138266_transfer_confidence_by_compartment.csv"),
  row.names = FALSE
)

write.csv(
  sample_confidence,
  file.path(TABLE_DIR, "GSE138266_transfer_confidence_by_sample.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Marker-program preservation
# ------------------------------------------------------------

marker_panels <- list(
  "IVS-01" = c("CCR7", "TCF7", "LEF1", "MAL"),
  "IVS-02" = c("ANKRD55", "LRRN3", "IGF1R", "NOG"),
  "IVS-03" = c("MAF", "LIMS1", "GPR183", "CD82"),
  "IVS-04" = c("TNFRSF4", "TNFRSF18", "FUT7", "LGALS3"),
  "IVS-05" = c("CX3CR1", "FGFBP2", "GZMH", "GZMB"),
  "IVS-06" = c("CD27", "CD8A", "CD8B", "LTB"),
  "IVS-07" = c("GZMK", "EOMES", "CCL5", "XCL2"),
  "IVS-08" = c("KLRB1", "SLC4A10", "ZBTB16", "IL18RAP"),
  "IVS-09" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "IVS-10" = c("KLRC1", "FCER1G", "TYROBP", "NKG7"),
  "IVS-11" = c("IFIT1", "IFIT3", "IFI44L", "ISG15")
)

marker_availability <- bind_rows(lapply(
  names(marker_panels),
  function(state_id) {
    data.frame(
      state_id = state_id,
      gene = marker_panels[[state_id]],
      available = marker_panels[[state_id]] %in% rownames(t_cells),
      stringsAsFactors = FALSE
    )
  }
))

write.csv(
  marker_availability,
  file.path(TABLE_DIR, "GSE138266_marker_panel_availability.csv"),
  row.names = FALSE
)

available_marker_panels <- lapply(
  marker_panels,
  intersect,
  y = rownames(t_cells)
)

available_marker_panels <- available_marker_panels[
  lengths(available_marker_panels) > 0L
]

Idents(t_cells) <- "predicted_state_id"

marker_dotplot <- DotPlot(
  t_cells,
  features = available_marker_panels,
  assay = "RNA",
  group.by = "predicted_state_id",
  dot.scale = 5
) +
  RotatedAxis() +
  labs(
    title = "Marker-program preservation after external transfer",
    subtitle = paste0(
      "Rows are transferred states; marker groups originate from ",
      "the curated GSE144744 reference"
    ),
    x = NULL,
    y = "Transferred immune state"
  ) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_external_transfer_marker_preservation.png"
  ),
  marker_dotplot,
  width = 16,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_external_transfer_marker_preservation.pdf"
  ),
  marker_dotplot,
  width = 16,
  height = 6.5,
  device = grDevices::pdf,
  bg = "white"
)

# ------------------------------------------------------------
# External-transfer summary figure
# ------------------------------------------------------------

state_colors <- c(
  "IVS-01" = "#3B6FB6",
  "IVS-02" = "#77A6D1",
  "IVS-03" = "#4E9F62",
  "IVS-04" = "#F28E2B",
  "IVS-05" = "#E15759",
  "IVS-06" = "#76B7B2",
  "IVS-07" = "#9C6BAE",
  "IVS-08" = "#E8B51A",
  "IVS-09" = "#C76E9A",
  "IVS-10" = "#8C564B",
  "IVS-11" = "#2F4858"
)

family_colors <- c(
  "CD4 T-cell states" = "#4C78A8",
  "Cytotoxic/effector CD8" = "#E15759",
  "Naive/memory CD8" = "#76B7B2",
  "Innate-like cytotoxic T" = "#B07AA1",
  "Regulatory T" = "#C76E9A",
  "Interferon-responsive T" = "#2F4858"
)

panel_A <- ggplot(
  cell_metadata,
  aes(
    x = prediction_score,
    color = compartment,
    fill = compartment
  )
) +
  geom_density(alpha = 0.20, linewidth = 0.8) +
  geom_vline(
    xintercept = 0.50,
    linetype = "dashed",
    color = "#B33A3A"
  ) +
  scale_color_manual(values = c(PBMC = "#4C78A8", CSF = "#E45756")) +
  scale_fill_manual(values = c(PBMC = "#4C78A8", CSF = "#E45756")) +
  labs(
    title = "Transfer confidence differs by compartment",
    subtitle = "The curated reference was constructed from peripheral blood",
    x = "Prediction score",
    y = "Density",
    color = "Compartment",
    fill = "Compartment"
  ) +
  theme_classic(base_size = 12)

panel_B <- ggplot(
  sample_confidence,
  aes(
    x = interaction(compartment, group, sep = " | "),
    y = median_score,
    color = compartment
  )
) +
  geom_boxplot(
    aes(group = interaction(compartment, group)),
    outlier.shape = NA,
    color = "#555555",
    fill = "white",
    width = 0.6
  ) +
  geom_point(size = 2.5, position = position_jitter(width = 0.08)) +
  geom_hline(
    yintercept = 0.50,
    linetype = "dashed",
    color = "#B33A3A"
  ) +
  scale_color_manual(values = c(PBMC = "#4C78A8", CSF = "#E45756")) +
  labs(
    title = "Sample-level mapping confidence",
    x = NULL,
    y = "Median prediction score",
    color = "Compartment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

panel_C <- DimPlot(
  t_cells,
  reduction = "umap_external",
  group.by = "predicted_state_id",
  cols = state_colors,
  raster = TRUE,
  pt.size = 0.15,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Transferred ImmunoVerse states") +
  NoLegend()

panel_D <- DimPlot(
  t_cells,
  reduction = "umap_external",
  group.by = "higher_order_family",
  cols = family_colors,
  raster = TRUE,
  pt.size = 0.15
) +
  ggtitle("Higher-order immune-state families") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

summary_figure <- (
  panel_A + panel_B
) / (
  panel_C + panel_D
) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_external_state_transfer_validation.png"
  ),
  summary_figure,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_external_state_transfer_validation.pdf"
  ),
  summary_figure,
  width = 16,
  height = 12,
  device = grDevices::pdf,
  bg = "white"
)

message("")
message("GSE138266 external-transfer validation completed.")
message("")
message("Confidence by compartment:")
print(compartment_confidence, row.names = FALSE)
message("")
message("State recurrence:")
print(state_recurrence, n = Inf)
message("")
message("Saved mapped object: ", MAPPED_FILE)
