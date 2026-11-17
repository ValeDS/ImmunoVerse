# ============================================================
# 32_extract_and_QC_Tcells_GSE138266.R
# Conservative extraction of T cells for external validation
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

INPUT_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_all_leukocytes_normalized_clustered.rds"
)

OUTPUT_DIR <- paste0(
  "data/processed/external_validation/GSE138266"
)
TABLE_DIR <- paste0(
  "results/tables/external_validation/GSE138266/",
  "Tcell_identification"
)
FIGURE_DIR <- paste0(
  "results/figures/external_validation/GSE138266"
)

OUTPUT_FILE <- file.path(
  OUTPUT_DIR,
  "GSE138266_Tcells_external_validation.rds"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input object not found: ",
    INPUT_FILE,
    ". Run script 31 first."
  )
}

message("Loading the GSE138266 broad-lineage object...")
external <- readRDS(INPUT_FILE)

required_metadata <- c(
  "sample",
  "donor",
  "group",
  "compartment",
  "external_leukocyte_cluster"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(external[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (!"umap_harmony" %in% Reductions(external)) {
  stop("The expected Harmony UMAP is absent from the input object.")
}

DefaultAssay(external) <- "RNA"

# ------------------------------------------------------------
# Cell-level evidence for canonical T-cell identity
# ------------------------------------------------------------

canonical_T_genes <- intersect(
  c("CD3D", "CD3E", "CD3G", "TRAC"),
  rownames(external)
)

if (length(canonical_T_genes) < 3L) {
  stop(
    "Fewer than three canonical T-cell genes are available."
  )
}

normalized_data <- LayerData(
  external,
  assay = "RNA",
  layer = "data"
)

external$canonical_T_gene_count <- as.integer(
  Matrix::colSums(
    normalized_data[
      canonical_T_genes,
      ,
      drop = FALSE
    ] > 0
  )
)

external$primary_T_inclusion <- (
  external$external_leukocyte_cluster %in% c("0", "1", "2") &
    external$canonical_T_gene_count >= 1L
) | (
  external$external_leukocyte_cluster == "4" &
    external$canonical_T_gene_count >= 2L
)

external$Tcell_selection_class <- case_when(
  external$external_leukocyte_cluster %in% c("0", "1", "2") &
    external$canonical_T_gene_count >= 1L ~
    "Included: conventional T-cell clusters",
  external$external_leukocyte_cluster == "4" &
    external$canonical_T_gene_count >= 2L ~
    "Included: CD3-supported T/NK cells",
  external$external_leukocyte_cluster == "4" ~
    "Excluded: predominantly NK-like cells",
  external$external_leukocyte_cluster == "11" ~
    "Excluded: platelet-associated T-cell cluster",
  external$external_leukocyte_cluster %in% c("0", "1", "2") ~
    "Excluded: insufficient canonical T evidence",
  TRUE ~
    "Excluded: non-T lineage"
)

rm(normalized_data)
gc(verbose = FALSE)

# ------------------------------------------------------------
# Selection audit
# ------------------------------------------------------------

selection_audit <- external[[]] %>%
  group_by(
    external_leukocyte_cluster,
    Tcell_selection_class,
    primary_T_inclusion
  ) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    n_donors = n_distinct(donor),
    median_canonical_T_gene_count =
      median(canonical_T_gene_count),
    .groups = "drop"
  ) %>%
  arrange(
    desc(primary_T_inclusion),
    as.integer(external_leukocyte_cluster),
    Tcell_selection_class
  )

included_metadata <- external[[]] %>%
  filter(primary_T_inclusion)

sample_audit <- included_metadata %>%
  count(
    sample,
    donor,
    group,
    compartment,
    name = "n_T_cells"
  ) %>%
  arrange(compartment, group, sample)

design_audit <- included_metadata %>%
  count(
    compartment,
    group,
    name = "n_T_cells"
  ) %>%
  left_join(
    included_metadata %>%
      distinct(
        sample,
        donor,
        group,
        compartment
      ) %>%
      group_by(compartment, group) %>%
      summarise(
        n_samples = n_distinct(sample),
        n_donors = n_distinct(donor),
        .groups = "drop"
      ),
    by = c("compartment", "group")
  )

selection_rules <- data.frame(
  rule = c(
    "Conventional T-cell inclusion",
    "T/NK boundary inclusion",
    "Platelet-associated exclusion"
  ),
  definition = c(
    paste0(
      "Broad clusters 0, 1, or 2 with expression of at least ",
      "one of CD3D, CD3E, CD3G, or TRAC"
    ),
    paste0(
      "Broad cluster 4 with expression of at least two of ",
      "CD3D, CD3E, CD3G, or TRAC"
    ),
    paste0(
      "Broad cluster 11 excluded because its leading markers ",
      "include PF4, PPBP, GNG11, and SDPR"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  selection_audit,
  file.path(
    TABLE_DIR,
    "GSE138266_Tcell_selection_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_audit,
  file.path(
    TABLE_DIR,
    "GSE138266_Tcell_counts_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  design_audit,
  file.path(
    TABLE_DIR,
    "GSE138266_Tcell_external_validation_design.csv"
  ),
  row.names = FALSE
)

write.csv(
  selection_rules,
  file.path(
    TABLE_DIR,
    "GSE138266_Tcell_selection_rules.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Diagnostic figure
# ------------------------------------------------------------

selection_colors <- c(
  "Included: conventional T-cell clusters" = "#4C78A8",
  "Included: CD3-supported T/NK cells" = "#7A5195",
  "Excluded: predominantly NK-like cells" = "#C7A9D6",
  "Excluded: platelet-associated T-cell cluster" = "#F2CF5B",
  "Excluded: insufficient canonical T evidence" = "#9D9D9D",
  "Excluded: non-T lineage" = "#E5E5E5"
)

external$Tcell_selection_class <- factor(
  external$Tcell_selection_class,
  levels = names(selection_colors)
)

panel_A <- DimPlot(
  external,
  reduction = "umap_harmony",
  group.by = "Tcell_selection_class",
  cols = selection_colors,
  raster = TRUE,
  pt.size = 0.12
) +
  ggtitle("Conservative T-cell selection") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

plot_design <- design_audit %>%
  mutate(
    condition = paste(compartment, group, sep = " | ")
  )

panel_B <- ggplot(
  plot_design,
  aes(
    x = condition,
    y = n_T_cells,
    fill = compartment
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label = paste0(
        scales::comma(n_T_cells),
        " cells\n",
        n_samples,
        " samples"
      )
    ),
    vjust = -0.25,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "PBMC" = "#4C78A8",
      "CSF" = "#E45756"
    )
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Selected T cells across the validation design",
    x = NULL,
    y = "Selected T cells",
    fill = "Compartment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

diagnostic_figure <- panel_A + panel_B +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_conservative_Tcell_selection.png"
  ),
  diagnostic_figure,
  width = 15,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_conservative_Tcell_selection.pdf"
  ),
  diagnostic_figure,
  width = 15,
  height = 6.5,
  device = grDevices::pdf,
  bg = "white"
)

# ------------------------------------------------------------
# Save the external-validation T-cell object
# ------------------------------------------------------------

selected_cells <- rownames(
  external[[]]
)[external$primary_T_inclusion]

if (length(selected_cells) == 0L) {
  stop("No cells satisfied the primary T-cell inclusion rules.")
}

t_cells <- subset(
  external,
  cells = selected_cells
)

t_cells$external_Tcell_selection <- "included"

# Retain the RNA counts and normalized values needed for label
# transfer while removing the broad-lineage reductions and graphs.
t_cells <- DietSeurat(
  t_cells,
  assays = "RNA",
  layers = c("counts", "data"),
  dimreducs = character(0),
  graphs = character(0),
  misc = TRUE
)

saveRDS(
  t_cells,
  OUTPUT_FILE,
  compress = FALSE
)

validation_summary <- data.frame(
  metric = c(
    "input_leukocytes",
    "selected_T_cells",
    "selected_samples",
    "selected_donors",
    "selected_PBMC_cells",
    "selected_CSF_cells",
    "canonical_T_genes_used"
  ),
  value = c(
    ncol(external),
    ncol(t_cells),
    n_distinct(t_cells$sample),
    n_distinct(t_cells$donor),
    sum(t_cells$compartment == "PBMC"),
    sum(t_cells$compartment == "CSF"),
    length(canonical_T_genes)
  )
)

write.csv(
  validation_summary,
  file.path(
    TABLE_DIR,
    "GSE138266_Tcell_extraction_validation_summary.csv"
  ),
  row.names = FALSE
)

message("")
message("GSE138266 conservative T-cell extraction completed.")
message("")
message("Validation summary:")
print(validation_summary, row.names = FALSE)
message("")
message("External-validation design:")
print(design_audit, row.names = FALSE)
message("")
message("Saved T-cell object: ", OUTPUT_FILE)

