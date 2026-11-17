# ============================================================
# 16_characterize_V3_harmony_clusters.R
# Biological characterization of the V3 Harmony reference
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    ".rds"
  )
)

if (!file.exists(in_rds)) {
  stop(
    "V3 Harmony reference not found: ",
    in_rds
  )
}

message("Loading V3 Harmony reference...")

tcells_v3 <- readRDS(in_rds)
DefaultAssay(tcells_v3) <- "RNA"

if (
  !"V3_state_cluster" %in%
    colnames(tcells_v3@meta.data)
) {
  stop("V3_state_cluster is missing.")
}

Idents(tcells_v3) <- "V3_state_cluster"

# ------------------------------------------------------------
# Canonical and diagnostic marker panel
# ------------------------------------------------------------

marker_panel <- list(
  Pan_T = c(
    "CD3D", "CD3E", "TRAC"
  ),
  Naive_central_memory = c(
    "CCR7", "SELL", "TCF7", "LEF1",
    "MAL", "LTB", "IL7R"
  ),
  CD4 = c(
    "CD4"
  ),
  CD8 = c(
    "CD8A", "CD8B"
  ),
  Cytotoxic_effector = c(
    "NKG7", "GNLY", "PRF1", "GZMB",
    "GZMH", "FGFBP2", "CX3CR1"
  ),
  Memory_effector = c(
    "GZMK", "CCL5", "CST7", "DUSP2"
  ),
  Innate_like_MAIT = c(
    "KLRB1", "SLC4A10", "ZBTB16",
    "TRAV1-2", "NCR3"
  ),
  Regulatory = c(
    "FOXP3", "IL2RA", "CTLA4",
    "TIGIT", "IKZF2", "TNFRSF1B"
  ),
  Activation = c(
    "CD69", "HLA-DRA", "HLA-DRB1",
    "NR4A1", "NR4A3", "BHLHE40"
  ),
  Immediate_early_stress = c(
    "FOS", "JUN", "JUNB", "DUSP1",
    "EGR1", "HSPA1A"
  ),
  Interferon_response = c(
    "ISG15", "IFIT1", "IFIT3",
    "MX1", "IFI6", "IFI44L"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "STMN1",
    "TYMS"
  ),
  Myeloid_contamination = c(
    "LYZ", "S100A8", "S100A9",
    "FCN1", "CTSS"
  ),
  B_cell_contamination = c(
    "MS4A1", "CD79A", "CD79B",
    "CD37"
  ),
  Platelet_contamination = c(
    "PPBP", "PF4", "TUBB1",
    "TREML1"
  )
)

marker_panel_present <- lapply(
  marker_panel,
  function(x) {
    x[x %in% rownames(tcells_v3)]
  }
)

# DotPlot requires every feature to occur only once.
seen_markers <- character()

marker_panel_present <- lapply(
  marker_panel_present,
  function(x) {
    unique_markers <- x[
      !x %in% seen_markers
    ]
    
    seen_markers <<- c(
      seen_markers,
      unique_markers
    )
    
    unique_markers
  }
)

marker_panel_present <- marker_panel_present[
  lengths(marker_panel_present) > 0
]

# ------------------------------------------------------------
# Diagnostic marker dot plot
# ------------------------------------------------------------

p_v3_marker_dot <- DotPlot(
  tcells_v3,
  features = marker_panel_present,
  group.by = "V3_state_cluster",
  assay = "RNA"
) +
  RotatedAxis() +
  labs(
    title =
      "V3 reference: biological and technical marker profiles",
    x = NULL,
    y = "V3 Harmony cluster",
    colour = "Average expression",
    size = "Percent expressed"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 8
    ),
    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_diagnostic_marker_dotplot.pdf"
    )
  ),
  p_v3_marker_dot,
  width = 22,
  height = 8
)

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_diagnostic_marker_dotplot.png"
    )
  ),
  p_v3_marker_dot,
  width = 22,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Positive markers for each V3 cluster
# ------------------------------------------------------------

message(
  "Finding positive markers for V3 Harmony clusters..."
)

v3_markers <- FindAllMarkers(
  tcells_v3,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  max.cells.per.ident = 3000,
  random.seed = 123,
  verbose = TRUE
)

write.csv(
  v3_markers,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_all_positive_markers.csv"
    )
  ),
  row.names = FALSE
)

v3_top10_markers <- v3_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(
      as.character(cluster)
    ),
    desc(avg_log2FC)
  )

write.csv(
  v3_top10_markers,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_top10_markers.csv"
    )
  ),
  row.names = FALSE
)

v3_marker_summary <- v3_top10_markers %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  )

write.csv(
  v3_marker_summary,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_marker_summary.csv"
    )
  ),
  row.names = FALSE
)

print(
  tibble::as_tibble(
    v3_marker_summary
  ),
  n = Inf,
  width = Inf
)

message(
  "V3 Harmony cluster characterization completed."
)
