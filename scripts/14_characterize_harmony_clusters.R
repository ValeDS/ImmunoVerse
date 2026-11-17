# ============================================================
# 14_characterize_harmony_clusters.R
# Biological and technical characterization of Harmony clusters
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
  "_Tcells_Harmony_reference"
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
    "Harmony reference not found: ",
    in_rds
  )
}

message("Loading Harmony-integrated reference...")

tcells_harmony <- readRDS(in_rds)
DefaultAssay(tcells_harmony) <- "RNA"

if (
  !"integrated_state_cluster" %in%
    colnames(tcells_harmony@meta.data)
) {
  stop(
    "integrated_state_cluster is missing."
  )
}

Idents(tcells_harmony) <-
  "integrated_state_cluster"

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

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
    "CD4", "IL7R", "LTB", "MAL"
  ),
  CD8 = c(
    "CD8A", "CD8B", "CCL5"
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

markers_present <- unique(
  unlist(marker_panel)
)

markers_present <- markers_present[
  markers_present %in%
    rownames(tcells_harmony)
]

marker_panel_present <- lapply(
  marker_panel,
  function(x) {
    x[x %in% markers_present]
  }
)

marker_panel_present <- marker_panel_present[
  lengths(marker_panel_present) > 0
]

# DotPlot does not allow the same feature to occur in more than
# one named feature group. Retain each gene in its first relevant
# category so that factor levels remain unique.
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

message(
  "Diagnostic markers found: ",
  length(
    unique(
      unlist(marker_panel_present)
    )
  )
)

# ------------------------------------------------------------
# Diagnostic dot plot
# ------------------------------------------------------------

p_marker_dot <- DotPlot(
  tcells_harmony,
  features = marker_panel_present,
  group.by =
    "integrated_state_cluster",
  assay = "RNA"
) +
  RotatedAxis() +
  labs(
    title =
      "Biological and technical marker profiles",
    x = NULL,
    y = "Harmony cluster",
    colour =
      "Average expression",
    size =
      "Percent expressed"
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
  p_marker_dot,
  width = 22,
  height = 9
)

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_diagnostic_marker_dotplot.png"
    )
  ),
  p_marker_dot,
  width = 22,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# Cluster marker discovery
# RNA expression is used; Harmony is used only for clustering.
# Downsampling within identity limits computational cost and
# prevents the largest clusters from dominating the test.
# ------------------------------------------------------------

message(
  "Finding positive markers for all Harmony clusters..."
)

harmony_markers <- FindAllMarkers(
  tcells_harmony,
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
  harmony_markers,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_all_positive_markers.csv"
    )
  ),
  row.names = FALSE
)

top10_markers <- harmony_markers %>%
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
  top10_markers,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_top10_markers.csv"
    )
  ),
  row.names = FALSE
)

marker_summary <- top10_markers %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  )

write.csv(
  marker_summary,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_marker_summary.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Focused review of residual RRMS/V2-dominated clusters
# ------------------------------------------------------------

residual_clusters <- c(
  "5", "7", "10", "13",
  "16", "17", "18", "19"
)

residual_marker_review <- top10_markers %>%
  filter(
    as.character(cluster) %in%
      residual_clusters
  )

write.csv(
  residual_marker_review,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_residual_RRMS_V2_cluster_markers.csv"
    )
  ),
  row.names = FALSE
)

print(
  tibble::as_tibble(
    marker_summary
  ),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(
    residual_marker_review
  ),
  n = Inf,
  width = Inf
)

message(
  "Harmony cluster characterization completed."
)
