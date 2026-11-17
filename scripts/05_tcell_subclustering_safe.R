# ============================================================
# 05_tcell_subclustering_safe.R
# Safe T-cell subclustering with stratified downsampling
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

set.seed(123)

message("Starting SAFE T-cell subclustering analysis...")

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

MAX_CELLS_PER_GROUP <- 20000
T_RESOLUTION <- 0.8
T_NPCS <- 30

# ------------------------------------------------------------
# Load global processed object
# ------------------------------------------------------------

in_rds <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_seurat_qc_clustered.rds"))

if (!file.exists(in_rds)) {
  stop("Processed Seurat object not found: ", in_rds)
}

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"

message("Loaded global object:")
print(seu)

# ------------------------------------------------------------
# Subset T cells
# ------------------------------------------------------------

if (!"basictype" %in% colnames(seu@meta.data)) {
  stop("Column 'basictype' not found.")
}

tcells <- subset(seu, subset = basictype == "t_cells")

message("Extracted T cells:")
print(tcells)

saveRDS(
  tcells,
  file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_tcells_raw.rds"))
)

# ------------------------------------------------------------
# Stratified downsampling by group
# ------------------------------------------------------------

if (!"group" %in% colnames(tcells@meta.data)) {
  stop("Column 'group' not found.")
}

cell_meta <- tcells@meta.data %>%
  as_tibble(rownames = "cell_id")

  selected_cells <- cell_meta %>%
    group_by(group) %>%
    group_modify(~ {
      n_to_sample <- min(nrow(.x), MAX_CELLS_PER_GROUP)
      slice_sample(.x, n = n_to_sample)
    }) %>%
    ungroup() %>%
    pull(cell_id)

message("Selected cells after stratified downsampling: ", length(selected_cells))

tcells_ds <- subset(tcells, cells = selected_cells)

message("Downsampled T-cell object:")
print(tcells_ds)

# Save composition tables
write.csv(
  as.data.frame(table(tcells_ds$group)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_downsampled_group_counts.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(tcells_ds$cohort)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_downsampled_cohort_counts.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------
# Clear memory
# ------------------------------------------------------------

rm(seu, tcells, cell_meta)
gc()

# ------------------------------------------------------------
# Recompute T-cell representation
# ------------------------------------------------------------

DefaultAssay(tcells_ds) <- "RNA"

message("Normalizing downsampled T cells...")

tcells_ds <- NormalizeData(
  tcells_ds,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

message("Finding variable features...")

tcells_ds <- FindVariableFeatures(
  tcells_ds,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES
)

message("Scaling data...")

tcells_ds <- ScaleData(
  tcells_ds,
  features = VariableFeatures(tcells_ds)
)

message("Running PCA...")

tcells_ds <- RunPCA(
  tcells_ds,
  features = VariableFeatures(tcells_ds),
  npcs = 50
)

p_elbow <- ElbowPlot(tcells_ds, ndims = 50)

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_elbow_plot.pdf")),
  p_elbow,
  width = 7,
  height = 5
)

# ------------------------------------------------------------
# Neighbors, clustering, UMAP
# ------------------------------------------------------------

message("Computing nearest neighbors...")

tcells_ds <- FindNeighbors(
  tcells_ds,
  dims = 1:T_NPCS
)

message("Clustering T cells...")

tcells_ds <- FindClusters(
  tcells_ds,
  resolution = T_RESOLUTION
)

message("Running UMAP...")

tcells_ds <- RunUMAP(
  tcells_ds,
  dims = 1:T_NPCS
)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

p_clusters <- DimPlot(
  tcells_ds,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("T-cell subclusters")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_UMAP_subclusters.pdf")),
  p_clusters,
  width = 9,
  height = 7
)

p_group <- DimPlot(
  tcells_ds,
  reduction = "umap",
  group.by = "group"
) +
  ggtitle("T cells by group")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_UMAP_group.pdf")),
  p_group,
  width = 9,
  height = 7
)

p_cohort <- DimPlot(
  tcells_ds,
  reduction = "umap",
  group.by = "cohort"
) +
  ggtitle("T cells by cohort")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_UMAP_cohort.pdf")),
  p_cohort,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# T-cell marker dotplot
# ------------------------------------------------------------

tcell_marker_panel <- list(
  Pan_T = c("CD3D", "CD3E", "TRAC", "CD247"),
  CD4_T = c("CD4", "IL7R", "CCR7", "LTB", "SELL", "LEF1"),
  CD8_T = c("CD8A", "CD8B"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "GZMK", "CCL5"),
  Activation = c("CD69", "IL2RA", "CD40LG", "HLA-DRA", "HLA-DRB1"),
  Exhaustion = c("PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2"),
  Treg = c("FOXP3", "IL2RA", "IKZF2", "CTLA4"),
  CNS_homing = c("KLRB1", "LTB", "CCR6", "CXCR3", "CXCR4", "ITGA4"),
  Memory = c("IL7R", "CCR7", "SELL", "TCF7", "LEF1")
)

tcell_markers <- unique(unlist(tcell_marker_panel))
tcell_markers <- tcell_markers[tcell_markers %in% rownames(tcells_ds)]

message("T-cell markers found:")
print(tcell_markers)

if (length(tcell_markers) > 0) {
  p_dot <- DotPlot(
    tcells_ds,
    features = tcell_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Canonical T-cell state markers")

  ggsave(
    file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_marker_dotplot.pdf")),
    p_dot,
    width = 18,
    height = 8
  )
}

# ------------------------------------------------------------
# Save cluster tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(table(tcells_ds$seurat_clusters)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_cluster_counts.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(table(tcells_ds$seurat_clusters, tcells_ds$group)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_clusters_vs_group.csv"))
)

write.csv(
  round(prop.table(table(tcells_ds$seurat_clusters, tcells_ds$group), margin = 2) * 100, 2),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_clusters_vs_group_percent.csv"))
)

write.csv(
  as.data.frame.matrix(table(tcells_ds$seurat_clusters, tcells_ds$cohort)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_clusters_vs_cohort.csv"))
)

write.csv(
  round(prop.table(table(tcells_ds$seurat_clusters, tcells_ds$cohort), margin = 2) * 100, 2),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_safe_clusters_vs_cohort_percent.csv"))
)

# ------------------------------------------------------------
# Save final object
# ------------------------------------------------------------

out_rds <- file.path(
  PROCESSED_DIR,
  paste0(GEO_ACCESSION, "_tcells_safe_clustered.rds")
)

saveRDS(tcells_ds, out_rds)

message("Saved safe T-cell clustered object: ", out_rds)
message("SAFE T-cell subclustering analysis completed.")
