# ============================================================
# 05_tcell_subclustering.R
# T-cell focused subclustering for ImmunoVerse MS PoC
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

message("Starting T-cell subclustering analysis...")

# ------------------------------------------------------------
# 1. Load global processed Seurat object
# ------------------------------------------------------------

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(GEO_ACCESSION, "_seurat_qc_clustered.rds")
)

if (!file.exists(in_rds)) {
  stop("Processed Seurat object not found: ", in_rds)
}

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"

message("Loaded global object:")
print(seu)

# ------------------------------------------------------------
# 2. Subset T cells
# ------------------------------------------------------------

if (!"basictype" %in% colnames(seu@meta.data)) {
  stop("Column 'basictype' not found in metadata.")
}

tcells <- subset(
  seu,
  subset = basictype == "t_cells"
)

message("T-cell object:")
print(tcells)

# Save T-cell composition tables
tcell_summary <- tibble(
  n_tcells = ncol(tcells),
  n_genes = nrow(tcells),
  n_donors = length(unique(tcells$donor)),
  n_samples = length(unique(tcells$sample)),
  groups = paste(names(table(tcells$group)), collapse = ";"),
  cohorts = paste(names(table(tcells$cohort)), collapse = ";")
)

write_csv(
  tcell_summary,
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_tcell_summary.csv"))
)

write.csv(
  as.data.frame(table(tcells$group)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_tcell_group_counts.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(tcells$cohort)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_tcell_cohort_counts.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Recompute variable features, scaling, PCA
# ------------------------------------------------------------

message("Finding variable features in T cells...")

tcells <- FindVariableFeatures(
  tcells,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES
)

message("Scaling T-cell data...")

tcells <- ScaleData(
  tcells,
  features = VariableFeatures(tcells)
)

message("Running PCA on T cells...")

tcells <- RunPCA(
  tcells,
  features = VariableFeatures(tcells),
  npcs = 50
)

p_elbow <- ElbowPlot(tcells, ndims = 50)

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_elbow_plot.pdf")),
  p_elbow,
  width = 7,
  height = 5
)

# ------------------------------------------------------------
# 4. T-cell clustering and UMAP
# ------------------------------------------------------------

message("Computing T-cell neighbors...")

tcells <- FindNeighbors(
  tcells,
  dims = 1:30
)

message("Clustering T cells...")

tcells <- FindClusters(
  tcells,
  resolution = 0.8
)

message("Running UMAP for T cells...")

tcells <- RunUMAP(
  tcells,
  dims = 1:30
)

# ------------------------------------------------------------
# 5. Save UMAP plots
# ------------------------------------------------------------

p_tcluster <- DimPlot(
  tcells,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("T-cell subclusters")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_UMAP_subclusters.pdf")),
  p_tcluster,
  width = 9,
  height = 7
)

p_tgroup <- DimPlot(
  tcells,
  reduction = "umap",
  group.by = "group"
) +
  ggtitle("T cells by group")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_UMAP_group.pdf")),
  p_tgroup,
  width = 9,
  height = 7
)

p_tcohort <- DimPlot(
  tcells,
  reduction = "umap",
  group.by = "cohort"
) +
  ggtitle("T cells by cohort")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_UMAP_cohort.pdf")),
  p_tcohort,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# 6. T-cell marker dotplot
# ------------------------------------------------------------

tcell_marker_panel <- list(
  Pan_T = c("CD3D", "CD3E", "TRAC", "CD247"),
  CD4_T = c("CD4", "IL7R", "CCR7", "LTB", "SELL", "LEF1"),
  CD8_T = c("CD8A", "CD8B"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "GZMK", "CCL5"),
  Activation = c("CD69", "IL2RA", "CD40LG", "HLA-DRA", "HLA-DRB1"),
  Exhaustion = c("PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2"),
  Treg = c("FOXP3", "IL2RA", "IKZF2", "CTLA4"),
  Th17_CNS_homing = c("KLRB1", "LTB", "CCR6", "CXCR3", "CXCR4", "ITGA4"),
  Memory = c("IL7R", "CCR7", "SELL", "TCF7", "LEF1")
)

tcell_markers <- unique(unlist(tcell_marker_panel))
tcell_markers <- tcell_markers[tcell_markers %in% rownames(tcells)]

message("Markers found for T-cell dotplot:")
print(tcell_markers)

if (length(tcell_markers) > 0) {
  p_tdot <- DotPlot(
    tcells,
    features = tcell_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Canonical T-cell state markers")

  ggsave(
    file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_marker_dotplot.pdf")),
    p_tdot,
    width = 18,
    height = 8
  )
}

# ------------------------------------------------------------
# 7. Save T-cell object
# ------------------------------------------------------------

out_rds <- file.path(
  PROCESSED_DIR,
  paste0(GEO_ACCESSION, "_tcells_clustered.rds")
)

saveRDS(tcells, out_rds)

message("Saved T-cell clustered object: ", out_rds)

# ------------------------------------------------------------
# 8. Save cluster tables
# ------------------------------------------------------------

write.csv(
  as.data.frame.matrix(table(tcells$seurat_clusters)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_cluster_counts.csv"))
)

write.csv(
  as.data.frame.matrix(table(tcells$seurat_clusters, tcells$group)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_clusters_vs_group.csv"))
)

write.csv(
  round(prop.table(table(tcells$seurat_clusters, tcells$group), margin = 2) * 100, 2),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_clusters_vs_group_percent.csv"))
)

write.csv(
  as.data.frame.matrix(table(tcells$seurat_clusters, tcells$cohort)),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_clusters_vs_cohort.csv"))
)

write.csv(
  round(prop.table(table(tcells$seurat_clusters, tcells$cohort), margin = 2) * 100, 2),
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_clusters_vs_cohort_percent.csv"))
)

message("T-cell subclustering analysis completed.")
