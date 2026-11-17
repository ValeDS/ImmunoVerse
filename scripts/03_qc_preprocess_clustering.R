# ============================================================
# 03_qc_preprocess_clustering.R
# QC, preprocessing, PCA, clustering, UMAP, markers
# ============================================================

source("scripts/config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

in_rds <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_seurat_raw.rds"))
if (!file.exists(in_rds)) stop("Raw Seurat object not found. Run 02_create_seurat_object.R first.")

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"

# Mitochondrial percentage. Human genes usually start with MT-.
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

qc_before <- seu@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  summarise(
    n_cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt)
  )
readr::write_csv(qc_before, file.path(TAB_DIR, paste0(GEO_ACCESSION, "_qc_before_filtering.csv")))

p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
ggsave(file.path(FIG_DIR, paste0(GEO_ACCESSION, "_QC_violin_before_filtering.pdf")), p_qc, width = 12, height = 5)

message("Filtering cells with nFeature_RNA > ", MIN_FEATURES,
        ", nFeature_RNA < ", MAX_FEATURES,
        ", percent.mt < ", MAX_MITO)

seu <- subset(
  seu,
  subset = nFeature_RNA > MIN_FEATURES &
    nFeature_RNA < MAX_FEATURES &
    percent.mt < MAX_MITO
)

qc_after <- seu@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  summarise(
    n_cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt)
  )
readr::write_csv(qc_after, file.path(TAB_DIR, paste0(GEO_ACCESSION, "_qc_after_filtering.csv")))

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = max(N_PCS, 50))

p_elbow <- ElbowPlot(seu, ndims = 50)
ggsave(file.path(FIG_DIR, paste0(GEO_ACCESSION, "_elbow_plot.pdf")), p_elbow, width = 7, height = 5)

seu <- FindNeighbors(seu, dims = 1:N_PCS)
seu <- FindClusters(seu, resolution = CLUSTER_RESOLUTION)
seu <- RunUMAP(seu, dims = 1:N_PCS)

p_umap_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
ggsave(file.path(FIG_DIR, paste0(GEO_ACCESSION, "_UMAP_clusters.pdf")), p_umap_cluster, width = 8, height = 6)

if ("sample_folder" %in% colnames(seu@meta.data)) {
  p_umap_sample <- DimPlot(seu, reduction = "umap", group.by = "sample_folder") + ggtitle("Samples")
  ggsave(file.path(FIG_DIR, paste0(GEO_ACCESSION, "_UMAP_samples.pdf")), p_umap_sample, width = 10, height = 7)
}

# Canonical immune marker panel for preliminary annotation
marker_panel <- list(
  T_cells = c("CD3D", "CD3E", "TRAC"),
  CD4_T = c("CD4", "IL7R", "CCR7"),
  CD8_T = c("CD8A", "CD8B", "GZMK", "NKG7"),
  B_cells = c("MS4A1", "CD79A", "CD79B"),
  Plasma_cells = c("MZB1", "JCHAIN", "XBP1"),
  NK_cells = c("GNLY", "NKG7", "KLRD1"),
  Monocytes = c("LYZ", "S100A8", "S100A9", "FCN1"),
  Dendritic_cells = c("FCER1A", "CST3", "CLEC10A"),
  Platelets = c("PPBP", "PF4")
)

markers_present <- unique(unlist(marker_panel))
markers_present <- markers_present[markers_present %in% rownames(seu)]

if (length(markers_present) > 0) {
  p_dot <- DotPlot(seu, features = markers_present) + RotatedAxis()
  ggsave(file.path(FIG_DIR, paste0(GEO_ACCESSION, "_canonical_marker_dotplot.pdf")), p_dot, width = 14, height = 7)
}

# message("Finding cluster markers...")
# cluster_markers <- FindAllMarkers(
#  seu,
#  only.pos = TRUE,
#  min.pct = 0.25,
#  logfc.threshold = 0.25
# )
# readr::write_csv(cluster_markers, file.path(TAB_DIR, paste0(GEO_ACCESSION, "_cluster_markers.csv")))

out_rds <- file.path(PROCESSED_DIR, paste0(GEO_ACCESSION, "_seurat_qc_clustered.rds"))
saveRDS(seu, out_rds)
message("Saved processed Seurat object: ", out_rds)
message("QC and clustering step completed.")
