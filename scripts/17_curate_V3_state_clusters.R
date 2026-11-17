# ============================================================
# 17_curate_V3_state_clusters.R
# Quantitative curation of the V3 Harmony clusters
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(OUT_PREFIX, ".rds")
)

if (!file.exists(in_rds)) {
  stop("V3 Harmony reference not found: ", in_rds)
}

message("Loading V3 Harmony reference...")

tcells_v3 <- readRDS(in_rds)
DefaultAssay(tcells_v3) <- "RNA"
Idents(tcells_v3) <- "V3_state_cluster"

if (!"percent.mt" %in% colnames(tcells_v3@meta.data)) {
  tcells_v3[["percent.mt"]] <- PercentageFeatureSet(
    tcells_v3,
    pattern = "^MT-"
  )
}

stopifnot(
  identical(
    rownames(tcells_v3@meta.data),
    colnames(tcells_v3)
  )
)

# ------------------------------------------------------------
# Normalized expression matrix
# ------------------------------------------------------------

rna_data <- LayerData(
  tcells_v3[["RNA"]],
  layer = "data"
)

# ------------------------------------------------------------
# Diagnostic expression programs
# ------------------------------------------------------------

diagnostic_programs <- list(
  Pan_T = c(
    "CD3D", "CD3E", "TRAC", "CD247"
  ),
  NK = c(
    "NCR1", "KLRD1", "KLRC1", "FCER1G",
    "TYROBP", "TRDC", "GNLY"
  ),
  Myeloid = c(
    "LYZ", "S100A8", "S100A9", "S100A12",
    "FCN1", "VCAN", "CLEC12A"
  ),
  Platelet = c(
    "PF4", "PF4V1", "PPBP", "GNG11",
    "CAVIN2", "TUBB1", "RGS18"
  ),
  Erythroid = c(
    "HBA1", "HBA2", "HBB", "HBD"
  )
)

mean_program_expression <- function(genes) {
  genes_present <- intersect(
    genes,
    rownames(rna_data)
  )
  
  if (length(genes_present) == 0) {
    return(rep(NA_real_, ncol(tcells_v3)))
  }
  
  Matrix::colMeans(
    rna_data[genes_present, , drop = FALSE]
  )
}

program_scores <- lapply(
  diagnostic_programs,
  mean_program_expression
)

program_scores <- as.data.frame(
  program_scores,
  check.names = FALSE
)

colnames(program_scores) <- paste0(
  "diagnostic_",
  colnames(program_scores)
)

tcells_v3@meta.data <- bind_cols(
  tcells_v3@meta.data,
  program_scores
)

# ------------------------------------------------------------
# Cluster-level quality-control summary
# ------------------------------------------------------------

curation_qc_summary <- tcells_v3@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    )
  ) %>%
  group_by(V3_state_cluster) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    n_donors = n_distinct(donor),
    median_nCount_RNA = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    median_nFeature_RNA = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    median_percent_mt = median(
      percent.mt,
      na.rm = TRUE
    ),
    median_pan_T = median(
      diagnostic_Pan_T,
      na.rm = TRUE
    ),
    median_NK = median(
      diagnostic_NK,
      na.rm = TRUE
    ),
    median_myeloid = median(
      diagnostic_Myeloid,
      na.rm = TRUE
    ),
    median_platelet = median(
      diagnostic_Platelet,
      na.rm = TRUE
    ),
    median_erythroid = median(
      diagnostic_Erythroid,
      na.rm = TRUE
    ),
    pct_pan_T_zero = 100 * mean(
      diagnostic_Pan_T == 0,
      na.rm = TRUE
    ),
    pct_NK_positive = 100 * mean(
      diagnostic_NK > 0,
      na.rm = TRUE
    ),
    pct_myeloid_positive = 100 * mean(
      diagnostic_Myeloid > 0,
      na.rm = TRUE
    ),
    pct_platelet_positive = 100 * mean(
      diagnostic_Platelet > 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    as.numeric(V3_state_cluster)
  )

write.csv(
  curation_qc_summary,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_curation_QC.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cluster centroid correlations
# ------------------------------------------------------------

centroid_features <- VariableFeatures(tcells_v3)

centroid_features <- intersect(
  centroid_features,
  rownames(rna_data)
)

centroid_features <- centroid_features[
  !grepl(
    "^(MT-|RPL|RPS|HBA|HBB|HBD)",
    centroid_features
  )
]

if (length(centroid_features) < 500) {
  stop(
    "Insufficient variable features for centroid analysis."
  )
}

cluster_levels <- sort(
  unique(
    as.integer(
      as.character(
        tcells_v3$V3_state_cluster
      )
    )
  )
)

cluster_factor <- factor(
  as.character(tcells_v3$V3_state_cluster),
  levels = as.character(cluster_levels)
)

cluster_design <- Matrix::sparse.model.matrix(
  ~ 0 + cluster_factor
)

colnames(cluster_design) <- levels(cluster_factor)

centroid_sums <- (
  rna_data[
    centroid_features,
    ,
    drop = FALSE
  ] %*%
    cluster_design
)

cluster_sizes <- as.numeric(
  table(cluster_factor)
)

centroid_means <- sweep(
  as.matrix(centroid_sums),
  2,
  cluster_sizes,
  "/"
)

centroid_correlation <- cor(
  centroid_means,
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  centroid_correlation,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_centroid_correlations.csv"
    )
  )
)

correlation_for_neighbors <- centroid_correlation
diag(correlation_for_neighbors) <- NA_real_

nearest_cluster_pairs <- as.data.frame(
  as.table(correlation_for_neighbors)
) %>%
  transmute(
    cluster = as.character(Var1),
    candidate_neighbor = as.character(Var2),
    spearman_correlation = Freq
  ) %>%
  filter(!is.na(spearman_correlation)) %>%
  group_by(cluster) %>%
  slice_max(
    spearman_correlation,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(as.numeric(cluster))

write.csv(
  nearest_cluster_pairs,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_nearest_cluster_pairs.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Focused review of cluster 10
# ------------------------------------------------------------

cluster10_neighbor <- nearest_cluster_pairs %>%
  filter(cluster == "10") %>%
  pull(candidate_neighbor)

if (length(cluster10_neighbor) != 1) {
  stop(
    "Unable to identify a unique neighbor for cluster 10."
  )
}

message(
  "Cluster 10 nearest centroid: cluster ",
  cluster10_neighbor
)

cluster10_markers <- FindMarkers(
  tcells_v3,
  ident.1 = "10",
  ident.2 = cluster10_neighbor,
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  max.cells.per.ident = 5000,
  random.seed = 123,
  verbose = FALSE
) %>%
  tibble::rownames_to_column("gene")

fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cluster10_markers)
)[1]

if (is.na(fc_column)) {
  stop("Fold-change column not found.")
}

cluster10_markers <- cluster10_markers %>%
  mutate(
    absolute_log2FC = abs(.data[[fc_column]]),
    absolute_detection_difference = abs(
      pct.1 - pct.2
    ),
    robust_difference =
      p_val_adj < 0.05 &
      absolute_log2FC >= 0.50 &
      absolute_detection_difference >= 0.10
  ) %>%
  arrange(
    desc(robust_difference),
    desc(absolute_log2FC)
  )

write.csv(
  cluster10_markers,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster10_vs_cluster",
      cluster10_neighbor,
      "_markers.csv"
    )
  ),
  row.names = FALSE
)

cluster10_marker_diagnostics <- tibble(
  cluster = "10",
  nearest_cluster = cluster10_neighbor,
  centroid_spearman_correlation =
    nearest_cluster_pairs %>%
    filter(cluster == "10") %>%
    pull(spearman_correlation),
  tested_genes = nrow(cluster10_markers),
  adjusted_significant_genes = sum(
    cluster10_markers$p_val_adj < 0.05,
    na.rm = TRUE
  ),
  robustly_different_genes = sum(
    cluster10_markers$robust_difference,
    na.rm = TRUE
  )
)

write.csv(
  cluster10_marker_diagnostics,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster10_curation_diagnostics.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Diagnostic heatmaps
# ------------------------------------------------------------

pdf(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_centroid_correlation_heatmap.pdf"
    )
  ),
  width = 9,
  height = 8
)

pheatmap(
  centroid_correlation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Similarity among V3 immune-state cluster centroids"
)

dev.off()

qc_heatmap_columns <- c(
  "median_nFeature_RNA",
  "median_percent_mt",
  "median_pan_T",
  "median_NK",
  "median_myeloid",
  "median_platelet",
  "median_erythroid"
)

qc_heatmap <- curation_qc_summary %>%
  select(
    V3_state_cluster,
    all_of(qc_heatmap_columns)
  ) %>%
  tibble::column_to_rownames(
    "V3_state_cluster"
  ) %>%
  as.matrix()

qc_heatmap_scaled <- scale(qc_heatmap)
qc_heatmap_scaled[!is.finite(qc_heatmap_scaled)] <- 0

pdf(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_curation_QC_heatmap.pdf"
    )
  ),
  width = 8,
  height = 7
)

pheatmap(
  qc_heatmap_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Quality and contamination profiles of V3 clusters"
)

dev.off()

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

message("\nQC profiles for clusters 9, 10, and 12:")

curation_focus <- curation_qc_summary %>%
  filter(
    V3_state_cluster %in% c("9", "10", "12")
  )

print(
  as.data.frame(curation_focus),
  row.names = FALSE
)

message("\nCluster 10 centroid review:")

print(
  as.data.frame(cluster10_marker_diagnostics),
  row.names = FALSE
)

message("\nTop robust differences for cluster 10:")

cluster10_robust_markers <- cluster10_markers %>%
  filter(robust_difference) %>%
  select(
    gene,
    all_of(fc_column),
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  head(20)

print(
  as.data.frame(cluster10_robust_markers),
  row.names = FALSE
)

message(
  "\nV3 cluster curation diagnostics completed."
)