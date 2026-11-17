# ============================================================
# 06_tcell_annotation.R
# T-cell cluster annotation and immune-state scoring
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

message("Starting T-cell annotation...")

# ------------------------------------------------------------
# 1. Load T-cell object
# ------------------------------------------------------------

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(GEO_ACCESSION, "_tcells_safe_clustered.rds")
)

if (!file.exists(in_rds)) {
  stop("T-cell clustered object not found: ", in_rds)
}

tcells <- readRDS(in_rds)
DefaultAssay(tcells) <- "RNA"

message("Loaded T-cell object:")
print(tcells)

# ------------------------------------------------------------
# 2. Define T-cell functional signatures
# ------------------------------------------------------------

signatures <- list(
  naive_memory = c("CCR7", "IL7R", "SELL", "LEF1", "TCF7", "LTB"),
  cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "CCL5"),
  activation = c("CD69", "IL2RA", "CD40LG", "HLA-DRA", "HLA-DRB1"),
  exhaustion = c("PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2"),
  treg = c("FOXP3", "IL2RA", "IKZF2", "CTLA4"),
  cns_homing = c("KLRB1", "LTB", "CCR6", "CXCR3", "CXCR4", "ITGA4"),
  cd4_like = c("CD4", "IL7R", "CCR7"),
  cd8_like = c("CD8A", "CD8B", "NKG7", "GZMK")
)

# Keep only genes present
signatures <- lapply(signatures, function(x) x[x %in% rownames(tcells)])
signatures <- signatures[lengths(signatures) > 0]

# ------------------------------------------------------------
# 3. Add module scores
# ------------------------------------------------------------

for (sig_name in names(signatures)) {
  message("Scoring signature: ", sig_name)
  tcells <- AddModuleScore(
    object = tcells,
    features = list(signatures[[sig_name]]),
    name = paste0(sig_name, "_score")
  )
}

# AddModuleScore creates columns like naive_memory_score1
score_cols <- grep("_score1$", colnames(tcells@meta.data), value = TRUE)

# ------------------------------------------------------------
# 4. Cluster-level mean signature scores
# ------------------------------------------------------------

cluster_scores <- tcells@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    across(all_of(score_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cluster_scores,
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_cluster_signature_scores.csv"))
)

# ------------------------------------------------------------
# 5. Cluster-level marker average expression
# ------------------------------------------------------------

marker_genes <- unique(unlist(signatures))
marker_genes <- marker_genes[marker_genes %in% rownames(tcells)]

avg_exp <- AverageExpression(
  tcells,
  features = marker_genes,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"
)$RNA

avg_exp_df <- as.data.frame(t(avg_exp)) %>%
  rownames_to_column("seurat_clusters")

write_csv(
  avg_exp_df,
  file.path(TAB_DIR, paste0(GEO_ACCESSION, "_Tcells_cluster_marker_average_expression.csv"))
)

# ------------------------------------------------------------
# 6. DotPlot for annotation
# ------------------------------------------------------------

p_dot <- DotPlot(
  tcells,
  features = marker_genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("T-cell annotation markers")

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_annotation_dotplot.pdf")),
  p_dot,
  width = 18,
  height = 8
)

# ------------------------------------------------------------
# 7. Feature plots for key genes
# ------------------------------------------------------------

key_genes <- c(
  "CCR7", "IL7R", "SELL", "LEF1", "TCF7",
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMK",
  "KLRB1", "LTB", "CXCR3", "CXCR4", "ITGA4",
  "CD69", "IL2RA", "FOXP3", "PDCD1", "TIGIT"
)

key_genes <- key_genes[key_genes %in% rownames(tcells)]

p_feat <- FeaturePlot(
  tcells,
  features = key_genes,
  reduction = "umap",
  ncol = 4,
  raster = TRUE
)

ggsave(
  file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_key_gene_featureplots.pdf")),
  p_feat,
  width = 16,
  height = ceiling(length(key_genes) / 4) * 4
)

# ------------------------------------------------------------
# 8. UMAP signature scores
# ------------------------------------------------------------

for (score_col in score_cols) {
  p_score <- FeaturePlot(
    tcells,
    features = score_col,
    reduction = "umap",
    raster = TRUE
  ) +
    ggtitle(score_col)

  ggsave(
    file.path(FIG_DIR, paste0(GEO_ACCESSION, "_Tcells_", score_col, "_UMAP.pdf")),
    p_score,
    width = 7,
    height = 6
  )
}

# ------------------------------------------------------------
# 9. Save annotated object
# ------------------------------------------------------------

out_rds <- file.path(
  PROCESSED_DIR,
  paste0(GEO_ACCESSION, "_tcells_annotated.rds")
)

saveRDS(tcells, out_rds)

message("Saved annotated T-cell object: ", out_rds)
message("T-cell annotation completed.")
