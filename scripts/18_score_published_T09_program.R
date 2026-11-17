# ============================================================
# 18_score_published_T09_program.R
# Relative enrichment of published T09-associated programs
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(OUT_PREFIX, ".rds")
)

message("Loading V3 Harmony reference...")

tcells_v3 <- readRDS(in_rds)
DefaultAssay(tcells_v3) <- "RNA"

curated_clusters <- c(
  "0", "1", "2", "3", "4", "5",
  "6", "7", "8", "9", "11"
)

tcells_curated_candidates <- subset(
  tcells_v3,
  subset = V3_state_cluster %in% curated_clusters
)

rna_data <- LayerData(
  tcells_curated_candidates[["RNA"]],
  layer = "data"
)

# Marker panels reported in the original publication.
# These are focused panels, not the complete supplementary signatures.

published_programs <- list(
  T09_marker_panel = c(
    "KLRB1", "LTB", "AQP3", "ICOS",
    "TNFRSF25", "MAF", "RORA", "NINJ1",
    "LIMS1", "ITGA4", "PDE4D", "ABCB1",
    "CCR1", "LGALS1"
  ),
  CNS_h_reported_panel = c(
    "ITGA4", "AQP3", "EZR", "ADGRE5",
    "RHOB", "RGS1", "PIK3CG", "KRAS",
    "SOCS3", "TOB1"
  )
)

published_programs <- lapply(
  published_programs,
  function(genes) {
    intersect(genes, rownames(rna_data))
  }
)

print(
  lapply(
    published_programs,
    function(x) {
      list(
        n_genes = length(x),
        genes = x
      )
    }
  )
)

all_program_genes <- unique(
  unlist(published_programs)
)

cluster_factor <- factor(
  as.character(
    tcells_curated_candidates$V3_state_cluster
  ),
  levels = curated_clusters
)

cluster_design <- Matrix::sparse.model.matrix(
  ~ 0 + cluster_factor
)

colnames(cluster_design) <- curated_clusters

cluster_sums <- (
  rna_data[
    all_program_genes,
    ,
    drop = FALSE
  ] %*%
    cluster_design
)

cluster_sizes <- as.numeric(
  table(cluster_factor)
)

cluster_means <- sweep(
  as.matrix(cluster_sums),
  2,
  cluster_sizes,
  "/"
)

# Standardize each gene across clusters so that highly expressed
# genes do not dominate the combined program scores.

cluster_means_scaled <- t(
  scale(
    t(cluster_means)
  )
)

cluster_means_scaled[
  !is.finite(cluster_means_scaled)
] <- 0

calculate_program_score <- function(genes) {
  
  genes_present <- intersect(
    genes,
    rownames(cluster_means_scaled)
  )
  
  colMeans(
    cluster_means_scaled[
      genes_present,
      ,
      drop = FALSE
    ]
  )
}

published_program_scores <- tibble(
  V3_state_cluster = curated_clusters,
  T09_marker_score = calculate_program_score(
    published_programs$T09_marker_panel
  ),
  CNS_h_reported_score = calculate_program_score(
    published_programs$CNS_h_reported_panel
  )
)

# Add the fraction derived from the original T09 cluster.

source_T09_fraction <- tcells_curated_candidates@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    )
  ) %>%
  count(
    V3_state_cluster,
    cluster_names,
    name = "n_cells"
  ) %>%
  group_by(V3_state_cluster) %>%
  mutate(
    source_percent =
      100 * n_cells / sum(n_cells)
  ) %>%
  filter(cluster_names == "T09") %>%
  select(
    V3_state_cluster,
    original_T09_percent = source_percent
  ) %>%
  ungroup()

published_program_scores <- published_program_scores %>%
  left_join(
    source_T09_fraction,
    by = "V3_state_cluster"
  ) %>%
  mutate(
    original_T09_percent = replace_na(
      original_T09_percent,
      0
    )
  ) %>%
  arrange(
    desc(T09_marker_score)
  )

write.csv(
  published_program_scores,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_published_T09_program_scores.csv"
    )
  ),
  row.names = FALSE
)

# Gene-level heatmap

annotation_col <- data.frame(
  V3_state_cluster = colnames(
    cluster_means_scaled
  )
)

rownames(annotation_col) <- colnames(
  cluster_means_scaled
)

pdf(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_published_T09_program_heatmap.pdf"
    )
  ),
  width = 10,
  height = 8
)

pheatmap(
  cluster_means_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  annotation_col = annotation_col,
  main = "Relative expression of published T09-associated programs"
)

dev.off()

message("\nPublished T09-associated program scores:")

print(
  as.data.frame(published_program_scores),
  row.names = FALSE
)

message(
  "\nPublished T09 program analysis completed."
)