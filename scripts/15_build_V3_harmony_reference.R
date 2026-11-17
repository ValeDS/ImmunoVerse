# ============================================================
# 15_build_V3_harmony_reference.R
# Chemistry-homogeneous V3 T-cell discovery reference
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

REFERENCE_CHEMISTRY <- "V3"
MAX_CELLS_PER_SAMPLE <- 2000
N_VARIABLE_FEATURES_REFERENCE <- 3000
N_PCS_REFERENCE <- 30
INTEGRATED_RESOLUTION <- 0.6

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

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
# Load full T-cell object
# ------------------------------------------------------------

in_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    GEO_ACCESSION,
    "_tcells_raw.rds"
  )
)

if (!file.exists(in_rds)) {
  stop(
    "Full T-cell object not found: ",
    in_rds
  )
}

message("Loading full T-cell object...")

tcells_full <- readRDS(in_rds)
DefaultAssay(tcells_full) <- "RNA"

required_metadata <- c(
  "sample",
  "donor",
  "group",
  "cohort",
  "V_10X",
  "exp_name"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(tcells_full@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

if (
  anyNA(tcells_full$V_10X) ||
  anyNA(tcells_full$exp_name)
) {
  stop(
    "V_10X or exp_name contains missing values."
  )
}

# Keep counts, normalized RNA data, and metadata while removing
# previous scaled data, reductions, and graphs.
tcells_full <- DietSeurat(
  tcells_full,
  assays = "RNA",
  layers = c(
    "counts",
    "data"
  ),
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)

gc()

# ------------------------------------------------------------
# Restrict discovery reference to V3 chemistry
# ------------------------------------------------------------

message(
  "Restricting discovery reference to ",
  REFERENCE_CHEMISTRY,
  " chemistry..."
)

v3_cells <- colnames(tcells_full)[
  tcells_full$V_10X ==
    REFERENCE_CHEMISTRY
]

if (length(v3_cells) == 0) {
  stop(
    "No cells found for chemistry ",
    REFERENCE_CHEMISTRY
  )
}

tcells_v3_full <- subset(
  tcells_full,
  cells = v3_cells
)

rm(
  tcells_full,
  v3_cells
)

gc()

# ------------------------------------------------------------
# Balanced downsampling by biological sample
# ------------------------------------------------------------

cell_metadata <- tcells_v3_full@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  )

selected_cells <- cell_metadata %>%
  group_by(sample) %>%
  group_modify(
    ~ {
      n_to_sample <- min(
        nrow(.x),
        MAX_CELLS_PER_SAMPLE
      )
      
      slice_sample(
        .x,
        n = n_to_sample
      )
    }
  ) %>%
  ungroup() %>%
  pull(cell_id)

sampling_summary <- cell_metadata %>%
  count(
    sample,
    donor,
    group,
    cohort,
    exp_name,
    name = "available_V3_tcells"
  ) %>%
  mutate(
    selected_V3_tcells = pmin(
      available_V3_tcells,
      MAX_CELLS_PER_SAMPLE
    ),
    sampling_fraction =
      selected_V3_tcells /
      available_V3_tcells
  )

write.csv(
  sampling_summary,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_sampling_summary.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    cell_id = selected_cells
  ),
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_selected_cells.csv"
    )
  ),
  row.names = FALSE
)

message(
  "Selected ",
  length(selected_cells),
  " V3 cells from ",
  dplyr::n_distinct(
    cell_metadata$sample
  ),
  " biological samples."
)

tcells_v3_reference <- subset(
  tcells_v3_full,
  cells = selected_cells
)

rm(
  tcells_v3_full,
  cell_metadata,
  selected_cells
)

gc()

# ------------------------------------------------------------
# Recompute RNA representation
# ------------------------------------------------------------

DefaultAssay(tcells_v3_reference) <- "RNA"

tcells_v3_reference <- NormalizeData(
  tcells_v3_reference,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

tcells_v3_reference <- FindVariableFeatures(
  tcells_v3_reference,
  selection.method = "vst",
  nfeatures =
    N_VARIABLE_FEATURES_REFERENCE,
  verbose = TRUE
)

tcells_v3_reference <- ScaleData(
  tcells_v3_reference,
  features = VariableFeatures(
    tcells_v3_reference
  ),
  verbose = TRUE
)

tcells_v3_reference <- RunPCA(
  tcells_v3_reference,
  features = VariableFeatures(
    tcells_v3_reference
  ),
  npcs = 50,
  reduction.name = "pca",
  verbose = TRUE
)

# ------------------------------------------------------------
# Uncorrected UMAP
# ------------------------------------------------------------

tcells_v3_reference <- RunUMAP(
  tcells_v3_reference,
  reduction = "pca",
  dims = 1:N_PCS_REFERENCE,
  reduction.name =
    "umap_unintegrated",
  reduction.key = "UMAPraw_",
  seed.use = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Harmony correction across V3 experimental batches
# ------------------------------------------------------------

message(
  "Running Harmony on V3 exp_name batches..."
)

tcells_v3_reference <- RunHarmony(
  object = tcells_v3_reference,
  group.by.vars = "exp_name",
  reduction.use = "pca",
  dims.use = 1:N_PCS_REFERENCE,
  reduction.save = "harmony",
  project.dim = FALSE,
  plot_convergence = TRUE,
  verbose = TRUE
)

if (
  !"harmony" %in%
    Reductions(tcells_v3_reference)
) {
  stop(
    "Harmony reduction was not created."
  )
}

# ------------------------------------------------------------
# Harmony-based neighborhood, clustering, and UMAP
# ------------------------------------------------------------

tcells_v3_reference <- FindNeighbors(
  tcells_v3_reference,
  reduction = "harmony",
  dims = 1:N_PCS_REFERENCE,
  k.param = 20,
  graph.name = c(
    "harmony_nn",
    "harmony_snn"
  ),
  verbose = TRUE
)

tcells_v3_reference <- FindClusters(
  tcells_v3_reference,
  graph.name = "harmony_snn",
  resolution =
    INTEGRATED_RESOLUTION,
  algorithm = 1,
  random.seed = 123,
  cluster.name =
    "V3_state_cluster",
  verbose = TRUE
)

tcells_v3_reference$V3_state_cluster <-
  factor(
    tcells_v3_reference$V3_state_cluster
  )

Idents(tcells_v3_reference) <-
  "V3_state_cluster"

tcells_v3_reference <- RunUMAP(
  tcells_v3_reference,
  reduction = "harmony",
  dims = 1:N_PCS_REFERENCE,
  reduction.name = "umap_harmony",
  reduction.key = "UMAPh_",
  seed.use = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Integration diagnostic figures
# ------------------------------------------------------------

plot_theme <- theme(
  plot.title = element_text(
    face = "bold",
    size = 12
  ),
  legend.title = element_text(
    face = "bold"
  )
)

p_before_cohort <- DimPlot(
  tcells_v3_reference,
  reduction = "umap_unintegrated",
  group.by = "cohort",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Before Harmony: cohort"
  ) +
  plot_theme

p_after_cohort <- DimPlot(
  tcells_v3_reference,
  reduction = "umap_harmony",
  group.by = "cohort",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "After Harmony: cohort"
  ) +
  plot_theme

p_before_batch <- DimPlot(
  tcells_v3_reference,
  reduction = "umap_unintegrated",
  group.by = "exp_name",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Before Harmony: experimental batch"
  ) +
  plot_theme

p_after_batch <- DimPlot(
  tcells_v3_reference,
  reduction = "umap_harmony",
  group.by = "exp_name",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "After Harmony: experimental batch"
  ) +
  plot_theme

p_integration_diagnostic <-
  (
    p_before_cohort |
      p_after_cohort
  ) /
  (
    p_before_batch |
      p_after_batch
  ) +
  plot_annotation(
    title = paste0(
      "GSE144744 V3 T-cell reference: ",
      "integration diagnostics"
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_before_after_integration.pdf"
    )
  ),
  p_integration_diagnostic,
  width = 16,
  height = 11
)

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_before_after_integration.png"
    )
  ),
  p_integration_diagnostic,
  width = 16,
  height = 11,
  dpi = 300
)

p_clusters <- DimPlot(
  tcells_v3_reference,
  reduction = "umap_harmony",
  group.by = "V3_state_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle(
    "V3 Harmony-integrated T-cell clusters"
  ) +
  NoLegend() +
  plot_theme

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_clusters.pdf"
    )
  ),
  p_clusters,
  width = 10,
  height = 8
)

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_clusters.png"
    )
  ),
  p_clusters,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Initial portability audit
# ------------------------------------------------------------

v3_meta <- tcells_v3_reference@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  )

cluster_counts <- v3_meta %>%
  count(
    V3_state_cluster,
    name = "n_cells"
  ) %>%
  mutate(
    percent_cells =
      100 * n_cells /
      sum(n_cells)
  )

cluster_cohort <- v3_meta %>%
  count(
    V3_state_cluster,
    cohort,
    name = "n_cells"
  ) %>%
  group_by(
    V3_state_cluster
  ) %>%
  mutate(
    percent_within_cluster =
      100 * n_cells /
      sum(n_cells)
  ) %>%
  ungroup()

cluster_donor <- v3_meta %>%
  count(
    V3_state_cluster,
    donor,
    name = "n_cells"
  ) %>%
  group_by(
    V3_state_cluster
  ) %>%
  mutate(
    percent_within_cluster =
      100 * n_cells /
      sum(n_cells)
  ) %>%
  ungroup()

cluster_audit_V3 <- cluster_counts %>%
  left_join(
    cluster_cohort %>%
      group_by(
        V3_state_cluster
      ) %>%
      summarise(
        n_cohorts = n_distinct(
          cohort[n_cells > 0]
        ),
        dominant_cohort =
          cohort[
            which.max(n_cells)
          ],
        maximum_cohort_contribution_percent =
          max(
            percent_within_cluster
          ),
        .groups = "drop"
      ),
    by = "V3_state_cluster"
  ) %>%
  left_join(
    cluster_donor %>%
      group_by(
        V3_state_cluster
      ) %>%
      summarise(
        n_donors = n_distinct(
          donor[n_cells > 0]
        ),
        dominant_donor =
          donor[
            which.max(n_cells)
          ],
        maximum_donor_contribution_percent =
          max(
            percent_within_cluster
          ),
        .groups = "drop"
      ),
    by = "V3_state_cluster"
  ) %>%
  arrange(
    V3_state_cluster
  )

write.csv(
  cluster_counts,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_counts.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cluster_cohort,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_by_cohort.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cluster_donor,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_cluster_by_donor.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cluster_audit_V3,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_initial_cluster_audit.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Save V3 discovery reference
# ------------------------------------------------------------

out_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    ".rds"
  )
)

saveRDS(
  tcells_v3_reference,
  out_rds,
  compress = TRUE
)

message(
  "Saved V3 Harmony reference: ",
  out_rds
)

message(
  "V3 reference cells: ",
  ncol(tcells_v3_reference)
)

message(
  "V3 integrated clusters: ",
  nrow(cluster_counts)
)

print(
  tibble::as_tibble(
    cluster_audit_V3
  ),
  n = Inf,
  width = Inf
)

message(
  "V3 Harmony reference completed."
)
