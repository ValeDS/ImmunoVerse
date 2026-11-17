# ============================================================
# 13_build_harmony_integrated_tcell_reference.R
# Balanced T-cell reference construction with Harmony correction
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

MAX_CELLS_PER_SAMPLE <- 2000
N_VARIABLE_FEATURES_REFERENCE <- 3000
N_PCS_REFERENCE <- 30
INTEGRATED_RESOLUTION <- 0.6

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_Harmony_reference"
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

# Remove previous dimensional reductions and scaled data before
# subsetting. Counts, normalized data, and metadata are retained.
message("Creating a lightweight starting object...")

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
# Balanced downsampling by biological sample
# ------------------------------------------------------------

cell_metadata <- tcells_full@meta.data %>%
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
    V_10X,
    exp_name,
    name = "available_tcells"
  ) %>%
  mutate(
    selected_tcells = pmin(
      available_tcells,
      MAX_CELLS_PER_SAMPLE
    ),
    sampling_fraction =
      selected_tcells /
      available_tcells
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
  " cells from ",
  dplyr::n_distinct(
    cell_metadata$sample
  ),
  " biological samples."
)

tcells_reference <- subset(
  tcells_full,
  cells = selected_cells
)

rm(
  tcells_full,
  cell_metadata,
  selected_cells
)

gc()

# ------------------------------------------------------------
# Recompute RNA representation
# ------------------------------------------------------------

DefaultAssay(tcells_reference) <- "RNA"

message("Normalizing RNA data...")

tcells_reference <- NormalizeData(
  tcells_reference,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

message("Selecting variable features...")

tcells_reference <- FindVariableFeatures(
  tcells_reference,
  selection.method = "vst",
  nfeatures =
    N_VARIABLE_FEATURES_REFERENCE,
  verbose = TRUE
)

message("Scaling variable features...")

tcells_reference <- ScaleData(
  tcells_reference,
  features = VariableFeatures(
    tcells_reference
  ),
  verbose = TRUE
)

message("Running PCA...")

tcells_reference <- RunPCA(
  tcells_reference,
  features = VariableFeatures(
    tcells_reference
  ),
  npcs = 50,
  reduction.name = "pca",
  verbose = TRUE
)

# ------------------------------------------------------------
# Uncorrected UMAP for before/after diagnostics
# ------------------------------------------------------------

message("Computing uncorrected UMAP...")

tcells_reference <- RunUMAP(
  tcells_reference,
  reduction = "pca",
  dims = 1:N_PCS_REFERENCE,
  reduction.name = "umap_unintegrated",
  reduction.key = "UMAPraw_",
  seed.use = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Harmony correction of technical covariates only
# ------------------------------------------------------------

message(
  "Running Harmony using V_10X and exp_name..."
)

tcells_reference <- RunHarmony(
  object = tcells_reference,
  group.by.vars = c(
    "V_10X",
    "exp_name"
  ),
  reduction.use = "pca",
  dims.use = 1:N_PCS_REFERENCE,
  reduction.save = "harmony",
  project.dim = FALSE,
  plot_convergence = TRUE,
  verbose = TRUE
)

if (!"harmony" %in% Reductions(tcells_reference)) {
  stop(
    "Harmony reduction was not created."
  )
}

# ------------------------------------------------------------
# Harmony-based neighborhood, clustering, and UMAP
# ------------------------------------------------------------

message("Computing Harmony neighborhood graph...")

tcells_reference <- FindNeighbors(
  tcells_reference,
  reduction = "harmony",
  dims = 1:N_PCS_REFERENCE,
  k.param = 20,
  graph.name = c(
    "harmony_nn",
    "harmony_snn"
  ),
  verbose = TRUE
)

message("Clustering the integrated reference...")

tcells_reference <- FindClusters(
  tcells_reference,
  graph.name = "harmony_snn",
  resolution = INTEGRATED_RESOLUTION,
  algorithm = 1,
  random.seed = 123,
  cluster.name =
    "integrated_state_cluster",
  verbose = TRUE
)

tcells_reference$integrated_state_cluster <-
  factor(
    tcells_reference$integrated_state_cluster
  )

Idents(tcells_reference) <-
  "integrated_state_cluster"

message("Computing Harmony UMAP...")

tcells_reference <- RunUMAP(
  tcells_reference,
  reduction = "harmony",
  dims = 1:N_PCS_REFERENCE,
  reduction.name = "umap_harmony",
  reduction.key = "UMAPh_",
  seed.use = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Before/after integration plots
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
  tcells_reference,
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
  tcells_reference,
  reduction = "umap_harmony",
  group.by = "cohort",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "After Harmony: cohort"
  ) +
  plot_theme

p_before_chemistry <- DimPlot(
  tcells_reference,
  reduction = "umap_unintegrated",
  group.by = "V_10X",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Before Harmony: 10x chemistry"
  ) +
  plot_theme

p_after_chemistry <- DimPlot(
  tcells_reference,
  reduction = "umap_harmony",
  group.by = "V_10X",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle(
    "After Harmony: 10x chemistry"
  ) +
  plot_theme

p_integration_diagnostic <-
  (
    p_before_cohort |
      p_after_cohort
  ) /
  (
    p_before_chemistry |
      p_after_chemistry
  ) +
  plot_annotation(
    title = paste0(
      "GSE144744 T-cell reference: ",
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
  width = 14,
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
  width = 14,
  height = 11,
  dpi = 300
)

p_clusters <- DimPlot(
  tcells_reference,
  reduction = "umap_harmony",
  group.by =
    "integrated_state_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Harmony-integrated T-cell clusters"
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
# Cluster representation audit
# ------------------------------------------------------------

integrated_meta <- tcells_reference@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  )

cluster_counts <- integrated_meta %>%
  count(
    integrated_state_cluster,
    name = "n_cells"
  ) %>%
  mutate(
    percent_cells =
      100 * n_cells /
      sum(n_cells)
  )

cluster_cohort <- integrated_meta %>%
  count(
    integrated_state_cluster,
    cohort,
    name = "n_cells"
  ) %>%
  group_by(
    integrated_state_cluster
  ) %>%
  mutate(
    percent_within_cluster =
      100 * n_cells /
      sum(n_cells)
  ) %>%
  ungroup()

cluster_donor <- integrated_meta %>%
  count(
    integrated_state_cluster,
    donor,
    name = "n_cells"
  ) %>%
  group_by(
    integrated_state_cluster
  ) %>%
  mutate(
    percent_within_cluster =
      100 * n_cells /
      sum(n_cells)
  ) %>%
  ungroup()

cluster_audit <- cluster_counts %>%
  left_join(
    cluster_cohort %>%
      group_by(
        integrated_state_cluster
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
    by = "integrated_state_cluster"
  ) %>%
  left_join(
    cluster_donor %>%
      group_by(
        integrated_state_cluster
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
    by = "integrated_state_cluster"
  ) %>%
  arrange(
    integrated_state_cluster
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
  cluster_audit,
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
# Save integrated reference
# ------------------------------------------------------------

out_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    ".rds"
  )
)

saveRDS(
  tcells_reference,
  out_rds,
  compress = TRUE
)

message(
  "Saved Harmony-integrated reference: ",
  out_rds
)

message(
  "Integrated clusters: ",
  nrow(cluster_counts)
)

print(
  tibble::as_tibble(
    cluster_audit
  ),
  n = Inf,
  width = Inf
)

message(
  "Harmony-integrated T-cell reference completed."
)
