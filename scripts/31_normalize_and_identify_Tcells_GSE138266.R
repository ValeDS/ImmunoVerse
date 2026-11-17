# ============================================================
# 31_normalize_and_identify_Tcells_GSE138266.R
# Global leukocyte representation and independent T-cell audit
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

INPUT_FILE <- paste0(
  "data/processed/external_validation/GSE138266/",
  "GSE138266_all_leukocytes_QC_counts.rds"
)

OUTPUT_DIR <- paste0(
  "data/processed/external_validation/GSE138266"
)
TABLE_DIR <- paste0(
  "results/tables/external_validation/GSE138266/",
  "Tcell_identification"
)
FIGURE_DIR <- paste0(
  "results/figures/external_validation/GSE138266"
)

OUTPUT_FILE <- file.path(
  OUTPUT_DIR,
  "GSE138266_all_leukocytes_normalized_clustered.rds"
)

REPRESENTATION_CHECKPOINT_FILE <- file.path(
  OUTPUT_DIR,
  "GSE138266_all_leukocytes_representation_checkpoint.rds"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input object not found: ",
    INPUT_FILE,
    ". Run script 30 first."
  )
}

representation_is_ready <- function(object) {
  inherits(object, "Seurat") &&
    ncol(object) == 75979L &&
    all(
      c(
        "pca",
        "harmony",
        "umap_harmony"
      ) %in% Reductions(object)
    ) &&
    "external_leukocyte_cluster" %in%
      colnames(object[[]])
}

reuse_in_memory <- exists(
  "external",
  envir = .GlobalEnv,
  inherits = FALSE
) &&
  representation_is_ready(
    get(
      "external",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  )

if (reuse_in_memory) {
  message(
    "Reusing the completed in-memory Harmony representation..."
  )

  external <- get(
    "external",
    envir = .GlobalEnv,
    inherits = FALSE
  )
} else if (
  file.exists(REPRESENTATION_CHECKPOINT_FILE)
) {
  message(
    "Loading the saved Harmony representation checkpoint..."
  )

  external <- readRDS(
    REPRESENTATION_CHECKPOINT_FILE
  )

  if (!representation_is_ready(external)) {
    stop(
      "The representation checkpoint is invalid or incomplete. ",
      "Move it before rerunning."
    )
  }
} else {
  message("Loading the QC-filtered GSE138266 object...")
  external <- readRDS(INPUT_FILE)

  if (
    !inherits(external, "Seurat") ||
      ncol(external) != 75979L ||
      n_distinct(external$sample) != 22L ||
      n_distinct(external$donor) != 12L
  ) {
    stop("The input object failed structural validation.")
  }

  DefaultAssay(external) <- "RNA"

  # ----------------------------------------------------------
  # Log normalization and unsupervised representation
  # ----------------------------------------------------------

  message("Normalizing expression values...")
  external <- NormalizeData(
    external,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  message("Selecting variable features...")
  external <- FindVariableFeatures(
    external,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  variable_features <- VariableFeatures(external)

  variable_features <- variable_features[
    !grepl(
      "^(MT-|RPL|RPS|HBA|HBB|HBD)",
      variable_features
    )
  ]

  if (length(variable_features) < 2000L) {
    stop("Too few non-technical variable features remain.")
  }

  VariableFeatures(external) <- variable_features

  message("Scaling variable features and running PCA...")
  external <- ScaleData(
    external,
    features = variable_features,
    verbose = FALSE
  )

  external <- RunPCA(
    external,
    features = variable_features,
    npcs = 30,
    verbose = FALSE
  )

  # Harmony is used only for broad-lineage clustering and
  # visualization. Raw RNA values remain unchanged.
  message("Running sample-aware Harmony representation...")
  external <- RunHarmony(
    object = external,
    group.by.vars = "sample",
    reduction.use = "pca",
    dims.use = 1:30,
    reduction.save = "harmony",
    project.dim = FALSE,
    verbose = FALSE
  )

  external <- RunUMAP(
    external,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap_harmony",
    reduction.key = "UMAPh_",
    seed.use = 123,
    verbose = FALSE
  )

  external <- FindNeighbors(
    external,
    reduction = "harmony",
    dims = 1:30,
    graph.name = c(
      "external_harmony_nn",
      "external_harmony_snn"
    ),
    verbose = FALSE
  )

  external <- FindClusters(
    external,
    graph.name = "external_harmony_snn",
    resolution = 0.4,
    cluster.name = "external_leukocyte_cluster",
    random.seed = 123,
    verbose = FALSE
  )

  external$external_leukocyte_cluster <- as.character(
    external$external_leukocyte_cluster
  )

  saveRDS(
    external,
    REPRESENTATION_CHECKPOINT_FILE,
    compress = FALSE
  )
}

DefaultAssay(external) <- "RNA"

# ------------------------------------------------------------
# Broad-lineage expression programs
# ------------------------------------------------------------

lineage_panels <- list(
  pan_T = c(
    "CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK"
  ),
  B_cell = c(
    "MS4A1", "CD79A", "CD79B", "CD37", "CD74", "CD22"
  ),
  myeloid = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CTSS",
    "S100A8", "S100A9", "FCN1"
  ),
  NK = c(
    "NKG7", "GNLY", "KLRD1", "NCR1", "FCER1G", "TYROBP"
  ),
  platelet = c(
    "PPBP", "PF4", "GNG11", "NRGN", "RGS18"
  ),
  erythroid = c(
    "HBA1", "HBA2", "HBB", "ALAS2", "AHSP"
  )
)

available_lineage_panels <- lapply(
  lineage_panels,
  intersect,
  y = rownames(external)
)

missing_panel_genes <- bind_rows(
  lapply(
    names(lineage_panels),
    function(program) {
      missing_genes <- setdiff(
        lineage_panels[[program]],
        available_lineage_panels[[program]]
      )

      if (length(missing_genes) == 0L) {
        return(
          data.frame(
            program = character(0),
            gene = character(0),
            stringsAsFactors = FALSE
          )
        )
      }

      data.frame(
        program = rep(
          program,
          length(missing_genes)
        ),
        gene = missing_genes,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  missing_panel_genes,
  file.path(
    TABLE_DIR,
    "GSE138266_missing_broad_lineage_marker_genes.csv"
  ),
  row.names = FALSE
)

if (any(lengths(available_lineage_panels) < 3L)) {
  stop(
    "One or more lineage panels contain fewer than three ",
    "available genes."
  )
}

normalized_data <- LayerData(
  external,
  assay = "RNA",
  layer = "data"
)

program_expression <- lapply(
  available_lineage_panels,
  function(genes) {
    Matrix::colMeans(
      normalized_data[
        genes,
        ,
        drop = FALSE
      ]
    )
  }
)

for (program in names(program_expression)) {
  external[[paste0(
    program,
    "_expression"
  )]] <- program_expression[[program]]
}

canonical_T_genes <- intersect(
  c("CD3D", "CD3E", "TRAC"),
  rownames(external)
)

canonical_T_detected <- Matrix::colSums(
  normalized_data[
    canonical_T_genes,
    ,
    drop = FALSE
  ] > 0
) > 0

external$canonical_T_detected <- canonical_T_detected

rm(normalized_data, program_expression)
gc(verbose = FALSE)

# ------------------------------------------------------------
# Cluster-level lineage audit
# ------------------------------------------------------------

cluster_audit <- external[[]] %>%
  group_by(external_leukocyte_cluster) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    n_donors = n_distinct(donor),
    n_compartments = n_distinct(compartment),
    pct_canonical_T_detected =
      100 * mean(canonical_T_detected),
    median_pan_T = median(pan_T_expression),
    median_B_cell = median(B_cell_expression),
    median_myeloid = median(myeloid_expression),
    median_NK = median(NK_expression),
    median_platelet = median(platelet_expression),
    median_erythroid = median(erythroid_expression),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    dominant_program = names(
      which.max(
        c(
          pan_T = median_pan_T,
          B_cell = median_B_cell,
          myeloid = median_myeloid,
          NK = median_NK,
          platelet = median_platelet,
          erythroid = median_erythroid
        )
      )
    ),
    provisional_lineage = case_when(
      dominant_program == "pan_T" ~ "T cell",
      dominant_program == "NK" &
        pct_canonical_T_detected >= 25 ~
        "T/NK candidate",
      dominant_program == "NK" ~ "NK cell",
      dominant_program == "B_cell" ~ "B cell",
      dominant_program == "myeloid" ~ "Myeloid",
      dominant_program == "platelet" ~ "Platelet",
      dominant_program == "erythroid" ~ "Erythroid",
      TRUE ~ "Unresolved"
    ),
    candidate_T_inclusion =
      provisional_lineage %in%
        c("T cell", "T/NK candidate")
  ) %>%
  ungroup() %>%
  arrange(
    as.integer(external_leukocyte_cluster)
  )

cluster_to_lineage <- setNames(
  cluster_audit$provisional_lineage,
  cluster_audit$external_leukocyte_cluster
)

external$provisional_lineage <- unname(
  cluster_to_lineage[
    external$external_leukocyte_cluster
  ]
)

# ------------------------------------------------------------
# Diagnostic cluster markers
# ------------------------------------------------------------

Idents(external) <- "external_leukocyte_cluster"

message("Finding broad-cluster diagnostic markers...")
cluster_markers <- FindAllMarkers(
  external,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  max.cells.per.ident = 2000,
  random.seed = 123,
  verbose = FALSE
)

fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cluster_markers)
)

if (length(fc_column) != 1L) {
  stop("Unable to identify the marker fold-change column.")
}

top_markers <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 10,
    with_ties = FALSE
  ) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  transmute(
    external_leukocyte_cluster = as.character(cluster),
    top_markers
  )

cluster_audit <- cluster_audit %>%
  left_join(
    top_markers,
    by = "external_leukocyte_cluster"
  )

write.csv(
  cluster_audit,
  file.path(
    TABLE_DIR,
    "GSE138266_broad_lineage_cluster_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  cluster_markers,
  file.path(
    TABLE_DIR,
    "GSE138266_broad_cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Diagnostic figure
# ------------------------------------------------------------

lineage_colors <- c(
  "T cell" = "#4C78A8",
  "T/NK candidate" = "#8E6C8A",
  "NK cell" = "#7A5195",
  "B cell" = "#59A14F",
  "Myeloid" = "#E45756",
  "Platelet" = "#F2CF5B",
  "Erythroid" = "#B279A2",
  "Unresolved" = "#9D9D9D"
)

panel_A <- DimPlot(
  external,
  reduction = "umap_harmony",
  group.by = "external_leukocyte_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.15
) +
  ggtitle("Unsupervised broad leukocyte clusters") +
  NoLegend()

panel_B <- DimPlot(
  external,
  reduction = "umap_harmony",
  group.by = "provisional_lineage",
  cols = lineage_colors,
  raster = TRUE,
  pt.size = 0.15
) +
  ggtitle("Marker-supported provisional lineages")

panel_C <- DimPlot(
  external,
  reduction = "umap_harmony",
  group.by = "compartment",
  cols = c(
    "PBMC" = "#4C78A8",
    "CSF" = "#E45756"
  ),
  raster = TRUE,
  pt.size = 0.15
) +
  ggtitle("Compartment")

diagnostic_markers <- unique(
  unlist(
    lapply(
      available_lineage_panels,
      head,
      n = 4L
    )
  )
)

panel_D <- DotPlot(
  external,
  features = diagnostic_markers,
  group.by = "external_leukocyte_cluster",
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle("Broad-lineage diagnostic markers") +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
  )

diagnostic_figure <- (
  panel_A + panel_B + panel_C
) /
  panel_D +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_broad_lineage_and_Tcell_audit.png"
  ),
  diagnostic_figure,
  width = 17,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE138266_broad_lineage_and_Tcell_audit.pdf"
  ),
  diagnostic_figure,
  width = 17,
  height = 11,
  device = grDevices::pdf,
  bg = "white"
)

# ------------------------------------------------------------
# Save checkpoint
# ------------------------------------------------------------

saveRDS(
  external,
  OUTPUT_FILE,
  compress = FALSE
)

candidate_T_clusters <- cluster_audit %>%
  filter(candidate_T_inclusion) %>%
  select(
    external_leukocyte_cluster,
    provisional_lineage,
    n_cells,
    n_samples,
    n_donors,
    pct_canonical_T_detected,
    median_pan_T,
    median_NK,
    top_markers
  )

message("")
message("GSE138266 broad-lineage analysis completed.")
message("")
message("Provisional lineage summary:")
print(
  as.data.frame(
    external[[]] %>%
      count(
        provisional_lineage,
        compartment,
        group,
        name = "n_cells"
      ) %>%
      arrange(
        provisional_lineage,
        compartment,
        group
      )
  ),
  row.names = FALSE
)
message("")
message("Candidate T-cell clusters requiring review:")
print(
  as.data.frame(candidate_T_clusters),
  row.names = FALSE
)
message("")
message("Saved checkpoint: ", OUTPUT_FILE)
