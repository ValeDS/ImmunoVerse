# ============================================================
# 19_build_curated_V3_immune_state_atlas.R
# Construction of the curated V3 Immune State Atlas
# ============================================================

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

OUT_PREFIX <- paste0(
  GEO_ACCESSION,
  "_Tcells_V3_Harmony_reference"
)

dir.create(
  PROCESSED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

reference_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    ".rds"
  )
)

marker_file <- file.path(
  TAB_DIR,
  paste0(
    OUT_PREFIX,
    "_marker_summary.csv"
  )
)

source_crosswalk_file <- file.path(
  TAB_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_state_source_crosswalk.csv"
  )
)

program_score_file <- file.path(
  TAB_DIR,
  paste0(
    OUT_PREFIX,
    "_published_T09_program_scores.csv"
  )
)

required_files <- c(
  reference_rds,
  marker_file,
  source_crosswalk_file,
  program_score_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# Load V3 Harmony reference
# ------------------------------------------------------------

message("Loading V3 Harmony reference...")

tcells_v3 <- readRDS(reference_rds)

DefaultAssay(tcells_v3) <- "RNA"

if (
  !"V3_state_cluster" %in%
  colnames(tcells_v3@meta.data)
) {
  stop(
    "Metadata column V3_state_cluster is missing."
  )
}

if (
  !"umap_harmony" %in%
  Reductions(tcells_v3)
) {
  stop(
    "Reduction umap_harmony is missing."
  )
}

tcells_v3$V3_state_cluster <- as.character(
  tcells_v3$V3_state_cluster
)

# ------------------------------------------------------------
# Curated state definitions
# ------------------------------------------------------------

state_definitions <- tibble::tribble(
  ~V3_state_cluster,
  ~immune_state_id,
  ~immune_state_name,
  ~state_family,
  ~confidence,
  ~interpretation,
  
  "0",
  "IVS-01",
  "CCR7+ naive/central-memory CD4 state",
  "CD4 naive/memory",
  "High",
  paste(
    "Quiescent CCR7-associated CD4 state",
    "derived predominantly from source cluster T02."
  ),
  
  "1",
  "IVS-02",
  "ANKRD55+ naive/central-memory CD4 state",
  "CD4 naive/memory",
  "High",
  paste(
    "Second T02-derived CD4 state characterized",
    "by ANKRD55, LRRN3, IGF1R and related",
    "quiescent-memory markers."
  ),
  
  "2",
  "IVS-03",
  "MAF+ LIMS1+ memory CD4 state",
  "CD4 memory",
  "Medium",
  paste(
    "Memory CD4 state containing T06- and",
    "T09-derived cells and expressing part of",
    "the published T09 marker program."
  ),
  
  "3",
  "IVS-04",
  "Activated/migratory memory CD4 state",
  "Activated CD4 memory",
  "High",
  paste(
    "FUT7+, TNFRSF4+ and TNFRSF18+ memory",
    "CD4 state showing the strongest combined",
    "evidence of T09 association and",
    "CNS-homing activity."
  ),
  
  "4",
  "IVS-05",
  "CX3CR1+ cytotoxic effector CD8/TEMRA state",
  "Cytotoxic CD8",
  "High",
  paste(
    "Highly cytotoxic effector state defined",
    "by CX3CR1, FGFBP2, GZMH, GZMB and ADGRG1."
  ),
  
  "5",
  "IVS-06",
  "Naive/central-memory CD8 state",
  "CD8 naive/memory",
  "High",
  paste(
    "CD8 state with low cytotoxic activity",
    "and a naive or central-memory",
    "transcriptional profile."
  ),
  
  "6",
  "IVS-07",
  "GZMK+ effector-memory CD8 state",
  "CD8 effector memory",
  "High",
  paste(
    "Effector-memory CD8 state characterized",
    "by GZMK, CCL5, EOMES and CRTAM."
  ),
  
  "7",
  "IVS-08",
  "KLRB1+ SLC4A10+ MAIT-like T state",
  "Innate-like T",
  "High",
  paste(
    "Innate-like T-cell state defined by",
    "KLRB1, SLC4A10, ZBTB16, IL18RAP and NCR3."
  ),
  
  "8",
  "IVS-09",
  "FOXP3+ regulatory T-cell state",
  "Regulatory T",
  "High",
  paste(
    "Canonical regulatory state defined by",
    "FOXP3, IL2RA, CTLA4, TIGIT and RTKN2."
  ),
  
  "9",
  "IVS-10",
  "KLRC1+ NK-like cytotoxic T-cell state",
  "NK-like T",
  "Medium",
  paste(
    "T-cell-associated population with a",
    "strong NK-like program including KLRC1,",
    "TYROBP, FCER1G and PTGDS."
  ),
  
  "11",
  "IVS-11",
  "Interferon-responsive T-cell state",
  "Interferon-responsive T",
  "High",
  paste(
    "Cross-lineage functional state defined by",
    "IFIT1, IFIT3, IFI44L, RSAD2, OAS1,",
    "MX1 and ISG15."
  )
)

excluded_definitions <- tibble::tribble(
  ~V3_state_cluster,
  ~exclusion_reason,
  
  "10",
  "Low-quality and low-complexity T-cell cluster",
  
  "12",
  "Platelet-myeloid contaminated cluster"
)

# Save the explicit curation definitions.

write.csv(
  state_definitions,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_curated_state_definitions.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  excluded_definitions,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_excluded_state_definitions.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Check cluster coverage
# ------------------------------------------------------------

observed_clusters <- sort(
  unique(
    tcells_v3$V3_state_cluster
  )
)

defined_clusters <- c(
  state_definitions$V3_state_cluster,
  excluded_definitions$V3_state_cluster
)

unmapped_clusters <- setdiff(
  observed_clusters,
  defined_clusters
)

if (length(unmapped_clusters) > 0) {
  stop(
    "The following V3 clusters have no curation decision: ",
    paste(
      unmapped_clusters,
      collapse = ", "
    )
  )
}

duplicated_definitions <- defined_clusters[
  duplicated(defined_clusters)
]

if (length(duplicated_definitions) > 0) {
  stop(
    "Duplicated cluster definitions found: ",
    paste(
      duplicated_definitions,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Add curation metadata
# ------------------------------------------------------------

state_match <- match(
  tcells_v3$V3_state_cluster,
  state_definitions$V3_state_cluster
)

excluded_match <- match(
  tcells_v3$V3_state_cluster,
  excluded_definitions$V3_state_cluster
)

tcells_v3$curation_status <- ifelse(
  !is.na(state_match),
  "Included",
  "Excluded"
)

tcells_v3$immune_state_id <-
  state_definitions$immune_state_id[
    state_match
  ]

tcells_v3$immune_state_name <-
  state_definitions$immune_state_name[
    state_match
  ]

tcells_v3$state_family <-
  state_definitions$state_family[
    state_match
  ]

tcells_v3$state_confidence <-
  state_definitions$confidence[
    state_match
  ]

tcells_v3$exclusion_reason <-
  excluded_definitions$exclusion_reason[
    excluded_match
  ]

# ------------------------------------------------------------
# Exclusion audit
# ------------------------------------------------------------

exclusion_audit <- tcells_v3@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  filter(
    curation_status == "Excluded"
  ) %>%
  count(
    V3_state_cluster,
    exclusion_reason,
    name = "n_cells"
  ) %>%
  mutate(
    percent_reference =
      100 * n_cells / ncol(tcells_v3)
  ) %>%
  arrange(
    as.numeric(V3_state_cluster)
  )

write.csv(
  exclusion_audit,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_excluded_cluster_audit.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Create curated reference object
# ------------------------------------------------------------

curated_cells <- colnames(tcells_v3)[
  tcells_v3$curation_status == "Included"
]

if (length(curated_cells) == 0) {
  stop(
    "No cells remained after curation."
  )
}

tcells_v3_curated <- subset(
  tcells_v3,
  cells = curated_cells
)

tcells_v3_curated$immune_state_id <- factor(
  tcells_v3_curated$immune_state_id,
  levels = state_definitions$immune_state_id
)

tcells_v3_curated$immune_state_name <- factor(
  tcells_v3_curated$immune_state_name,
  levels = state_definitions$immune_state_name
)

tcells_v3_curated$state_family <- factor(
  tcells_v3_curated$state_family,
  levels = unique(
    state_definitions$state_family
  )
)

tcells_v3_curated$state_confidence <- factor(
  tcells_v3_curated$state_confidence,
  levels = c(
    "High",
    "Medium",
    "Low"
  )
)

Idents(tcells_v3_curated) <-
  "immune_state_id"

# ------------------------------------------------------------
# Atlas statistics
# ------------------------------------------------------------

state_statistics <- tcells_v3_curated@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    ),
    immune_state_id = as.character(
      immune_state_id
    )
  ) %>%
  group_by(
    V3_state_cluster,
    immune_state_id
  ) %>%
  summarise(
    n_cells = n(),
    percent_curated_reference =
      100 * n() / ncol(tcells_v3_curated),
    n_samples = n_distinct(sample),
    n_donors = n_distinct(donor),
    n_cohorts = n_distinct(cohort),
    .groups = "drop"
  ) %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    ),
    immune_state_id = as.character(
      immune_state_id
    )
  )

# ------------------------------------------------------------
# Read and standardize marker summary
# ------------------------------------------------------------

marker_summary <- read.csv(
  marker_file,
  stringsAsFactors = FALSE
) %>%
  transmute(
    V3_state_cluster = as.character(
      cluster
    ),
    top_markers = as.character(
      top_markers
    )
  )

# ------------------------------------------------------------
# Read and standardize source-cluster crosswalk
# ------------------------------------------------------------

source_crosswalk <- read.csv(
  source_crosswalk_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    ),
    cluster_names = as.character(
      cluster_names
    )
  )

source_summary <- source_crosswalk %>%
  group_by(
    V3_state_cluster
  ) %>%
  arrange(
    desc(percent_within_V3_cluster),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 3
  ) %>%
  summarise(
    main_source_clusters = paste0(
      cluster_names,
      " (",
      round(
        percent_within_V3_cluster,
        1
      ),
      "%)",
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    )
  )

# ------------------------------------------------------------
# Read and standardize T09-associated program scores
# ------------------------------------------------------------

program_scores <- read.csv(
  program_score_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    )
  )

# ------------------------------------------------------------
# Build final curated atlas table
# ------------------------------------------------------------

curated_atlas_table <- state_definitions %>%
  mutate(
    V3_state_cluster = as.character(
      V3_state_cluster
    ),
    immune_state_id = as.character(
      immune_state_id
    )
  ) %>%
  left_join(
    state_statistics,
    by = c(
      "V3_state_cluster",
      "immune_state_id"
    )
  ) %>%
  left_join(
    marker_summary,
    by = "V3_state_cluster"
  ) %>%
  left_join(
    source_summary,
    by = "V3_state_cluster"
  ) %>%
  left_join(
    program_scores,
    by = "V3_state_cluster"
  ) %>%
  arrange(
    match(
      immune_state_id,
      state_definitions$immune_state_id
    )
  )

if (
  anyNA(
    curated_atlas_table$n_cells
  )
) {
  stop(
    "At least one curated state was not matched ",
    "to the state statistics."
  )
}

write.csv(
  curated_atlas_table,
  file.path(
    TAB_DIR,
    paste0(
      OUT_PREFIX,
      "_curated_immune_state_atlas.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Publication-oriented palettes
# ------------------------------------------------------------

state_palette <- c(
  "IVS-01" = "#4C78A8",
  "IVS-02" = "#72A0C1",
  "IVS-03" = "#59A14F",
  "IVS-04" = "#F28E2B",
  "IVS-05" = "#E15759",
  "IVS-06" = "#76B7B2",
  "IVS-07" = "#B07AA1",
  "IVS-08" = "#EDC948",
  "IVS-09" = "#AF7AA1",
  "IVS-10" = "#8C564B",
  "IVS-11" = "#FF7F0E"
)

family_palette <- c(
  "CD4 naive/memory" =
    "#4C78A8",
  "CD4 memory" =
    "#59A14F",
  "Activated CD4 memory" =
    "#F28E2B",
  "Cytotoxic CD8" =
    "#E15759",
  "CD8 naive/memory" =
    "#76B7B2",
  "CD8 effector memory" =
    "#B07AA1",
  "Innate-like T" =
    "#EDC948",
  "Regulatory T" =
    "#AF7AA1",
  "NK-like T" =
    "#8C564B",
  "Interferon-responsive T" =
    "#FF7F0E"
)

publication_theme <- theme(
  plot.title = element_text(
    face = "bold",
    size = 13
  ),
  legend.title = element_text(
    face = "bold"
  ),
  legend.text = element_text(
    size = 8
  )
)

# ------------------------------------------------------------
# UMAP by curated immune state
# ------------------------------------------------------------

p_curated_states <- DimPlot(
  tcells_v3_curated,
  reduction = "umap_harmony",
  group.by = "immune_state_id",
  cols = state_palette,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  labs(
    title = "Curated ImmunoVerse T-cell states",
    colour = "Immune state"
  ) +
  publication_theme

# ------------------------------------------------------------
# UMAP by higher-order state family
# ------------------------------------------------------------

p_curated_families <- DimPlot(
  tcells_v3_curated,
  reduction = "umap_harmony",
  group.by = "state_family",
  cols = family_palette,
  label = FALSE,
  raster = TRUE
) +
  labs(
    title = "Higher-order immune-state families",
    colour = "State family"
  ) +
  publication_theme

p_curated_atlas <- (
  p_curated_states +
    p_curated_families
) +
  plot_layout(
    widths = c(
      1,
      1.12
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  )

# ------------------------------------------------------------
# Save figures
# ------------------------------------------------------------

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_curated_immune_state_atlas_UMAP.pdf"
    )
  ),
  p_curated_atlas,
  width = 16,
  height = 7.5
)

ggsave(
  file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX,
      "_curated_immune_state_atlas_UMAP.png"
    )
  ),
  p_curated_atlas,
  width = 16,
  height = 7.5,
  dpi = 300
)

# ------------------------------------------------------------
# Save curated reference
# ------------------------------------------------------------

curated_rds <- file.path(
  PROCESSED_DIR,
  paste0(
    OUT_PREFIX,
    "_curated_atlas.rds"
  )
)

saveRDS(
  tcells_v3_curated,
  curated_rds,
  compress = FALSE
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

expected_cells <- sum(
  tcells_v3$curation_status == "Included"
)

if (
  ncol(tcells_v3_curated) !=
  expected_cells
) {
  stop(
    "Curated cell count does not match expectation."
  )
}

if (
  anyNA(
    tcells_v3_curated$immune_state_id
  )
) {
  stop(
    "Missing immune-state identifiers ",
    "remain in the curated object."
  )
}

if (
  n_distinct(
    tcells_v3_curated$immune_state_id
  ) !=
  nrow(state_definitions)
) {
  stop(
    "Unexpected number of immune states ",
    "in the curated reference."
  )
}

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

message("\nCurated atlas summary:")

curated_summary_console <- curated_atlas_table %>%
  select(
    immune_state_id,
    immune_state_name,
    state_family,
    confidence,
    n_cells,
    n_samples,
    n_donors
  )

print(
  as.data.frame(
    curated_summary_console
  ),
  row.names = FALSE
)

message("\nExcluded clusters:")

print(
  as.data.frame(
    exclusion_audit
  ),
  row.names = FALSE
)

message(
  "\nCurated object dimensions: ",
  paste(
    dim(tcells_v3_curated),
    collapse = " × "
  )
)

message(
  "Number of curated immune states: ",
  n_distinct(
    tcells_v3_curated$immune_state_id
  )
)

message(
  "Saved curated reference: ",
  curated_rds
)

message(
  "Curated V3 Immune State Atlas completed."
)