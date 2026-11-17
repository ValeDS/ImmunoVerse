# ============================================================
# 30_import_and_QC_GSE138266.R
# Import, cell calling, and per-sample QC for GSE138266
# ============================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "DropletUtils",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install Bioconductor packages with ",
    "BiocManager::install(\"DropletUtils\") and CRAN packages ",
    "with install.packages() before rerunning this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DropletUtils)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

GEO_ACCESSION_EXTERNAL <- "GSE138266"

AUDIT_DIR <- file.path(
  "results",
  "tables",
  "external_validation",
  GEO_ACCESSION_EXTERNAL
)
QC_TABLE_DIR <- file.path(
  AUDIT_DIR,
  "QC"
)
QC_FIGURE_DIR <- file.path(
  "results",
  "figures",
  "external_validation",
  GEO_ACCESSION_EXTERNAL
)
PROCESSED_EXTERNAL_DIR <- file.path(
  "data",
  "processed",
  "external_validation",
  GEO_ACCESSION_EXTERNAL
)
EMPTYDROPS_DIR <- file.path(
  PROCESSED_EXTERNAL_DIR,
  "emptyDrops"
)

dir.create(
  QC_TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  QC_FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  PROCESSED_EXTERNAL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  EMPTYDROPS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

SAMPLE_MANIFEST_FILE <- file.path(
  AUDIT_DIR,
  paste0(
    GEO_ACCESSION_EXTERNAL,
    "_sample_manifest_provisional.csv"
  )
)
FILE_MANIFEST_FILE <- file.path(
  AUDIT_DIR,
  paste0(
    GEO_ACCESSION_EXTERNAL,
    "_processed_file_manifest.csv"
  )
)

OUTPUT_OBJECT_FILE <- file.path(
  PROCESSED_EXTERNAL_DIR,
  paste0(
    GEO_ACCESSION_EXTERNAL,
    "_all_leukocytes_QC_counts.rds"
  )
)

required_input_files <- c(
  SAMPLE_MANIFEST_FILE,
  FILE_MANIFEST_FILE
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Required audit file(s) not found:\n",
    paste(missing_input_files, collapse = "\n"),
    "\nRun scripts/29_download_and_audit_GSE138266.R first."
  )
}

# ------------------------------------------------------------
# QC thresholds
# These match the discovery-workflow thresholds where possible.
# ------------------------------------------------------------

MIN_FEATURES <- 200L
MAX_FEATURES <- 6000L
MIN_COUNTS <- 500L
MAX_COUNTS <- 50000L
MAX_MITO_PERCENT <- 15

# A matrix with this many barcodes is treated as an unfiltered
# droplet matrix and subjected to emptyDrops cell calling.
RAW_DROPLET_BARCODE_THRESHOLD <- 100000L
EMPTYDROPS_FDR <- 0.01

# ------------------------------------------------------------
# Load and validate manifests
# ------------------------------------------------------------

sample_manifest <- read.csv(
  SAMPLE_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  transmute(
    geo_accession = as.character(geo_accession),
    sample = as.character(sample),
    donor = as.character(donor),
    compartment = as.character(compartment),
    group = as.character(provisional_group),
    characteristics = as.character(characteristics)
  )

file_manifest <- read.csv(
  FILE_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    geo_accession = as.character(geo_accession),
    file_role = as.character(file_role),
    file_path = as.character(file_path)
  )

if (
  nrow(sample_manifest) != 22L ||
    anyDuplicated(sample_manifest$geo_accession) ||
    anyDuplicated(sample_manifest$sample)
) {
  stop(
    "Expected 22 uniquely identified samples in the manifest."
  )
}

required_roles <- c(
  "matrix",
  "barcodes",
  "features"
)

input_triplets <- file_manifest %>%
  filter(
    file_role %in% required_roles,
    !is.na(geo_accession)
  ) %>%
  select(
    geo_accession,
    file_role,
    file_path
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = file_role,
    values_from = file_path
  ) %>%
  left_join(
    sample_manifest,
    by = "geo_accession"
  ) %>%
  arrange(
    factor(
      compartment,
      levels = c("PBMC", "CSF")
    ),
    group,
    donor
  )

if (
  nrow(input_triplets) != 22L ||
    anyNA(input_triplets[
      ,
      required_roles,
      drop = FALSE
    ])
) {
  stop(
    "Unable to reconstruct 22 complete matrix/barcode/feature ",
    "triplets."
  )
}

all_matrix_files <- unlist(
  input_triplets[
    ,
    required_roles,
    drop = FALSE
  ],
  use.names = FALSE
)

if (!all(file.exists(all_matrix_files))) {
  stop(
    "One or more processed matrix files are missing."
  )
}

# ------------------------------------------------------------
# Per-sample import, optional emptyDrops, and QC
# ------------------------------------------------------------

sample_objects <- vector(
  mode = "list",
  length = nrow(input_triplets)
)
names(sample_objects) <- input_triplets$sample

qc_summary_list <- vector(
  mode = "list",
  length = nrow(input_triplets)
)

qc_quantile_list <- vector(
  mode = "list",
  length = nrow(input_triplets)
)

for (i in seq_len(nrow(input_triplets))) {

  sample_info <- input_triplets[i, , drop = FALSE]
  sample_id <- sample_info$sample[[1]]
  geo_id <- sample_info$geo_accession[[1]]

  message(
    "[",
    i,
    "/",
    nrow(input_triplets),
    "] Importing ",
    sample_id,
    " (",
    geo_id,
    ")..."
  )

  counts <- Seurat::ReadMtx(
    mtx = sample_info$matrix[[1]],
    cells = sample_info$barcodes[[1]],
    features = sample_info$features[[1]],
    cell.column = 1L,
    feature.column = 2L,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  if (
    !inherits(counts, "dgCMatrix") ||
      nrow(counts) == 0L ||
      ncol(counts) == 0L
  ) {
    stop(
      "The imported count matrix is empty or invalid for ",
      sample_id,
      "."
    )
  }

  raw_barcodes <- ncol(counts)
  raw_droplet_matrix <-
    raw_barcodes >= RAW_DROPLET_BARCODE_THRESHOLD

  cell_call_method <- if (
    raw_droplet_matrix
  ) {
    "emptyDrops_FDR_0.01"
  } else {
    "GEO_supplied_filtered_matrix"
  }

  if (raw_droplet_matrix) {

    message(
      "  ",
      format(
        raw_barcodes,
        big.mark = ",",
        scientific = FALSE
      ),
      " barcodes detected; running emptyDrops..."
    )

    emptydrops_file <- file.path(
      EMPTYDROPS_DIR,
      paste0(
        geo_id,
        "_",
        sample_id,
        "_emptyDrops.rds"
      )
    )

    if (file.exists(emptydrops_file)) {
      emptydrops_result <- readRDS(
        emptydrops_file
      )

      if (
        nrow(emptydrops_result) != ncol(counts) ||
          !identical(
            rownames(emptydrops_result),
            colnames(counts)
          )
      ) {
        stop(
          "The cached emptyDrops result does not match ",
          sample_id,
          ". Move it before rerunning."
        )
      }

      message(
        "  Validated cached emptyDrops result loaded."
      )
    } else {
      emptydrops_result <- DropletUtils::emptyDrops(
        counts
      )

      saveRDS(
        emptydrops_result,
        emptydrops_file,
        compress = TRUE
      )
    }

    cell_called <- which(
      !is.na(emptydrops_result$FDR) &
        emptydrops_result$FDR <= EMPTYDROPS_FDR
    )

    if (length(cell_called) < 100L) {
      stop(
        "emptyDrops identified fewer than 100 cells in ",
        sample_id,
        ". Inspect the cell-calling output."
      )
    }
  } else {
    cell_called <- seq_len(ncol(counts))
  }

  counts_called <- counts[
    ,
    cell_called,
    drop = FALSE
  ]

  n_count <- Matrix::colSums(
    counts_called
  )
  n_feature <- Matrix::colSums(
    counts_called > 0
  )

  mitochondrial_features <- grep(
    "^MT-",
    rownames(counts_called),
    value = TRUE,
    ignore.case = FALSE
  )

  if (length(mitochondrial_features) > 0L) {
    percent_mt <- 100 *
      Matrix::colSums(
        counts_called[
          mitochondrial_features,
          ,
          drop = FALSE
        ]
      ) /
      pmax(n_count, 1)
  } else {
    percent_mt <- rep(
      0,
      ncol(counts_called)
    )
  }

  passes_qc <-
    n_feature >= MIN_FEATURES &
    n_feature <= MAX_FEATURES &
    n_count >= MIN_COUNTS &
    n_count <= MAX_COUNTS &
    percent_mt <= MAX_MITO_PERCENT

  retained_cells <- which(passes_qc)

  if (length(retained_cells) < 50L) {
    stop(
      "Fewer than 50 cells passed QC in ",
      sample_id,
      ". Inspect the per-sample QC metrics."
    )
  }

  counts_qc <- counts_called[
    ,
    retained_cells,
    drop = FALSE
  ]

  original_barcodes <- colnames(counts_qc)
  unique_cell_ids <- paste(
    sample_id,
    original_barcodes,
    sep = "__"
  )
  colnames(counts_qc) <- unique_cell_ids

  cell_metadata <- data.frame(
    row.names = unique_cell_ids,
    cell_id = unique_cell_ids,
    original_barcode = original_barcodes,
    geo_accession = geo_id,
    sample = sample_id,
    donor = sample_info$donor[[1]],
    group = sample_info$group[[1]],
    compartment = sample_info$compartment[[1]],
    dataset = GEO_ACCESSION_EXTERNAL,
    chemistry = "10x_3prime_v2",
    cell_call_method = cell_call_method,
    nCount_RNA = as.numeric(
      n_count[retained_cells]
    ),
    nFeature_RNA = as.numeric(
      n_feature[retained_cells]
    ),
    percent.mt = as.numeric(
      percent_mt[retained_cells]
    ),
    stringsAsFactors = FALSE
  )

  sample_object <- Seurat::CreateSeuratObject(
    counts = counts_qc,
    assay = "RNA",
    project = GEO_ACCESSION_EXTERNAL,
    meta.data = cell_metadata,
    min.cells = 0L,
    min.features = 0L
  )

  sample_objects[[sample_id]] <- sample_object

  called_cells <- ncol(counts_called)
  retained_qc_cells <- ncol(counts_qc)

  qc_summary_list[[i]] <- data.frame(
    geo_accession = geo_id,
    sample = sample_id,
    donor = sample_info$donor[[1]],
    group = sample_info$group[[1]],
    compartment = sample_info$compartment[[1]],
    raw_barcodes = raw_barcodes,
    raw_droplet_matrix = raw_droplet_matrix,
    cell_call_method = cell_call_method,
    cell_called_barcodes = called_cells,
    retained_QC_cells = retained_qc_cells,
    percent_retained_after_cell_call =
      100 * retained_qc_cells / called_cells,
    median_nCount_RNA = median(
      n_count[retained_cells]
    ),
    median_nFeature_RNA = median(
      n_feature[retained_cells]
    ),
    median_percent_mt = median(
      percent_mt[retained_cells]
    ),
    stringsAsFactors = FALSE
  )

  quantile_probabilities <- c(
    0,
    0.01,
    0.25,
    0.50,
    0.75,
    0.99,
    1
  )

  qc_quantile_list[[i]] <- bind_rows(
    data.frame(
      sample = sample_id,
      metric = "nCount_RNA",
      probability = quantile_probabilities,
      value = as.numeric(
        quantile(
          n_count,
          probs = quantile_probabilities,
          na.rm = TRUE
        )
      )
    ),
    data.frame(
      sample = sample_id,
      metric = "nFeature_RNA",
      probability = quantile_probabilities,
      value = as.numeric(
        quantile(
          n_feature,
          probs = quantile_probabilities,
          na.rm = TRUE
        )
      )
    ),
    data.frame(
      sample = sample_id,
      metric = "percent.mt",
      probability = quantile_probabilities,
      value = as.numeric(
        quantile(
          percent_mt,
          probs = quantile_probabilities,
          na.rm = TRUE
        )
      )
    )
  )

  rm(
    counts,
    counts_called,
    counts_qc,
    sample_object
  )

  if (exists("emptydrops_result")) {
    rm(emptydrops_result)
  }

  gc(verbose = FALSE)

  message(
    "  Retained ",
    format(
      retained_qc_cells,
      big.mark = ",",
      scientific = FALSE
    ),
    " cells after QC."
  )
}

# ------------------------------------------------------------
# Merge the 22 QC-filtered objects
# ------------------------------------------------------------

if (
  length(sample_objects) != 22L ||
    any(vapply(sample_objects, is.null, logical(1)))
) {
  stop("One or more sample objects are missing.")
}

message(
  "Merging 22 QC-filtered sample objects..."
)

external_qc <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  merge.data = FALSE,
  merge.dr = FALSE
)

DefaultAssay(external_qc) <- "RNA"

if (
  anyDuplicated(colnames(external_qc)) ||
    n_distinct(external_qc$sample) != 22L ||
    n_distinct(external_qc$donor) != 12L ||
    anyNA(external_qc$group) ||
    anyNA(external_qc$compartment)
) {
  stop(
    "The merged external-validation object failed metadata ",
    "validation."
  )
}

# Join the per-sample count layers into one counts layer.
external_qc <- SeuratObject::JoinLayers(
  external_qc,
  assay = "RNA"
)

if (
  !identical(
    SeuratObject::Layers(external_qc[["RNA"]]),
    "counts"
  )
) {
  stop(
    "Expected one counts layer after joining sample layers."
  )
}

saveRDS(
  external_qc,
  OUTPUT_OBJECT_FILE,
  compress = FALSE
)

# ------------------------------------------------------------
# Save QC tables
# ------------------------------------------------------------

qc_summary <- bind_rows(
  qc_summary_list
) %>%
  arrange(
    factor(
      compartment,
      levels = c("PBMC", "CSF")
    ),
    group,
    donor
  )

qc_quantiles <- bind_rows(
  qc_quantile_list
)

qc_thresholds <- data.frame(
  parameter = c(
    "MIN_FEATURES",
    "MAX_FEATURES",
    "MIN_COUNTS",
    "MAX_COUNTS",
    "MAX_MITO_PERCENT",
    "RAW_DROPLET_BARCODE_THRESHOLD",
    "EMPTYDROPS_FDR"
  ),
  value = c(
    MIN_FEATURES,
    MAX_FEATURES,
    MIN_COUNTS,
    MAX_COUNTS,
    MAX_MITO_PERCENT,
    RAW_DROPLET_BARCODE_THRESHOLD,
    EMPTYDROPS_FDR
  )
)

merged_object_summary <- external_qc[[]] %>%
  count(
    group,
    compartment,
    name = "retained_cells"
  ) %>%
  arrange(
    group,
    compartment
  )

write.csv(
  qc_summary,
  file.path(
    QC_TABLE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_per_sample_QC_summary.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  qc_quantiles,
  file.path(
    QC_TABLE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_pre_filter_QC_quantiles.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  qc_thresholds,
  file.path(
    QC_TABLE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_QC_thresholds.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  merged_object_summary,
  file.path(
    QC_TABLE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_post_QC_study_design.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# QC figure
# ------------------------------------------------------------

qc_long <- qc_summary %>%
  select(
    sample,
    group,
    compartment,
    raw_barcodes,
    cell_called_barcodes,
    retained_QC_cells
  ) %>%
  pivot_longer(
    cols = c(
      raw_barcodes,
      cell_called_barcodes,
      retained_QC_cells
    ),
    names_to = "stage",
    values_to = "cells"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "raw_barcodes",
        "cell_called_barcodes",
        "retained_QC_cells"
      ),
      labels = c(
        "Raw barcodes",
        "Cell-called",
        "Post-QC"
      )
    )
  )

compartment_colors <- c(
  "PBMC" = "#4C78A8",
  "CSF" = "#E45756"
)

panel_A <- ggplot(
  qc_long,
  aes(
    x = sample,
    y = cells,
    color = compartment,
    group = interaction(sample, compartment)
  )
) +
  geom_line(
    aes(group = sample),
    color = "grey75",
    linewidth = 0.35
  ) +
  geom_point(size = 1.8) +
  facet_grid(
    stage ~ .,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = compartment_colors
  ) +
  scale_y_continuous(
    labels = scales::label_comma()
  ) +
  labs(
    title = "GSE138266 cell recovery across import stages",
    subtitle = paste0(
      "PST83775_CSF is cell-called with emptyDrops; ",
      "all samples then undergo uniform QC"
    ),
    x = NULL,
    y = "Cells or barcodes",
    color = "Compartment"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 6
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

panel_B <- qc_summary %>%
  select(
    sample,
    group,
    compartment,
    median_nCount_RNA,
    median_nFeature_RNA,
    median_percent_mt
  ) %>%
  pivot_longer(
    cols = starts_with("median_"),
    names_to = "metric",
    values_to = "median_value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "median_nCount_RNA",
        "median_nFeature_RNA",
        "median_percent_mt"
      ),
      labels = c(
        "Median UMI counts",
        "Median detected genes",
        "Median mitochondrial %"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = interaction(
        compartment,
        group,
        sep = " / "
      ),
      y = median_value,
      color = compartment
    )
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    color = "grey45"
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1
  ) +
  scale_color_manual(
    values = compartment_colors
  ) +
  labs(
    title = "Post-QC sample characteristics",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

qc_figure <- panel_A + panel_B +
  plot_layout(
    widths = c(1.8, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  filename = file.path(
    QC_FIGURE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_import_and_QC.png"
    )
  ),
  plot = qc_figure,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    QC_FIGURE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_import_and_QC.pdf"
    )
  ),
  plot = qc_figure,
  width = 14,
  height = 10,
  device = grDevices::pdf,
  bg = "white"
)

# ------------------------------------------------------------
# Final validation report
# ------------------------------------------------------------

validation_summary <- data.frame(
  metric = c(
    "input_samples",
    "input_donors",
    "raw_droplet_matrices",
    "cell_called_barcodes",
    "post_QC_cells",
    "post_QC_PBMC_cells",
    "post_QC_CSF_cells",
    "genes",
    "RNA_layers"
  ),
  value = c(
    nrow(qc_summary),
    n_distinct(qc_summary$donor),
    sum(qc_summary$raw_droplet_matrix),
    sum(qc_summary$cell_called_barcodes),
    ncol(external_qc),
    sum(external_qc$compartment == "PBMC"),
    sum(external_qc$compartment == "CSF"),
    nrow(external_qc),
    length(
      SeuratObject::Layers(external_qc[["RNA"]])
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  validation_summary,
  file.path(
    QC_TABLE_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_import_QC_validation_summary.csv"
    )
  ),
  row.names = FALSE
)

message("")
message("GSE138266 import and QC completed.")
message("")
message("Validation summary:")
print(
  validation_summary,
  row.names = FALSE
)
message("")
message("Post-QC study design:")
print(
  as.data.frame(merged_object_summary),
  row.names = FALSE
)
message("")
message("Per-sample QC summary:")
print(
  as.data.frame(qc_summary),
  row.names = FALSE
)
message("")
message(
  "Saved counts-only Seurat object: ",
  OUTPUT_OBJECT_FILE
)
