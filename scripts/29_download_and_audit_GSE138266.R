# ============================================================
# 29_download_and_audit_GSE138266.R
# External-validation dataset acquisition and structural audit
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

GEO_ACCESSION_EXTERNAL <- "GSE138266"

RAW_BASE_DIR <- "data/raw"
RAW_DATASET_DIR <- file.path(
  RAW_BASE_DIR,
  GEO_ACCESSION_EXTERNAL
)
EXTRACTED_DIR <- file.path(
  RAW_DATASET_DIR,
  "extracted"
)
AUDIT_DIR <- file.path(
  "results",
  "tables",
  "external_validation",
  GEO_ACCESSION_EXTERNAL
)

dir.create(
  RAW_BASE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  RAW_DATASET_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  EXTRACTED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Download and save GEO sample metadata
# ------------------------------------------------------------

message(
  "Downloading GEO metadata for ",
  GEO_ACCESSION_EXTERNAL,
  "..."
)

gse_objects <- GEOquery::getGEO(
  GEO_ACCESSION_EXTERNAL,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse_objects) < 1L) {
  stop("No GEO series-matrix object was returned.")
}

pheno <- Biobase::pData(gse_objects[[1]]) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
  rownames_to_column("geo_rowname")

required_pheno_columns <- c(
  "geo_accession",
  "title"
)

missing_pheno_columns <- setdiff(
  required_pheno_columns,
  colnames(pheno)
)

if (length(missing_pheno_columns) > 0L) {
  stop(
    "The GEO phenotype table is missing column(s): ",
    paste(missing_pheno_columns, collapse = ", ")
  )
}

write.csv(
  pheno,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_GEO_sample_metadata_raw.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Construct a normalized provisional sample manifest
# ------------------------------------------------------------

characteristic_columns <- grep(
  "^characteristics",
  colnames(pheno),
  value = TRUE,
  ignore.case = TRUE
)

if (length(characteristic_columns) > 0L) {
  characteristic_text <- apply(
    pheno[
      ,
      characteristic_columns,
      drop = FALSE
    ],
    1L,
    function(x) {
      paste(
        x[
          !is.na(x) &
            nzchar(as.character(x))
        ],
        collapse = " | "
      )
    }
  )
} else {
  characteristic_text <- rep(
    "",
    nrow(pheno)
  )
}

sample_manifest <- pheno %>%
  transmute(
    geo_accession = as.character(geo_accession),
    sample = as.character(title),
    donor = sub(
      "_(CSF|PBMC)$",
      "",
      as.character(title),
      ignore.case = TRUE
    ),
    compartment = case_when(
      grepl(
        "_CSF$",
        as.character(title),
        ignore.case = TRUE
      ) ~ "CSF",
      grepl(
        "_PBMC$",
        as.character(title),
        ignore.case = TRUE
      ) ~ "PBMC",
      TRUE ~ NA_character_
    ),
    provisional_group = case_when(
      grepl(
        "^MS",
        as.character(title),
        ignore.case = TRUE
      ) ~ "MS",
      grepl(
        "^(PTC|PST)",
        as.character(title),
        ignore.case = TRUE
      ) ~ "Control",
      TRUE ~ "Unresolved"
    ),
    characteristics = characteristic_text
  ) %>%
  arrange(
    donor,
    factor(
      compartment,
      levels = c("PBMC", "CSF")
    )
  )

if (
  nrow(sample_manifest) != 22L ||
    anyDuplicated(sample_manifest$geo_accession) ||
    anyNA(sample_manifest$compartment)
) {
  stop(
    "Unexpected GSE138266 metadata structure: expected ",
    "22 uniquely identified CSF/PBMC samples."
  )
}

write.csv(
  sample_manifest,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_sample_manifest_provisional.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Download the processed count-matrix archive
# ------------------------------------------------------------

archive_candidates <- list.files(
  RAW_DATASET_DIR,
  pattern = "RAW\\.tar$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archive_candidates) == 0L) {
  message(
    "Downloading GEO supplementary files. ",
    "The processed archive is approximately 252 MB..."
  )

  GEOquery::getGEOSuppFiles(
    GEO_ACCESSION_EXTERNAL,
    baseDir = RAW_BASE_DIR,
    makeDirectory = TRUE,
    fetch_files = TRUE
  )

  archive_candidates <- list.files(
    RAW_DATASET_DIR,
    pattern = "RAW\\.tar$",
    full.names = TRUE,
    ignore.case = TRUE
  )
} else {
  message(
    "Existing supplementary archive found; ",
    "download skipped."
  )
}

if (length(archive_candidates) != 1L) {
  stop(
    "Expected one GSE138266 RAW tar archive but found ",
    length(archive_candidates),
    "."
  )
}

archive_file <- archive_candidates[[1]]

# ------------------------------------------------------------
# Extract the MTX/TSV files without decompressing their .gz files
# ------------------------------------------------------------

existing_extracted_files <- list.files(
  EXTRACTED_DIR,
  recursive = TRUE,
  full.names = TRUE
)

if (length(existing_extracted_files) == 0L) {
  message("Extracting processed MTX/TSV files...")

  utils::untar(
    archive_file,
    exdir = EXTRACTED_DIR
  )
} else {
  message(
    "Extracted files already exist; extraction skipped."
  )
}

processed_files <- list.files(
  EXTRACTED_DIR,
  recursive = TRUE,
  full.names = TRUE
)

if (length(processed_files) == 0L) {
  stop("No files were found after archive extraction.")
}

file_basenames <- basename(processed_files)

file_manifest <- data.frame(
  file_path = processed_files,
  file_name = file_basenames,
  size_MB = round(
    file.info(processed_files)$size / 1024^2,
    3
  ),
  file_role = case_when(
    grepl(
      "matrix\\.mtx(\\.gz)?$",
      file_basenames,
      ignore.case = TRUE
    ) ~ "matrix",
    grepl(
      "barcodes\\.tsv(\\.gz)?$",
      file_basenames,
      ignore.case = TRUE
    ) ~ "barcodes",
    grepl(
      "(genes|features)\\.tsv(\\.gz)?$",
      file_basenames,
      ignore.case = TRUE
    ) ~ "features",
    TRUE ~ "other"
  ),
  geo_accession = sub(
    ".*(GSM[0-9]+).*",
    "\\1",
    file_basenames
  ),
  stringsAsFactors = FALSE
)

file_manifest$geo_accession[
  !grepl(
    "GSM[0-9]+",
    file_manifest$file_name
  )
] <- NA_character_

write.csv(
  file_manifest,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_processed_file_manifest.csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Audit one matrix/barcode/feature triplet per GEO sample
# ------------------------------------------------------------

triplet_audit <- file_manifest %>%
  filter(
    file_role %in%
      c("matrix", "barcodes", "features"),
    !is.na(geo_accession)
  ) %>%
  count(
    geo_accession,
    file_role,
    name = "n_files"
  ) %>%
  tidyr::complete(
    geo_accession =
      sample_manifest$geo_accession,
    file_role =
      c("matrix", "barcodes", "features"),
    fill = list(n_files = 0L)
  ) %>%
  tidyr::pivot_wider(
    names_from = file_role,
    values_from = n_files,
    names_prefix = "n_"
  ) %>%
  left_join(
    sample_manifest,
    by = "geo_accession"
  ) %>%
  mutate(
    complete_10x_triplet =
      n_matrix == 1L &
      n_barcodes == 1L &
      n_features == 1L
  ) %>%
  arrange(
    donor,
    factor(
      compartment,
      levels = c("PBMC", "CSF")
    )
  )

write.csv(
  triplet_audit,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_10x_triplet_audit.csv"
    )
  ),
  row.names = FALSE
)

if (
  nrow(triplet_audit) != 22L ||
    !all(triplet_audit$complete_10x_triplet)
) {
  stop(
    "The processed-file audit did not recover one complete ",
    "10x matrix/barcode/feature triplet for every sample. ",
    "Inspect the saved file manifest before continuing."
  )
}

# ------------------------------------------------------------
# Pairing and study-design audit
# ------------------------------------------------------------

donor_pairing <- sample_manifest %>%
  count(
    donor,
    provisional_group,
    compartment,
    name = "n_samples"
  ) %>%
  tidyr::pivot_wider(
    names_from = compartment,
    values_from = n_samples,
    values_fill = 0L
  ) %>%
  mutate(
    has_paired_PBMC_CSF =
      PBMC == 1L &
      CSF == 1L
  ) %>%
  arrange(
    provisional_group,
    donor
  )

write.csv(
  donor_pairing,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_donor_pairing_audit.csv"
    )
  ),
  row.names = FALSE
)

study_design_summary <- sample_manifest %>%
  count(
    provisional_group,
    compartment,
    name = "n_samples"
  ) %>%
  arrange(
    provisional_group,
    compartment
  )

pairing_summary <- donor_pairing %>%
  summarise(
    total_donors = n(),
    paired_donors = sum(has_paired_PBMC_CSF),
    CSF_only_donors = sum(CSF == 1L & PBMC == 0L),
    PBMC_only_donors = sum(PBMC == 1L & CSF == 0L)
  )

write.csv(
  study_design_summary,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_study_design_summary.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  pairing_summary,
  file.path(
    AUDIT_DIR,
    paste0(
      GEO_ACCESSION_EXTERNAL,
      "_pairing_summary.csv"
    )
  ),
  row.names = FALSE
)

message("")
message("GSE138266 acquisition and audit completed.")
message("")
message("Study-design summary:")
print(
  as.data.frame(study_design_summary),
  row.names = FALSE
)
message("")
message("Pairing summary:")
print(
  as.data.frame(pairing_summary),
  row.names = FALSE
)
message("")
message("Processed-file audit:")
print(
  data.frame(
    total_files = nrow(file_manifest),
    complete_10x_samples = sum(
      triplet_audit$complete_10x_triplet
    ),
    incomplete_10x_samples = sum(
      !triplet_audit$complete_10x_triplet
    )
  ),
  row.names = FALSE
)
message("")
message(
  "Audit tables saved under: ",
  AUDIT_DIR
)
